/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2022 The EXP-T developers.
 *
 *  This file is part of EXP-T.
 *
 *  EXP-T is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EXP-T is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with EXP-T.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  E-mail:        exp-t-program@googlegroups.com
 *  Google Groups: https://groups.google.com/d/forum/exp-t-program
 */

/*
 * Subroutines for reordering of multidimensional tensors.
 * Tensors are stored in a linear array row-wise (C-like style).
 *
 * Algorithm notes:
 * 1. Tensors are splitted into symmetry blocks (blocks) by irrep. Every
 *    symblock "knows" its dimensions and ranges of spinor indices for each
 *    dimension.
 * 2. Using metainfo you can easily calculate sizes of tensor for each dimension
 * 3. Tensor elements are stored linearly, row-wise. Sequential index of the
 *    element can be transformed to the "compound index" in the multidimensional
 *    array. Algorithm for this transformation is a straightforward generaliza-
 *    tion of indexing of 2D matrix elements.
 *    Example 1. 2-dim matrix MxN. For element with "compound" index (i,j) its
 *      "linear" index is equal to M*(i-1)+j (-1 since array indices start with
 *      1 in Fortran)
 *    Example 2. 4-dimensional tensor with dimensions [1:L,1:M,1:N,1:K].
 *      For element (l,i,j,k) linear index is equal to
 *      (((M*(l-1)+(i-1)))*N+(j-1))*K+k ==
 *      == MNK*(l-1) + NK*(i-1) + K*(j-1) + (k-1) + 1
 *      (+1 since we use this address to index Fortran array)
 *      The latter expression can be rewritten as a scalar product. Multipliers
 *      MNK,NK,etc are calculated only once at start.
 * 4. Algorithm from pt 3 can be reversed to calculate compound index from linear
 * 5. This algorithm represents a very simple loop over tensor elements without
 *    any data dependence => simply can be parallelized! However, parallelisation
 *    inside of symblock will be efficient only for extremely large blocks
 *    (like those in triples and quadruples)
 * 6. Asymptotic behaviour: O(N^{2n}) for the n-particle diagram.
 *
 * 2017-2021 Alexander Oleynichenko
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include "omp.h"

#include "comdef.h"
#include "cuda_code.h"
#include "engine.h"
#include "datamodel.h"
#include "templates.h"
#include "options.h"

void reverse_perm(int n, int *a, int *ainv);


/**
 * Reorders diagram 'src'; puts result into 'target'
 *
 * Arguments:
 *   src        symbolic name of the diagram to be reordered
 *   target     symbolic name of the new diagram -- result
 *   perm_str   permutation, smth like "124365" (reorders dim-s 3<->4 and 5<->6)
 */
void reorder(char *src, char *target, char *perm_str)
{
    diagram_t *dg_src, *dg_tgt;
    int i;
    int perm[32];

    timer_new_entry("reorder", "Multidim transposition (reorder)");
    timer_start("reorder");

    dg_src = diagram_stack_find(src);
    if (dg_src == NULL) {
        errquit("reorder(): diagram '%s' not found", src);
    }

    // parse "permutation string"
    for (i = 0; i < dg_src->rank; i++) {
        perm[i] = perm_str[i] - '0';
    }

    // reorder
    dg_tgt = diagram_reorder(dg_src, perm, dg_src->only_unique /*(target[0] == '$') ? 1 : 0*/);
    strcpy(dg_tgt->name, target);

    // save new diagram to stack
    // (replace the old diagram named 'target' if needed)
    if (diagram_stack_find(target) != NULL) {
        diagram_stack_replace(target, dg_tgt);
    }
    else {
        diagram_stack_push(dg_tgt);
    }

    timer_stop("reorder");
}

/* include "template function" reorder_block_<TYPENAME>
 * for real and complex numbers
 */
#define TYPENAME double_complex_t
#include "reorder_block.c"
#undef TYPENAME

#define TYPENAME double
#include "reorder_block.c"
#undef TYPENAME


diagram_t *diagram_reorder(diagram_t *dg, int *perm, int perm_unique)
{
    int rank;
    diagram_t *tgt;
    char qparts[CC_DIAGRAM_MAX_RANK];
    char valence[CC_DIAGRAM_MAX_RANK];
    char t3space[CC_DIAGRAM_MAX_RANK];
    char order[CC_DIAGRAM_MAX_RANK];

    rank = dg->rank;

    // dirty code!
    // create strings from in arrays
    for (int i = 0; i < rank; i++) {
        perm[i] = perm[i] - 1;
        qparts[i] = dg->qparts[perm[i]];
        valence[i] = dg->valence[perm[i]] + '0';
        t3space[i] = dg->t3space[perm[i]] + '0';
        order[i] = dg->order[perm[i]] + '0';
    }
    qparts[rank] = '\0';
    valence[rank] = '\0';
    t3space[rank] = '\0';
    order[rank] = '\0';

    tgt = diagram_new("rdr " /*dg->name*/, qparts, valence, t3space, order, perm_unique);

    // name of the new diagram consist of the old one + "_copy_" + clone's ID
    sprintf(tgt->name, "%s_rdr_%ld", dg->name, tgt->dg_id);

    // почему нет перестановки order, valence, qparts для блока?
    // эти характеристики для блоков вообще используются?
    // copy data: block-by-block
    for (size_t ib1 = 0; ib1 < dg->n_blocks; ib1++) {
        block_t *sb_src = dg->blocks[ib1];

        int src_spib_rdr[CC_DIAGRAM_MAX_RANK];
        for (int idim = 0; idim < rank; idim++) {
            src_spib_rdr[idim] = sb_src->spinor_blocks[perm[idim]];
        }

        block_t *sb_tgt = diagram_get_block(tgt, src_spib_rdr);

        if (sb_src->is_unique != sb_tgt->is_unique) {
            printf("reorder(): %d %d\n", sb_src->is_unique, sb_tgt->is_unique);
            exit(0);
        }

        if (sb_src->is_unique == 0) {
            continue;
        }

        if (arith == CC_ARITH_COMPLEX) {
            TEMPLATE(reorder_block, double_complex_t)(sb_src, sb_tgt, perm);
        }
        else {
            TEMPLATE(reorder_block, double)(sb_src, sb_tgt, perm);
        }

    }  // end loop over blocks

    return tgt;
}


/* helper functions */

typedef struct {
    int i;
    int j;
} pair_t;


static int pair_t_less(pair_t *a, pair_t *b)
{
    return a->j - b->j;
}


void reverse_perm(int n, int *a, int *ainv)
{
    pair_t p[CC_DIAGRAM_MAX_RANK];

    for (int i = 0; i < n; i++) {
        p[i].i = i;
        p[i].j = a[i];
    }
    qsort(p, n, sizeof(pair_t), (int (*)(const void *, const void *)) pair_t_less);
    for (int i = 0; i < n; i++) {
        ainv[i] = p[i].i;
    }
}


/**
 * restores permutationally non-unique block by reordering of its unique counterpart
 * @param dg diagram to which the block belongs
 * @param b  block to be restored
 */
void restore_block(diagram_t *dg, block_t *b)
{
    int uniq_spinor_blocks[CC_DIAGRAM_MAX_RANK];

    void transform(int n, int *idx, int *out, int *perm, int shift);

    transform(b->rank, b->spinor_blocks, uniq_spinor_blocks, b->perm_to_unique, 0);

    block_t *uniq_block = diagram_get_block(dg, uniq_spinor_blocks);

    if (b->storage_type == CC_DIAGRAM_IN_MEM) {
        printf("restore_block(): already in mem!\n");
        return;
    }

    b->storage_type = CC_DIAGRAM_IN_MEM;
    b->buf = (double complex *) cc_calloc(b->size, SIZEOF_WORKING_TYPE);
    block_store(b);

    if (arith == CC_ARITH_COMPLEX) {
        TEMPLATE(reorder_block, double_complex_t)(uniq_block, b, b->perm_from_unique);

        block_load(b);
        double complex *zbuf = b->buf;
        for (size_t i = 0; i < b->size; i++) {
            zbuf[i] *= b->sign;
        }
        block_store(b);
    }
    else {
        TEMPLATE(reorder_block, double)(uniq_block, b, b->perm_from_unique);

        block_load(b);
        double *dbuf = (double *) b->buf;
        for (size_t i = 0; i < b->size; i++) {
            dbuf[i] *= b->sign;
        }
        block_store(b);
    }
}
