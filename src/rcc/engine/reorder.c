/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2025 The EXP-T developers.
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
 * 1. Tensors are split into symmetry blocks (blocks) by irrep. Every
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
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include "omp.h"

#include "comdef.h"
#include "cuda_code.h"
#include "engine.h"
#include "c_templates.h"
#include "options.h"
#include "utils.h"
#include "linalg.h"


void reverse_perm(int n, const int *direct_perm, int *inv_perm);

void reorder_block(block_t *source_block, block_t *target_block, int *perm);


/**
 * Reorders diagram 'src_diargam_name'; puts result into 'target_diagram_name'
 * (tensor transposition)
 *
 * Arguments:
 *   src_diargam_name        symbolic name of the diagram to be reordered
 *   target_diagram_name     symbolic name of the new diagram -- result
 *   perm_str                permutation, smth like "124365" (reorders dim-s 3<->4 and 5<->6)
 */
void reorder(char *src_diargam_name, char *target_diagram_name, char *perm_str)
{
    timer_new_entry("reorder", "Tensor transposition (reorder)");
    timer_start("reorder");

    // parse "permutation string"
    int perm[32];
    int success = (str_to_int_array(perm_str, perm) == 0);
    if (!success) {
        errquit("wrong permutation string in reorder: \"%s\"", perm_str);
    }

    // get operand
    assert_diagram_exists(src_diargam_name);
    diagram_t *dg_src = diagram_stack_find(src_diargam_name);

    // perform reordering
    diagram_t *dg_tgt = diagram_reorder(dg_src, perm);

    // save new diagram to stack
    // (replace the old diagram named 'target_diagram_name' if needed)
    strcpy(dg_tgt->name, target_diagram_name);
    if (diagram_stack_find(target_diagram_name) != NULL) {
        diagram_stack_replace(target_diagram_name, dg_tgt);
    }
    else {
        diagram_stack_push(dg_tgt);
    }

    timer_stop("reorder");
}


/*
 * include "template function" reorder_block_<TYPENAME>
 * for real and complex numbers
 */
#define TYPENAME double_complex_t

#include "tensor_transpose.c"

#undef TYPENAME

#define TYPENAME double

#include "tensor_transpose.c"

#undef TYPENAME


/**
 * Apply tensor transposition to a diagram.
 * Returns new diagram with the dimensions reordered (transposed) according to
 * a permutation 'perm'
 */
diagram_t *diagram_reorder(diagram_t *diag, int *perm)
{
    int rank = diag->rank;
    char qparts[CC_DIAGRAM_MAX_RANK];
    char valence[CC_DIAGRAM_MAX_RANK];
    char t3space[CC_DIAGRAM_MAX_RANK];
    char order[CC_DIAGRAM_MAX_RANK];

    // prepare characteristics of a new (transposed) diagram
    for (int i = 0; i < rank; i++) {
        perm[i] = perm[i] - 1;
        qparts[i] = diag->qparts[perm[i]];
        valence[i] = diag->valence[perm[i]] + '0';
        t3space[i] = diag->t3space[perm[i]] + '0';
        order[i] = diag->order[perm[i]] + '0';
    }
    qparts[rank] = '\0';
    valence[rank] = '\0';
    t3space[rank] = '\0';
    order[rank] = '\0';

    // name of the new diagram consist of the old one + "_rdr"
    diagram_t *target_diag = diagram_new("", qparts, valence, t3space, order, diag->only_unique, diag->symmetry);
    sprintf(target_diag->name, "%s_rdr", diag->name);

    int nthreads = 1;
    if (cc_opts->openmp_algorithm == CC_OPENMP_ALGORITHM_EXTERNAL) {
        nthreads = cc_opts->nthreads;
    }
    else {
        nthreads = 1;
    }

#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (size_t iblock = 0; iblock < diag->n_blocks; iblock++) {
        block_t *src_block = diag->blocks[iblock];
        if (src_block->is_unique == 0) {
            continue;
        }

        int spinor_blocks_transposed[CC_DIAGRAM_MAX_RANK];
        for (int idim = 0; idim < rank; idim++) {
            spinor_blocks_transposed[idim] = src_block->spinor_blocks[perm[idim]];
        }

        block_t *target_block = diagram_get_block(target_diag, spinor_blocks_transposed);
        if (src_block->is_unique != target_block->is_unique) {
            printf("reorder(): %d %d\n", src_block->is_unique, target_block->is_unique);
            exit(0);
        }

        reorder_block(src_block, target_block, perm);
    }  // end of loop over blocks

    return target_diag;
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

    reorder_block(uniq_block, b, b->perm_from_unique);

    /*
     * multiply by a sign factor according to a parity of a permutation
     */
    block_load(b);
    if (arith == CC_ARITH_COMPLEX) {
        double complex sign = b->sign + 0.0 * I;
        xscale(CC_DOUBLE_COMPLEX, b->size, b->buf, &sign);
    }
    else {
        double sign = b->sign;
        xscale(CC_DOUBLE, b->size, b->buf, &sign);
    }
    block_store(b);
}


/**
 * Tensor transposition of a block 'source_block' (according to the 'perm' permutation).
 * Result is stored in the 'target_block' block.
 */
void reorder_block(block_t *source_block, block_t *target_block, int *perm)
{
    block_load(source_block);
    block_load(target_block);

    if (arith == CC_ARITH_COMPLEX) {
        TEMPLATE(tensor_transpose, double_complex_t)(source_block->rank, source_block->buf,
                                                     source_block->shape, target_block->shape, perm,
                                                     target_block->buf, cc_opts->nthreads);
    }
    else {
        TEMPLATE(tensor_transpose, double)(source_block->rank, (const double *) source_block->buf,
                                           source_block->shape, target_block->shape, perm,
                                           (double *) target_block->buf, cc_opts->nthreads);
    }

    block_unload(source_block);
    block_store(target_block);
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


void reverse_perm(int n, const int *direct_perm, int *inv_perm)
{
    pair_t p[CC_DIAGRAM_MAX_RANK];

    for (int i = 0; i < n; i++) {
        p[i].i = i;
        p[i].j = direct_perm[i];
    }

    qsort(p, n, sizeof(pair_t), (int (*)(const void *, const void *)) pair_t_less);

    for (int i = 0; i < n; i++) {
        inv_perm[i] = p[i].i;
    }
}



