/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2021 The EXP-T developers.
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

/*******************************************************************************
 * symblock.c
 *
 * Contains implementation of functions and data structures defined in
 * diagrams.h.
 * SymBlock = symmetry block of integrals -- a "building block" for constructing
 * diagrams (recall: the DPD scheme is implemented, so all the diagrams consist
  * of "independent" symmetry blocks).
 * See also symmetry.{c,h} for the evaluation of the direct product of irreps.
 *
 * To learn more about DPD, please see the following papers:
 *                             [basic papers]
 * -> P. Čársky, L.J. Schaad, B.A. Hess, M. Urban, J. Noga.
 *    J. Chem. Phys., V. 87, P. 411 (1987)
 * -> J.F. Stanton, J. Gauss, J.D. Watts, R.J. Bartlett.
 *    J. Chem. Phys., V. 94, P. 4334 (1991)
 *                           [relativistic DPD]
 * -> A. Shee, L. Visscher, T. Saue. J. Chem. Phys., V. 145, P. 184107 (2016)
 * -> L. Visscher, T. Lee, K. Dyall. J. Chem. Phys., V. 105, P. 8769 (1996)
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "platform.h"
#include "datamodel.h"
#include "error.h"
#include "io.h"
#include "memory.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"

// locally used functions
int as_linear_index(int n /* rank */, int *dims, int *idx);
void block_unique(block_t *b, int *qparts, int *valence, int *order);
int64_t block_get_unique_id();


/*******************************************************************************
 * symblock_new
 *
 * constructor of blocks.
 * creates a template for the symblock with dummies for all fields (except ID).
 * arguments:
 *   rank                 rank of the tensor
 *   spinor_blocks_nums   seq numbers of spinor blocks; dimensions will be
 *                        composed from these spinors
 *                        (array of seq int numbers)
*******************************************************************************/
block_t *symblock_new(int rank, int *spinor_blocks_nums, int *qparts, int *valence, int *t3space, int *order, int storage_type,
                      int only_unique)
{
    size_t i, j;
    block_t *block = (block_t *) cc_malloc(1 * sizeof(block_t));

    // id numbers
    block->id = block_get_unique_id();

    // tensor properties
    block->size = 1;
    block->rank = rank;
    block->indices = (int **) cc_malloc(sizeof(int *) * rank);
    for (i = 0; i < rank; i++) {
        size_t iblock = spinor_blocks_nums[i];
        block->indices[i] = (int *) cc_malloc(sizeof(int) * (spinor_blocks[iblock].size + 1));
        block->spinor_blocks[i] = iblock;

        size_t bsize = 0;
        for (j = 0; j < spinor_blocks[iblock].size; j++) {
            // это - цикл по всем спинорам в блоке ! крайне дорого !
            // могут быть готовые блочки по h/p,0/1 в каждом блоке спиноров
            // и анализируя готовый блочок мы сможем сразу отбросить этот симблок еще на старте функции
            // до всех аллокаций

            // *** hole/particle & active/inactive ***

            int rate = 0;
            int ispinor = spinor_blocks[iblock].indices[j];
            int act = is_active(ispinor);
            int occ = is_hole(ispinor);

            if (qparts[i] == 'h' && occ == 1) {
                rate++;
            }
            else if (qparts[i] == 'p' && occ == 0) {
                rate++;
            }

            if (valence[i] == 1 && act == 1) {
                rate++;
            }
            else if (valence[i] == 0 && (act == 0 || act == 1)) {
                rate++;
            }

            if (cc_opts->do_restrict_t3) {
                if (t3space[i] == 1 && is_t3_space_spinor(ispinor)) {
                    rate++;
                }
                else if (t3space[i] == 0) {
                    rate++;
                }
            }

            if ((!cc_opts->do_restrict_t3 && rate == 2) || (cc_opts->do_restrict_t3 && rate == 3)) {
            //if (rate == 2) {
                block->indices[i][1 + bsize] = ispinor;
                bsize++;
            }

        }
        block->indices[i][0] = bsize;   // и тут вот можно получить ноль! и на этом можно закончить
        block->size *= bsize;
        if (block->size == 0) { break; }
    }

    block->is_unique = 1;
    block->sign = 1;
    block->n_equal_perms = 1;
    if (only_unique) {
        block_unique(block, qparts, valence, order);
    }

    // storage of data

    if (block->is_unique == 0) {
        block->storage_type = CC_DIAGRAM_DUMMY;
    }
    else {
        block->storage_type = storage_type;
    }

    if (block->storage_type == CC_DIAGRAM_IN_MEM) {
        block->file_name = NULL;
    }
    if (block->storage_type == CC_DIAGRAM_ON_DISK) {
        block->file_name = (char *) cc_malloc(256);
        sprintf(block->file_name, "block-%ld-%ld.sb", run_id(), block->id);   // будет долго работать!
    }
    if (block->storage_type == CC_DIAGRAM_DUMMY) {
        block->file_name = NULL;
        block->buf = NULL;
    }

    // alloc memory for the buffer
    if (block->storage_type != CC_DIAGRAM_DUMMY) {
        block->buf = (double complex *) cc_calloc(block->size, SIZEOF_WORKING_TYPE);
    }

    // store block on disk if required
    symblock_store(block);

    return block;
}


void transform(int n, int *idx, int *out, int *perm, int shift);


int is_ascending_order(int n, int *a)
{
    for (int i = 1; i < n; i++) {
        if (a[i] < a[i - 1]) {
            return 0;
        }
    }
    return 1;
}


void destroy_block(block_t *b)
{
    if (b->is_unique) {
        printf("destroying unique block!\n");
    }
    b->storage_type = CC_DIAGRAM_DUMMY;
    cc_free(b->buf);
    b->buf = NULL;
}


// order: starts from 1, e.g. 1234, 123456, etc
void block_unique(block_t *b, int *qparts, int *valence, int *order)
{
    int reverse_order[CC_DIAGRAM_MAX_RANK];
    int rank = b->rank;
    reverse_perm(rank, order, reverse_order);

    if (rank == 2) {
        b->is_unique = 1;
        b->sign = 1;
        b->n_equal_perms = 1;
        return;
    }

    // приводим все характеристики к "нормальному порядку"
    int norm_qparts[CC_DIAGRAM_MAX_RANK];
    int norm_valence[CC_DIAGRAM_MAX_RANK];
    int norm_order[CC_DIAGRAM_MAX_RANK];
    int norm_spinor_blocks[CC_DIAGRAM_MAX_RANK];
    for (int i = 0; i < rank; i++) {
        norm_qparts[i] = qparts[reverse_order[i]];
        norm_valence[i] = valence[reverse_order[i]];
        norm_order[i] = order[reverse_order[i]];
        norm_spinor_blocks[i] = b->spinor_blocks[reverse_order[i]];
    }

    // нужно приписать теперь тип каждому индексу
    int norm_types[CC_DIAGRAM_MAX_RANK];
    for (int i = 0; i < rank; i++) {
        if (norm_qparts[i] == 'h' && norm_valence[i] == 0) {
            norm_types[i] = 'h';
        }
        if (norm_qparts[i] == 'p' && norm_valence[i] == 0) {
            norm_types[i] = 'p';
        }
        if (norm_qparts[i] == 'h' && norm_valence[i] == 1) {
            norm_types[i] = 'g';
        }
        if (norm_qparts[i] == 'p' && norm_valence[i] == 1) {
            norm_types[i] = 'v';
        }
    }

    // и в бра, и в кет ВСЕ индексы должны быть одного типа
    int n = rank / 2;
    int one_type_bra = 1;
    int one_type_ket = 1;
    for (int i = 1; i < n; i++) {
        if (norm_types[i] != norm_types[i - 1]) {
            one_type_bra = 0;
            break;
        }
    }
    for (int i = 1; i < n; i++) {
        if (norm_types[n + i] != norm_types[n + i - 1]) {
            one_type_ket = 0;
            break;
        }
    }

    b->is_unique = 1;
    b->sign = 1;
    b->n_equal_perms = 1;

    int bra_sign = 1;
    int ket_sign = 1;
    int bra_unique = 1;
    int ket_unique = 1;
    int bra_n_eq_perm = 1;
    int ket_n_eq_perm = 1;

    typedef struct {
        int p[CC_DIAGRAM_MAX_RANK];
        int sign;
    } perm_t;

    perm_t perms_2[] = {
            {{0, 1}, 1},
            {{1, 0}, -1}
    };
    perm_t perms_3[] = {
            {{0, 1, 2}, 1},
            {{0, 2, 1}, -1},
            {{1, 2, 0}, 1},
            {{1, 0, 2}, -1},
            {{2, 0, 1}, 1},
            {{2, 1, 0}, -1}
    };

    perm_t *perm_list;
    int len_perm_list;
    int P[CC_DIAGRAM_MAX_RANK] = {0, 1, 2, 3, 4, 5, 6, 7};

    if (rank == 4) {
        perm_list = perms_2;
        len_perm_list = 2;
    }
    else if (rank == 6) {
        perm_list = perms_3;
        len_perm_list = 6;
    }
    else {
        errquit("please, add array of permutations for rank-8+ diagrams");
    }

    // bra
    if (one_type_bra) {

        int iperm;
        int *spb = norm_spinor_blocks; // spinor blocks

        for (iperm = 0; iperm < len_perm_list; iperm++) {
            int buf[CC_DIAGRAM_MAX_RANK];
            for (int i = 0; i < n; i++) {
                buf[i] = spb[perm_list[iperm].p[i]];
            }
            if (is_ascending_order(n, buf)) {
                break;
            }
        }
        if (iperm == 0) {
            bra_unique = 1;
            bra_sign = 1;
        }
        else {
            bra_unique = 0;
            bra_sign = perm_list[iperm].sign;
        }
        for (int i = 0; i < n; i++) {
            P[i] = perm_list[iperm].p[i];
        }

        // count the number of blocks which are equal to this one
        assert(rank == 4 || rank == 6);
        if (rank == 4) {
            if (spb[0] == spb[1]) {
                bra_n_eq_perm = 1;
            }
            else {
                bra_n_eq_perm = 2;
            }
        }
        else if (rank == 6) {
            if (spb[0] == spb[1] && spb[1] == spb[2]) {
                bra_n_eq_perm = 1;
            }
            else if (spb[0] == spb[1] || spb[1] == spb[2] || spb[0] == spb[2]) {
                bra_n_eq_perm = 3;
            }
            else {
                bra_n_eq_perm = 6;
            }
        }
    }

    // ket
    if (one_type_ket) {
        int iperm;
        int *spb = norm_spinor_blocks + n; // spinor blocks

        for (iperm = 0; iperm < len_perm_list; iperm++) {
            int buf[CC_DIAGRAM_MAX_RANK];
            for (int i = 0; i < n; i++) {
                buf[i] = spb[perm_list[iperm].p[i]];
            }
            if (is_ascending_order(n, buf)) {
                break;
            }
        }

        if (iperm == 0) {
            ket_unique = 1;
            ket_sign = 1;
        }
        else {
            ket_unique = 0;
            ket_sign = perm_list[iperm].sign;
        }
        for (int i = 0; i < n; i++) {
            P[n + i] = n + perm_list[iperm].p[i];
        }

        // count the number of blocks which are equal to this one
        assert(rank == 4 || rank == 6);
        if (rank == 4) {
            if (spb[0] == spb[1]) {
                ket_n_eq_perm = 1;
            }
            else {
                ket_n_eq_perm = 2;
            }
        }
        else if (rank == 6) {
            if (spb[0] == spb[1] && spb[1] == spb[2]) {
                ket_n_eq_perm = 1;
            }
            else if (spb[0] == spb[1] || spb[1] == spb[2] || spb[0] == spb[2]) {
                ket_n_eq_perm = 3;
            }
            else {
                ket_n_eq_perm = 6;
            }
        }
    }

    if (bra_unique && ket_unique) {
        b->sign = 1;
        b->is_unique = 1;
        b->n_equal_perms = 1;
        b->n_equal_perms = bra_n_eq_perm * ket_n_eq_perm;
    }
    else {
        b->sign = bra_sign * ket_sign;
        b->is_unique = 0;
        b->n_equal_perms = bra_n_eq_perm * ket_n_eq_perm;

        int P0[] = {0, 1, 2, 3, 4, 5, 6, 7};
        transform(rank, P0, P0, reverse_order, 0);
        transform(rank, P0, P0, P, 0);
        transform(rank, P0, b->perm_to_unique, order, 1);
        reverse_perm(rank, b->perm_to_unique, b->perm_from_unique);
    }
}


// применяет перестановку perm к элементам массива idx
// perm может содержать числа от 0+shift до n+shift-1
// т.е. если 0123 => shift = 0
// если 1234 => shift = 1
void transform(int n, int *idx, int *out, int *perm, int shift)
{
    int buf[CC_DIAGRAM_MAX_RANK];

    for (int i = 0; i < n; i++) {
        buf[i] = idx[perm[i] - shift];
    }

    for (int i = 0; i < n; i++) {
        out[i] = buf[i];
    }
}


/*******************************************************************************
 * symblock_gen_indices
 *
 * Generates "compound indices" for each matrix element stored in block.
 * Arguments:
 *   block    block of matrix elements
 *   indices  compound indices for each matrix element (stored linearly!).
 *            it is array of size block->size * block->rank.
 *            [indices_1] [indices_2] ... [indices_size]
 *            NOTE: this function assumes that the 'indices' array HAS ALREADY
 *                  BEEN allocated on heap or on stack
 ******************************************************************************/
void symblock_gen_indices(block_t *block, int *indices)
{
    int i, j;
    int rank = block->rank;
    int *dims = (int *) cc_malloc(sizeof(int) * rank);
    int *comb = (int *) cc_malloc(sizeof(int) * rank);

    // dimension length for each dimension of the tensor
    for (i = 0; i < rank; i++) {
        dims[i] = block->indices[i][0];
    }

    // begin loop: set compound index to the minimal possible
    for (i = 0; i < rank; i++) {
        comb[i] = 0;
    }

    for (i = 0; i < block->size; i++) {

        // transform combination to the compound index & save it
        for (j = 0; j < rank; j++) {
            indices[rank * i + j] = block->indices[j][1 + comb[j]];
        }

        // go to the next combination
        for (j = rank - 1; j >= 0; j--) {
            if (comb[j] == dims[j] - 1) {
                comb[j] = 0;
            }
            else {
                comb[j]++;
                break;
            }
        }
    }

    cc_free(dims);
    cc_free(comb);
}


/*******************************************************************************
 * symblock_delete
 *
 * Destructor. Deallocates memory, closes and removes files.
 ******************************************************************************/
void symblock_delete(block_t *block)
{
    for (size_t i = 0; i < block->rank; i++) {
        cc_free(block->indices[i]);
    }
    cc_free(block->indices);

    if (block->storage_type != CC_DIAGRAM_DUMMY) {
        cc_free(block->buf);
    }

    // delete file with buf:
    if (block->storage_type == CC_DIAGRAM_ON_DISK) {
        io_remove(block->file_name);
        if (block->file_name) {
            cc_free(block->file_name);
        }
    }

    cc_free(block);
}


/*******************************************************************************
 * symblock_copy_data
 *
 * Copies all matrix elements from the 'src' symblock to the 'dst'.
 * NOTE: NO CHECKS ARE PERFORMED! (for better performance).
 * TODO: copy data in case of they stored on disk.
 * possible problems: overlapping memory regions
 ******************************************************************************/
void symblock_copy_data(block_t *dst, block_t *src)
{
    symblock_load(dst);
    symblock_load(src);
    if (dst->storage_type == CC_DIAGRAM_DUMMY) {

    }
    else {
        memcpy(dst->buf, src->buf, SIZEOF_WORKING_TYPE * src->size);
    }
    memcpy(dst->spinor_blocks, src->spinor_blocks, sizeof(int) * CC_DIAGRAM_MAX_RANK);
    symblock_store(dst);
    symblock_unload(src);
}


/**
 * clears symmetry block: set all elements to zero
 */
void symblock_clear(block_t *block)
{
    if (block->storage_type == CC_DIAGRAM_DUMMY) {
        return;
    }

    symblock_load(block);
    memset(block->buf, 0, block->size * SIZEOF_WORKING_TYPE);
    symblock_store(block);
}


/**
 * returns dimensions of the block (array of dim = rank)
 */
void symblock_get_dims(block_t *block, int *dims)
{
    for (int i = 0; i < block->rank; i++) {
        dims[i] = block->indices[i][0];
    }
}


/*******************************************************************************
 * as_linear_index
 *
 * converts compound index into linear one
 * all indices must be >= 0
 * n - number of dimensions in the tensor
 *
 * [Fortran version comment]
 * Example: 4-dimensional tensor with dimensions [1:L,1:M,1:N,1:K].
 * For element (l,i,j,k) linear index is equal to
 * (((M*(l-1)+(i-1)))*N+(j-1))*K+k ==
 *  == MNK*(l-1) + NK*(i-1) + K*(j-1) + (k-1) + 1
 * (+1 since we use this address to index Fortran array)
 *
 * [C version comment]
 * Example: 4-dimensional tensor with dimensions [0:L-1,0:M-1,0:N-1,0:K-1].
 * For element (l,i,j,k) linear index is equal to
 * (((M*l+i))*N+j*K+k ==
 *  == MNK*l + NK*i + K*j + k
 ******************************************************************************/
int as_linear_index(int n /* rank */, int *dims, int *idx)
{
    int i, j;
    int lin_index;
    int coef;

    lin_index = 0;
    for (i = 0; i < n; i++) {
        coef = 1;
        for (j = i + 1; j < n; j++) {
            coef *= dims[j];
        }
        lin_index += coef * idx[i];
    }

    return lin_index;
}


/*******************************************************************************
 * as_compound_index
 *
 * converts linear index into compound one (is inverse of the 'as_linear_index')
 * all indices must be >= 0
 * n - number of dimensions in the tensor
 * lin_idx -- linear index (input)
 * idx[n] -- compound index (output)
 ******************************************************************************/
void as_compound_index(int n /* rank */, int *dims, int lin_idx, int *idx)
{
    int i, j;
    int coefs[CC_DIAGRAM_MAX_RANK];

    // calculate coefficients
    for (i = 0; i < n; i++) {
        coefs[i] = 1;
        for (j = i + 1; j < n; j++) {
            coefs[i] *= dims[j];
        }
    }

    for (i = 0; i < n; i++) {
        idx[i] = lin_idx / coefs[i];
        lin_idx = lin_idx % coefs[i];
    }
}


/*******************************************************************************
 * symblock_set
 *
 * set the value of the matrix element with compound index 'idx'
 * NOTE: len(idx) must be equal to diagram's rank! (NO CHECK IS PERFORMED)
 ******************************************************************************/
void symblock_set(block_t *block, double complex val, int *idx)
{
    int dims[CC_DIAGRAM_MAX_RANK];
    int rel_idx[CC_DIAGRAM_MAX_RANK];
    int i, j;

    for (i = 0; i < block->rank; i++) {
        dims[i] = block->indices[i][0];
        for (j = 0; j < dims[i]; j++) {
            if (idx[i] == block->indices[i][1 + j]) {
                rel_idx[i] = j;
                break;
            }
        }
        // no such element
        if (j == dims[i]) {
            return;
        }
    }

    // obtain the linear index of the element in the buffer
    int li = as_linear_index(block->rank, dims, rel_idx);

    // set value
    if (carith) {
        block->buf[li] = val;
    }
    else {
        ((double *) block->buf)[li] = creal(val);
    }
}


/*******************************************************************************
 * symblock_get
 *
 * get the value of the matrix element with compound index 'idx'
 * NOTE: len(idx) must be equal to diagram's rank! (NO CHECK IS PERFORMED)
 ******************************************************************************/
double complex symblock_get(block_t *block, int *idx)
{
    int dims[CC_DIAGRAM_MAX_RANK];
    int rel_idx[CC_DIAGRAM_MAX_RANK];
    int i, j;

    for (i = 0; i < block->rank; i++) {
        dims[i] = block->indices[i][0];
        for (j = 0; j < dims[i]; j++) {
            if (idx[i] == block->indices[i][1 + j]) {
                rel_idx[i] = j;
                break;
            }
        }
        // no such element
        if (j == dims[i]) {
            return 0.0 + 0.0 * I;
        }
    }

    int li = as_linear_index(block->rank, dims, rel_idx);

    if (carith) {
        return block->buf[li];
    }
    else {
        return ((double *) block->buf)[li] + 0.0 * I;
    }
}


void symblock_load(block_t *block)
{
    int f;

    if (block->storage_type != CC_DIAGRAM_ON_DISK) {
        return;
    }

    block->buf = (double complex *) cc_malloc(block->size * SIZEOF_WORKING_TYPE);

    f = io_open(block->file_name, "r");
    if (f == -1) {
        errquit("-1 in load, unable to open block file %s\n", block->file_name);
    }

    io_read_compressed(f, block->buf, block->size * SIZEOF_WORKING_TYPE);

    io_close(f);
}


void symblock_unload(block_t *block)
{
    if (block->storage_type != CC_DIAGRAM_ON_DISK) {
        return;
    }

    cc_free(block->buf);
    block->buf = NULL;
}


void symblock_store(block_t *block)
{
    int f;

    if (block->storage_type != CC_DIAGRAM_ON_DISK) {
        return;
    }

    f = io_open(block->file_name, "w");
    if (f == -1) {
        errquit("-1 in store name = %s\n", block->file_name);
    }

    io_write_compressed(f, block->buf, block->size * SIZEOF_WORKING_TYPE);
    io_close(f);

    cc_free(block->buf);
    block->buf = NULL;
}


void symblock_print(block_t *block)
{
    size_t i, j;
    int *indices = (int *) cc_malloc(block->size * block->rank * sizeof(int));

    symblock_load(block);

    printf("  --------------------------------------\n");
    printf("  SYMBLOCK\n");
    printf("    ID      %ld\n", block->id);
    printf("    rank    %d\n", block->rank);
    printf("    unique  %d\n", block->is_unique);
    printf("    #eqperm %d\n", block->n_equal_perms);
    printf("    sign    %d\n", block->sign);
    if (block->is_unique == 0) {
        printf("    P->uniq ");
        for (int i = 0; i < block->rank; i++) {
            printf("%d ", block->perm_to_unique[i]);
        }
        printf("\n");
        printf("    P<-uniq ");
        for (int i = 0; i < block->rank; i++) {
            printf("%d ", block->perm_from_unique[i]);
        }
        printf("\n");
    }

    printf("    spi blk ");
    for (i = 0; i < block->rank; i++) {
        printf("%5d", block->spinor_blocks[i]);
    }
    printf("\n");
    printf("    indices:\n");
    for (i = 0; i < block->rank; i++) {
        printf("      [%d] ", block->indices[i][0]);
        for (j = 0; j < block->indices[i][0]; j++) {
            printf("%d ", block->indices[i][j + 1]);
        }
        printf("\n");
    }
    printf("    size = %ld\n", block->size);

    symblock_gen_indices(block, indices);
    for (i = 0; i < block->size; i++) {
        for (j = 0; j < block->rank; j++) {
            printf("%2d", indices[block->rank * i + j]);
        }
        if (!carith) {
            double val = ((double *) block->buf)[i];
            printf("%12.4e ", (fabs(val) > 1e-10) ? val : 0.0);
            if ((i + 1) % 8 == 0) {
                printf("\n");
            }
            continue;
        }
        // 5 numbers if complex arithmetics and 8 if real
        if (arith == CC_ARITH_COMPLEX) {
            printf("%12.8f%12.8f   ", creal(block->buf[i]), cimag(block->buf[i]));
            if ((i + 1) % 5 == 0) {
                printf("\n");
            }
        }
        else {
            printf("%12.4e ", creal(block->buf[i]));
            if ((i + 1) % 8 == 0) {
                printf("\n");
            }
        }

    }
    printf("\n");

    symblock_unload(block);

    cc_free(indices);
}


/* each entry on disk will contain 64 x int64_t for future applications */
#define LEN_RAW_BUF 16


/**
 * symblock_write
 *
 * Write symblock data to the binary file.
 * Argument 'fd' = file descriptor.
 *
 */
void symblock_write(int fd, block_t *block)
{
    int i;
    int n;
    static int64_t raw_buf[LEN_RAW_BUF];

    io_write_compressed(fd, &block->id, sizeof(block->id));
    io_write_compressed(fd, &block->rank, sizeof(block->rank));

    // information about uniqueness
    io_write_compressed(fd, &block->is_unique, sizeof(block->is_unique));
    io_write_compressed(fd, &block->sign, sizeof(block->sign));
    io_write_compressed(fd, &block->n_equal_perms, sizeof(block->n_equal_perms));
    io_write_compressed(fd, block->perm_from_unique, sizeof(block->perm_from_unique));
    io_write_compressed(fd, block->perm_to_unique, sizeof(block->perm_to_unique));

    // raw buffer -- for future applications and backward binary compatibility
    memset(raw_buf, 0, sizeof(raw_buf));
    io_write_compressed(fd, raw_buf, sizeof(raw_buf));

    io_write_compressed(fd, block->spinor_blocks, sizeof(int) * CC_DIAGRAM_MAX_RANK);

    // indices
    for (i = 0; i < block->rank; i++) {
        n = block->indices[i][0];
        io_write_compressed(fd, &n, sizeof(n));
        io_write_compressed(fd, &block->indices[i][1], sizeof(int) * n);
    }

    // size & data
    io_write_compressed(fd, &block->size, sizeof(block->size));
    io_write_compressed(fd, &block->storage_type, sizeof(block->storage_type));
    if (block->storage_type == CC_DIAGRAM_IN_MEM) {
        io_write_compressed(fd, block->buf, SIZEOF_WORKING_TYPE * block->size);
    }
    else if (block->storage_type == CC_DIAGRAM_DUMMY) {
        // nothing
    }
    else {   // CC_DIAGRAM_ON_DISK
        n = strlen(block->file_name);
        io_write_compressed(fd, &n, sizeof(n));
        io_write_compressed(fd, block->file_name, n * sizeof(char)); // NOTE: without '\0'!
    }
}


/*******************************************************************************
 * symblock_read
 *
 * Reads symblock from the binary file. Returns (newly allocated) symblock.
 * Argument 'fd' = file descriptor.
 ******************************************************************************/
block_t *symblock_read(int fd)
{
    int i;
    int n;
    block_t *block;
    static int64_t raw_buf[LEN_RAW_BUF];

    block = (block_t *) cc_malloc(1 * sizeof(block_t));

    io_read_compressed(fd, &block->id, sizeof(block->id));
    io_read_compressed(fd, &block->rank, sizeof(block->rank));

    // information about uniqueness
    io_read_compressed(fd, &block->is_unique, sizeof(block->is_unique));
    io_read_compressed(fd, &block->sign, sizeof(block->sign));
    io_read_compressed(fd, &block->n_equal_perms, sizeof(block->n_equal_perms));
    io_read_compressed(fd, block->perm_from_unique, sizeof(block->perm_from_unique));
    io_read_compressed(fd, block->perm_to_unique, sizeof(block->perm_to_unique));

    // raw buffer -- for future applications and backward binary compatibility
    io_read_compressed(fd, raw_buf, sizeof(raw_buf));

    // info about spinor blocks
    io_read_compressed(fd, block->spinor_blocks, sizeof(int) * CC_DIAGRAM_MAX_RANK);

    // set new unique block's ID
    block->id = block_get_unique_id();

    // indices
    block->indices = (int **) cc_malloc(sizeof(int *) * block->rank);
    for (i = 0; i < block->rank; i++) {
        io_read_compressed(fd, &n, 1 * sizeof(int));
        block->indices[i] = (int *) cc_malloc(sizeof(int) * (n + 1));
        block->indices[i][0] = n;
        io_read_compressed(fd, &block->indices[i][1], sizeof(int) * n);
    }

    // size & data
    io_read_compressed(fd, &block->size, sizeof(block->size));
    io_read_compressed(fd, &block->storage_type, sizeof(block->storage_type));
    if (block->storage_type == CC_DIAGRAM_IN_MEM) {
        block->buf = (double complex *) cc_malloc(SIZEOF_WORKING_TYPE * block->size);
        io_read_compressed(fd, block->buf, SIZEOF_WORKING_TYPE * block->size);
    }
    else if (block->storage_type == CC_DIAGRAM_DUMMY) {
        block->buf = NULL;
    }
    else {   // CC_DIAGRAM_ON_DISK
        io_read_compressed(fd, &n, sizeof(int));
        block->file_name = (char *) cc_malloc(sizeof(char) * (n + 1));
        io_read_compressed(fd, block->file_name, n * sizeof(char)); // NOTE: without '\0'!
        block->file_name[n] = '\0';
    }

    return block;
}


/*
 * Helper functions
 */
int64_t block_get_unique_id()
{
    // counter of symmetry blocks (for unique ID's)
    static int64_t blocks_count = 0;

    int64_t id = blocks_count;
    blocks_count++;
    return id;
}
