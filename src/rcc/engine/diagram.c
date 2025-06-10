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
 * Contains implementation of functions and data structures defined in
 * diagrams.h.
 * Diagram -- main data structure.
 * Implements the DPD scheme, consists of "independent" symmetry blocks.
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
 */

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io.h"
#include "dgstack.h"
#include "diagram.h"
#include "error.h"
#include "memory.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"

static int cmp_pair_t(const void *op1, const void *op2);

static void inverse_perm(int n, int *perm, int *out);

int guess_rank(char *qparts, char *valence, char *order);

void parse_qp_string(int rank, char *qp_str, int *qp);

void parse_valence_string(int rank, char *val_str, int *val);

void parse_order_string(int rank, char *order_str, int *order);

int guess_storage_class(char *name, int rank, char *qparts, char *valence);

void diagram_init_inverse_index(diagram_t *dg, size_t n_blocks, block_t **block_list);

void diagram_bind_blocks(diagram_t *dg, size_t n_blocks, block_t **block_list);


/**
 * Constructor of diagrams.
 *
 *
 * @param name     name of the diagram
 * @param qparts   list of quasiparticles in "natural form": 'hp', 'pphh', etc
 * @param valence  valence/non-valence flags in "natural form": '00', '0100', etc
 * @param order    reordering from "natural form". in fact no reordering is
 *                 performed since all integrals in the template of the diagram
 *                 are zero.
 *                 examples: '12' (natural), '1234' (nat), '3412', '123465', etc
 *
 * @note len(qparts) == len(valence) == len(order) == rank of the diagram
 * @note new diagram is not binded to the singly-linked list "dg_stack"
 */
diagram_t *diagram_new(char *name, char *qparts, char *valence, char *t3space, char *order, int perm_unique, int irrep)
{
    int qparts_arr[CC_DIAGRAM_MAX_RANK];
    int valence_arr[CC_DIAGRAM_MAX_RANK];
    int t3space_arr[CC_DIAGRAM_MAX_RANK];
    int order_arr[CC_DIAGRAM_MAX_RANK];
    size_t i;
    block_t **block_list;  // temporary storage for the created symmetry blocks
    int max_sbs;            // max number of sym blocks = (# spinor blocks)^rank
    int blocks_counter = 0;          // counter of created blocks
    int ijkl[CC_DIAGRAM_MAX_RANK];  // indices for the arbitrary-nested loop
    int only_unique = perm_unique;

    // check arguments for correctness and pre-process them
    int rank = guess_rank(qparts, valence, order);
    parse_qp_string(rank, qparts, qparts_arr);
    parse_valence_string(rank, valence, valence_arr);
    parse_valence_string(rank, t3space, t3space_arr);
    parse_order_string(rank, order, order_arr);

    // allocate memory for all the metainfo and the diagram itself
    diagram_t *dg = (diagram_t *) cc_malloc(1 * sizeof(diagram_t));

    // fill diagram's fields
    //dg->dg_id = diagram_get_unique_id();
    strncpy(dg->name, name, CC_DIAGRAM_MAX_NAME);
    dg->name[CC_DIAGRAM_MAX_NAME - 1] = '\0';
    dg->rank = rank;
    dg->symmetry = irrep;
    dg->only_unique = only_unique;
    intcpy(dg->qparts, qparts_arr, rank);
    intcpy(dg->valence, valence_arr, rank);
    intcpy(dg->t3space, t3space_arr, rank);
    intcpy(dg->order, order_arr, rank);

    // permutation which is inverse to the 'order' perm-n
    // is required for the subsequent application of the DPD scheme
    int reverse_order[CC_DIAGRAM_MAX_RANK];
    inverse_perm(rank, order_arr, reverse_order);

    int storage_type = guess_storage_class(name, rank, qparts, valence);

    // create symmetry blocks
    max_sbs = (int) pow(n_spinor_blocks, rank);
    block_list = (block_t **) cc_malloc(sizeof(block_t *) * max_sbs);

    int is_fully_sym = (irrep == get_totally_symmetric_irrep()) ? 1 : 0;

    // dynamically-nested loop over tensor's dimensions
    for (i = 0; i < rank; i++) {
        ijkl[i] = 0;
    }
    while (ijkl[0] < n_spinor_blocks) {

        // check symmetry for the block (direct product decomposition scheme)
        int sym[CC_DIAGRAM_MAX_RANK];
        for (size_t i = 0; i < rank; i++) {
            size_t iblock = ijkl[reverse_order[i] - 1];
            sym[i] = spinor_blocks[iblock].repno;
        }

        if (is_fully_sym) {
            // check if this block is symmetry-allowed
            if (!dpd_contains_totsym_rank[rank / 2 - 1](sym)) {
                goto next_symblock;
            }
        }
        else {
            // check if this block is symmetry-allowed
            if (!dpd_contains_totsym_rank_op_sym[rank / 2 - 1](sym, irrep)) {
                goto next_symblock;
            }
        }

        int is_zero = is_symblock_zero(rank, ijkl, dg->qparts, dg->valence, dg->t3space);
        if (is_zero) {
            goto next_symblock;
        }

        // create empty block
        block_t *block = block_new(rank, ijkl, dg->qparts, dg->valence, dg->t3space, dg->order, storage_type,
                                   only_unique);
        if (block == NULL) {
            goto next_symblock;
        }
        /*else {
            if (rank == 6)
            printf("nonzero!\n");
        }*/

        // save this block (and count it)
        block_list[blocks_counter] = block;
        blocks_counter++;

        // proceed to the next iteration: new set of indices
        next_symblock:
        ijkl[rank - 1]++;
        for (i = rank - 1; i > 0; i--) {
            if (ijkl[i] >= n_spinor_blocks) {
                ijkl[i] = 0;
                ijkl[i - 1]++;
            }
            else {
                break;
            }
        }
        // new unique combination of spinor blocks is constructed.
        // (next set of indices)
    }

    diagram_bind_blocks(dg, blocks_counter, block_list);

    // cleanup
    cc_free(block_list);

    return dg;
}


/**
 * Binds blocks to the diagram.
 */
void diagram_bind_blocks(diagram_t *dg, size_t n_blocks, block_t **block_list)
{
    dg->n_blocks = n_blocks;
    dg->blocks = (block_t **) cc_malloc(sizeof(block_t *) * n_blocks);

    for (size_t i = 0; i < n_blocks; i++) {
        dg->blocks[i] = block_list[i];
    }

    diagram_init_inverse_index(dg, dg->n_blocks, dg->blocks);
}


/**
 * Adds blocks to the inverted index.
 */
void diagram_init_inverse_index(diagram_t *dg, size_t n_blocks, block_t **block_list)
{
    int dims[CC_DIAGRAM_MAX_RANK];
    size_t rank = dg->rank;

    for (int i = 0; i < CC_DIAGRAM_MAX_RANK; i++) {
        dims[i] = n_spinor_blocks;
    }

    size_t inv_index_size = int_pow(n_spinor_blocks, dg->rank);
    dg->inv_index = (size_t *) cc_malloc(inv_index_size * sizeof(size_t));
    memset(dg->inv_index, 0, inv_index_size * sizeof(size_t));

    for (size_t i = 0; i < n_blocks; i++) {
        block_t *block = block_list[i];
        size_t ii = tensor_index_to_linear(rank, dims, block->spinor_blocks);
        dg->inv_index[ii] = i;
    }
}


/**
 * returns list of valence(1)/general(0) labels for each
 * dimension of the diagram
 */
void diagram_get_valence(diagram_t *diag, char *valence)
{
    for (size_t i = 0; i < diag->rank; i++) {
        valence[i] = '0' + diag->valence[i];
    }
    valence[diag->rank] = '\0';
}


/**
 * returns list of t3space(1)/general(0) labels for each
 * dimension of the diagram
 */
void diagram_get_t3space(diagram_t *diag, char *t3space)
{
    for (size_t i = 0; i < diag->rank; i++) {
        t3space[i] = '0' + diag->t3space[i];
    }
    t3space[diag->rank] = '\0';
}


/**
 * returns list of hole(h)/particle(p) labels for each
 * dimension of the diagram
 */
void diagram_get_quasiparticles(diagram_t *diag, char *qparts)
{
    for (size_t i = 0; i < diag->rank; i++) {
        qparts[i] = diag->qparts[i];
    }
    qparts[diag->rank] = '\0';
}


/**
 * for the given diagram returns "permutation from the normal order"
 */
void diagram_get_order(diagram_t *diag, char *order)
{
    for (size_t i = 0; i < diag->rank; i++) {
        order[i] = '0' + diag->order[i];
    }
    order[diag->rank] = '\0';
}


/**
 * returns memory (bytes) used by the diagram (excluding metainfo, only integrals)
 */
void diagram_get_memory_used(diagram_t *dg, size_t *ram_used, size_t *disk_used)
{
    *ram_used = 0;
    *disk_used = 0;

    for (size_t ib = 0; ib < dg->n_blocks; ib++) {
        block_t *block = dg->blocks[ib];
        if (block->storage_type == CC_DIAGRAM_IN_MEM) {
            *ram_used += block->size * SIZEOF_WORKING_TYPE;
        }
        else { // CC_DIAGRAM_ON_DISK
            *disk_used += block->size * SIZEOF_WORKING_TYPE;
        }
    }
}


/**
 * get number of a block (its index)
 * using inverted index
 *
 * BUG: returns wrong pointer! to the moment only the block_index parameter
 * can be used safely!
 */
block_t *diagram_get_block(diagram_t *dg, int *spinor_blocks_nums)
{
    static int dims[CC_DIAGRAM_MAX_RANK];
    static int initialized = 0;

    if (initialized == 0) {
        initialized = 1;
        for (int i = 0; i < CC_DIAGRAM_MAX_RANK; i++) {
            dims[i] = n_spinor_blocks;
        }
    }

    /*
     * Some diagrams can contain no integrals due to symmetry reasons
     * (sometimes for modest-size problems)
     */
    if (dg->n_blocks == 0) {
        return NULL;
    }

    size_t ii = tensor_index_to_linear(dg->rank, dims, spinor_blocks_nums);
    size_t block_index = dg->inv_index[ii];
    block_t *ret = dg->blocks[block_index];

    for (int i = 0; i < dg->rank; i++) {
        if (ret->spinor_blocks[i] != spinor_blocks_nums[i]) {
            return NULL;
        }
    }

    return ret;
}


/**
 * destructor -- deallocates all memory associated with this object.
 */
void diagram_delete(diagram_t *dg)
{
    for (size_t i = 0; i < dg->n_blocks; i++) {
        block_delete(dg->blocks[i]);
    }
    cc_free(dg->blocks);
    cc_free(dg->inv_index);
    cc_free(dg);
}


/**
 * "Copy constructor".
 * Creates a copy of a diagram 'dg' (the same h/p, valence, order, indices, data
 * BUT the other IDs, names of files with integrals, etc)
 */
diagram_t *diagram_copy(diagram_t *dg)
{
    diagram_t *clone;
    char qparts[CC_DIAGRAM_MAX_RANK];
    char valence[CC_DIAGRAM_MAX_RANK];
    char t3space[CC_DIAGRAM_MAX_RANK];
    char order[CC_DIAGRAM_MAX_RANK];

    diagram_get_quasiparticles(dg, qparts);
    diagram_get_valence(dg, valence);
    diagram_get_t3space(dg, t3space);
    diagram_get_order(dg, order);

    clone = diagram_new(dg->name, qparts, valence, t3space, order, dg->only_unique, dg->symmetry);

    // name of the new diagram consist of the old one + "_copy_" + clone's ID
    sprintf(clone->name, "%s_copy", dg->name);//, clone->dg_id);

    // copy data: block-by-block
    for (size_t i = 0; i < dg->n_blocks; i++) {
        block_copy_data(clone->blocks[i], dg->blocks[i]);
    }

    return clone;
}


/**
 * Storage type: in memory or on disk.
 * Diagram is assumed to be stored on disk even in case only one block is
 * stored on disk.
 */
storage_type_t diagram_get_storage_type(diagram_t *dg)
{
    for (int i = 0; i < dg->n_blocks; i++) {
        if (dg->blocks[i]->storage_type == CC_DIAGRAM_ON_DISK) {
            return CC_DIAGRAM_ON_DISK;
        }
    }

    return CC_DIAGRAM_IN_MEM;
}


/**
 * TODO: error handling (only errquit is now implemented)
 */
double complex diagram_get(diagram_t *dg, int *idx)
{
    int rank = dg->rank;
    int sp_blocks[CC_DIAGRAM_MAX_RANK];
    int i, is;

    for (i = 0; i < rank; i++) {
        sp_blocks[i] = spinor_info[idx[i]].blockno;
    }

    // find required block

    // TODO: tests for this new version

    block_t *b = diagram_get_block(dg, sp_blocks);
    if (b == NULL) {
        return 0.0 + 0.0 * I;
    }

    if (b->is_unique == 0) {
        int idx_uniq[CC_DIAGRAM_MAX_RANK];
        int uniq_spinor_blocks[CC_DIAGRAM_MAX_RANK];

        transform(b->rank, b->spinor_blocks, uniq_spinor_blocks, b->perm_to_unique, 0);
        transform(b->rank, idx, idx_uniq, b->perm_to_unique, 0);

        block_t *uniq_block = diagram_get_block(dg, uniq_spinor_blocks);

        return b->sign * block_get_element(uniq_block, idx_uniq);
    }
    else {
        return block_get_element(b, idx);
    }
}


double complex diagram_get_2(diagram_t *dg, int i, int j)
{
    int idx[2];
    idx[0] = i;
    idx[1] = j;

    return diagram_get(dg, idx);
}


double complex diagram_get_4(diagram_t *dg, int i, int j, int k, int l)
{
    int idx[4];
    idx[0] = i;
    idx[1] = j;
    idx[2] = k;
    idx[3] = l;

    return diagram_get(dg, idx);
}


/**
 * set value of the matrix element with the compound index 'idx'
 *
 * @param new val value of the matrix element
 * @param idx compound index of matrix element
 *
 * @note len(idx) must be equal to diagram's rank! (NO CHECK IS PERFORMED)
 */
void diagram_set(diagram_t *dg, double complex val, int *idx)
{
    int rank = dg->rank;
    int spinor_blocks[CC_DIAGRAM_MAX_RANK];

    // extract symmetry of all indices
    for (int i = 0; i < rank; i++) {
        spinor_blocks[i] = spinor_info[idx[i]].blockno;
    }

    block_t *block = diagram_get_block(dg, spinor_blocks);
    if (block != NULL && block->is_unique) {
        block_set_element(block, val, idx);
    }
}


void diagram_set_4(diagram_t *dg, double complex val, int i, int j, int k, int l)
{
    int idx[4];
    idx[0] = i;
    idx[1] = j;
    idx[2] = k;
    idx[3] = l;

    diagram_set(dg, val, idx);
}


/**
 * Prints diagram 'dg' to stdout (both metainfo & all integrals).
 * May cause extremely huge output; only for debug purposes.
 */
void diagram_debug_print(diagram_t *dg)
{
    printf("=================\n");
    printf("DIAGRAM %s\n", dg->name);
    printf("  rank = %d\n", dg->rank);
    printf("  symmetry = %d (%s)\n", dg->symmetry, get_irrep_name(dg->symmetry));

    printf("  quasiparticles: ");
    for (int i = 0; i < dg->rank; i++) {
        printf("%c", dg->qparts[i]);
    }
    printf("\n");

    printf("  valence: ");
    for (int i = 0; i < dg->rank; i++) {
        printf("%1d", dg->valence[i]);
    }
    printf("\n");

    printf("  t3space: ");
    for (int i = 0; i < dg->rank; i++) {
        printf("%1d", dg->t3space[i]);
    }
    printf("\n");

    printf("  order: ");
    for (int i = 0; i < dg->rank; i++) {
        printf("%1d", dg->order[i]);
    }
    printf("\n");

    printf("number of non-zero blocks = %ld\n", dg->n_blocks);
    for (int i = 0; i < dg->n_blocks; i++) {
        printf("number = %d/%d\n", i + 1, dg->n_blocks);
        if (dg->blocks[i]->is_unique == 0) {
            printf("non unique, restoring...\n");
            restore_block(dg, dg->blocks[i]);
        }
        block_debug_print(dg->blocks[i]);
    }
    printf("=================\n");
}


/**
 * prints brief information about the diagram
 */
void diagram_summary(diagram_t *diag)
{
    char str_qparts[CC_DIAGRAM_MAX_RANK + 1];
    char str_valence[CC_DIAGRAM_MAX_RANK + 1];
    char str_t3space[CC_DIAGRAM_MAX_RANK + 1];
    char str_order[CC_DIAGRAM_MAX_RANK + 1];
    size_t n_unique = 0;

    diagram_get_quasiparticles(diag, str_qparts);
    diagram_get_valence(diag, str_valence);
    diagram_get_t3space(diag, str_t3space);
    diagram_get_order(diag, str_order);

    for (size_t i = 0; i < diag->n_blocks; i++) {
        block_t *b = diag->blocks[i];
        if (b->is_unique == 1) {
            n_unique++;
        }
    }

    printf(" diagram %s: irrep=%d(%s) %s %s %s %s %d/%d\n", diag->name, diag->symmetry,
           get_irrep_name(diag->symmetry), str_qparts, str_valence, str_t3space, str_order, n_unique, diag->n_blocks);
}


/**
 * Sets all matrix elements to zero.
 * NOTE: All diagrams must be allocated and initialized!
 */
diagram_t *diagram_clear(diagram_t *dg)
{
    assert(dg != NULL);

    for (size_t isb = 0; isb < dg->n_blocks; isb++) {
        block_clear(dg->blocks[isb]);
    }

    return dg;
}


void set_order(char *dg_name, char *new_order)
{
    diagram_t *dg = diagram_stack_find(dg_name);
    if (dg == NULL) {
        errquit("set_order(): diagram '%s' not found", dg_name);
    }

    assert(strlen(new_order) >= dg->rank);

    for (int i = 0; i < dg->rank; i++) {
        dg->order[i] = new_order[i] - '0';
    }
}


/**
 * diagram_write_binary
 *
 * Writes diagram to disk (binary format). The name of the file is 'dg->name.dg'
 * The platform-independent I/O interface is used (prefixes 'io_')
 * (see io.h and io.c for details)
 */
void diagram_write_binary(diagram_t *dg, char *file_name)
{
    int f = io_open(file_name, "w");

    // diagram's name (len = CC_DIAGRAM_MAX_NAME)
    io_write(f, dg->name, CC_DIAGRAM_MAX_NAME * sizeof(char));

    // diagram's rank
    io_write_compressed(f, &dg->rank, sizeof(dg->rank));

    // symmetry of the operator represented by the diagram
    io_write_compressed(f, &dg->symmetry, sizeof(dg->symmetry));

    // about uniqueness
    io_write_compressed(f, &dg->only_unique, sizeof(dg->only_unique));

    // quasiparticle type list (h/p), valence, order
    io_write(f, dg->qparts, sizeof(int) * dg->rank);
    io_write(f, dg->valence, sizeof(int) * dg->rank);
    io_write(f, dg->t3space, sizeof(int) * dg->rank);
    io_write(f, dg->order, sizeof(int) * dg->rank);

    // inverted index
    io_write(f, dg->inv_index, sizeof(size_t) * int_pow(n_spinor_blocks, dg->rank));

    // total number of symmetry blocks in this diagram
    io_write(f, &dg->n_blocks, sizeof(dg->n_blocks));

    for (int i = 0; i < dg->n_blocks; i++) {
        block_write_binary(f, dg->blocks[i]);
    }

    io_close(f);
}


/**
 * diagram_read_binary
 *
 * Reads diagram from disk (from binary file).
 * The platform-independent I/O interface is used (prefixes 'io_')
 * (see io.h and io.c for details)
 */
diagram_t *diagram_read_binary(char *file_name)
{
    int f = io_open(file_name, "r");

    if (f == -1) {
        printf(" Unable to open diagram file %s\n", file_name);
        return NULL;
    }

    diagram_t *dg = (diagram_t *) cc_malloc(1 * sizeof(diagram_t));

    // diagram's name (len = CC_DIAGRAM_MAX_NAME)
    io_read(f, dg->name, CC_DIAGRAM_MAX_NAME * sizeof(char));

    // diagram's rank
    io_read_compressed(f, &dg->rank, sizeof(dg->rank));

    // symmetry of the operator represented by the diagram
    io_read_compressed(f, &dg->symmetry, sizeof(dg->symmetry));

    // about uniqueness
    io_read_compressed(f, &dg->only_unique, sizeof(dg->only_unique));

    // quasiparticle type list (h/p), valence, order
    io_read(f, dg->qparts, sizeof(int) * dg->rank);
    io_read(f, dg->valence, sizeof(int) * dg->rank);
    io_read(f, dg->t3space, sizeof(int) * dg->rank);
    io_read(f, dg->order, sizeof(int) * dg->rank);

    // inverted index
    dg->inv_index = (size_t *) cc_malloc(sizeof(size_t) * int_pow(n_spinor_blocks, dg->rank));
    io_read(f, dg->inv_index, sizeof(size_t) * int_pow(n_spinor_blocks, dg->rank));

    // total number of symmetry blocks in this diagram
    io_read(f, &dg->n_blocks, sizeof(dg->n_blocks));

    // alloc memory for blocks pointers & then read blocks
    dg->blocks = (block_t **) cc_malloc(dg->n_blocks * sizeof(block_t *));
    for (int i = 0; i < dg->n_blocks; i++) {
        dg->blocks[i] = block_read_binary(f);
    }

    io_close(f);

    // try to find in the stack diagram with the same name. if found -- replace it
    // with the new diagram, if not found -- add new diagram to the stack
    int idx = diagram_stack_find_index(dg->name);
    if (idx == -1) {
        diagram_stack_push(dg);
    }
    else {
        diagram_stack_replace(dg->name, dg);
    }

    return dg;
}


/**
 * writes amplitudes (or matrix elements) which are contained in the given diagram
 * to the formatted file.
 * for example,
 *
 * singles:
 * <i> <a> <t_ia>
 * doubles:
 * <i> <j> <a> <b> <t_ijab>
 */
void diagram_write_formatted(diagram_t *dg, char *file_name)
{
    FILE *txt_file;

    txt_file = fopen(file_name, "w");
    if (txt_file == NULL) {
        errquit("cannot open txt file '%s'", file_name);
    }

    // loop over blocks
    for (int iblock = 0; iblock < dg->n_blocks; iblock++) {

        block_t *block = dg->blocks[iblock];

        if (block->is_unique == 0) {
            restore_block(dg, block);
        }

        block_write_file_formatted(txt_file, block);

        // non-unique blocks allocated and restored previously must be de-allocated
        if (block->is_unique == 0) {
            destroy_block(block);
        }
    }

    fclose(txt_file);
}


/*********************** helper functions ********************/

// locally used functions

typedef struct {
    int a;
    int b;
} pair_t;


static int cmp_pair_t(const void *op1, const void *op2)
{
    return ((pair_t *) op1)->b - ((pair_t *) op2)->b;
}


static void inverse_perm(int n, int *perm, int *out)
{
    pair_t tmp[CC_DIAGRAM_MAX_RANK];
    int i;
    for (i = 0; i < n; i++) {
        tmp[i].a = i + 1;
        tmp[i].b = perm[i];
    }
    qsort(tmp, n, sizeof(pair_t), cmp_pair_t);
    for (i = 0; i < n; i++) {
        out[i] = tmp[i].a;
    }
}


void parse_qp_string(int rank, char *qp_str, int *qp)
{
    for (int i = 0; i < rank; i++) {
        int c = qp_str[i];
        if (c == 'h' || c == 'p') {
            qp[i] = c;
        }
        else {
            errquit("wrong quasiparticle symbol: %c (allowed are: h,p)", c);
        }
    }
}


void parse_valence_string(int rank, char *val_str, int *val)
{
    for (int i = 0; i < rank; i++) {
        int c = val_str[i];
        if (c == '0' || c == '1') {
            val[i] = c - '0';
        }
        else {
            errquit("wrong valence/inactive flag: %c (allowed are: 0,1)", c);
        }
    }
}


void parse_order_string(int rank, char *order_str, int *order)
{
    for (int i = 0; i < rank; i++) {
        int c = order_str[i];
        if (isdigit(c)) {
            order[i] = c - '0';
        }
        else {
            errquit("wrong order symbol: %c (only digits are allowed)", c);
        }
    }
}


int guess_rank(char *qparts, char *valence, char *order)
{
    int rk1 = strlen(qparts);
    int rk2 = strlen(valence);
    int rk3 = strlen(order);

    if (rk1 != rk2 || rk2 != rk3 || rk1 != rk3) {
        printf("strlen(qparts='%s') = %d\n", qparts, rk1);
        printf("strlen(valence='%s') = %d\n", valence, rk2);
        printf("strlen(order='%s') = %d\n", order, rk3);
        errquit("lengths of arguments 'qparts', 'valence' and 'order' must "
                "coincide and be equal to 2, 4, 6 ...");
    }

    if (rk1 == 0 || rk1 % 2 == 1) {
        errquit("lengths of arguments 'qparts', 'valence' and 'order' must "
                "be equal to 2, 4, 6 ... (actual length == %d)", rk1);
    }

    return rk1;
}


// storage class: RAM or DISK
int guess_storage_class(char *name, int rank, char *qparts, char *valence)
{
    int storage_type = CC_DIAGRAM_IN_MEM;
    if (strcmp(name, "pppp") == 0 || strcmp(name, "ppppr") == 0) {
        if (cc_opts->disk_usage_level >= 2) {
            storage_type = CC_DIAGRAM_ON_DISK;
        }
        else {
            storage_type = CC_DIAGRAM_IN_MEM;
        }
    }
    else if (rank == 4 && cc_opts->disk_usage_level >= 3) { // three inactive 'p' indices
        int np = 0;
        for (int i = 0; i < rank; i++) {
            if (qparts[i] == 'p' && valence[i] == '0') {
                np++;
            }
        }
        if (np >= 3) {
            storage_type = CC_DIAGRAM_ON_DISK;
        }
    }
    else if (rank >= 6 && cc_opts->disk_usage_level >= 1) { // triples+ diagrams: always on disk
        storage_type = CC_DIAGRAM_ON_DISK;
    }
    else {
        storage_type = CC_DIAGRAM_IN_MEM;
    }
    return storage_type;
}
