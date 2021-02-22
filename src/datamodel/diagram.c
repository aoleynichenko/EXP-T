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
 * diagram.c
 *
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
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io.h"
#include "datamodel.h"
#include "error.h"
#include "memory.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"

static int cmp_pair_t(const void *op1, const void *op2);
static void inverse_perm(int n, int *perm, int *out);

// counter of diagrams (for unique ID's)
int64_t diagrams_count = 0;


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
diagram_t *diagram_new(char *name, char *qparts, char *valence, char *order, int perm_unique)
{
    int rk1, rk2, rk3, rank;
    int qparts_arr[CC_DIAGRAM_MAX_RANK];
    int valence_arr[CC_DIAGRAM_MAX_RANK];
    int order_arr[CC_DIAGRAM_MAX_RANK];
    int *spinor_blocks_list;
    block_t *sb;
    size_t i;
    char c;
    block_t **sbs_temp;  // temporary storage for the created symmetry blocks
    int max_sbs;            // max number of sym blocks = (# spinor blocks)^rank
    int isymb = 0;          // counter of created blocks
    int *ijkl;              // indices for the arbitrary-nested loop (len = rank)
    int only_unique = perm_unique;

    // check arguments for correctness and pre-process them
    rk1 = strlen(qparts);
    rk2 = strlen(valence);
    rk3 = strlen(order);
    if (rk1 != rk2 || rk2 != rk3 || rk1 != rk3) {
        printf("name=%s\n", name);
        printf("strlen(qparts='%s') = %d\n", qparts, rk1);
        printf("strlen(valence='%s') = %d\n", valence, rk2);
        printf("strlen(order='%s') = %d\n", order, rk3);
        errquit("lengths of arguments 'qparts', 'valence' and 'order' must "
                "coincide and be equal to 2, 4, 6 ...");
    }
    rank = rk1;
    if (rank == 0 || rank % 2 == 1) {
        errquit("lengths of arguments 'qparts', 'valence' and 'order' must "
                "be equal to 2, 4, 6 ... (actual length == %d)", rank);
    }
    // strings -> integer arrays
    for (i = 0; i < rank; i++) {
        c = qparts[i];
        if (c == 'h' || c == 'p') {
            qparts_arr[i] = c;
        }
        else {
            errquit("wrong quasiparticle symbol: %c (allowed are: h,p)", c);
        }
        c = valence[i];
        if (c == '0' || c == '1') {
            valence_arr[i] = c - '0';
        }
        else {
            errquit("wrong valence/inactive flag: %c (allowed are: 0,1)", c);
        }
        // TODO: careful check of the 'order' argument here
        c = order[i];
        if (isdigit(c)) {
            order_arr[i] = c - '0';
        }
        else {
            errquit("wrong order symbol: %c (only digits are allowed)", c);
        }
    }

    // if the diagram's name starts with '$', only permutationally unique blocks will be treated
    if (name[0] == '$') {
        only_unique = 1;
    }

    // allocate memory for all the metainfo and the diagram itself
    diagram_t *dg = (diagram_t *) cc_malloc(1 * sizeof(diagram_t));

    // fill diagram's fields
    dg->dg_id = diagrams_count++;
    strncpy(dg->name, name, CC_DIAGRAM_MAX_NAME);
    dg->name[CC_DIAGRAM_MAX_NAME - 1] = '\0';
    dg->rank = rank;
    dg->only_unique = only_unique;
    intcpy(dg->qparts, qparts_arr, rank);
    intcpy(dg->valence, valence_arr, rank);
    intcpy(dg->order, order_arr, rank);

    // create inverted index (numbers of spinor blocks -> number of block)
    size_t inv_index_size = int_pow(n_spinor_blocks, rank);
    dg->inv_index = (size_t *) cc_malloc(inv_index_size * sizeof(size_t));
    memset(dg->inv_index, 0, inv_index_size * sizeof(size_t));

    // permutation which is inverse to the 'order' perm-n
    // is required for the subsequent application of the DPD scheme
    int reverse_order[CC_DIAGRAM_MAX_RANK];
    inverse_perm(rank, order_arr, reverse_order);

    // storage class: RAM or DISK
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

    // create symmetry blocks
    max_sbs = (int) pow(n_spinor_blocks, rank);
    sbs_temp = (block_t **) cc_malloc(sizeof(block_t *) * max_sbs);
    spinor_blocks_list = (int *) cc_malloc(sizeof(int) * rank);

    // dynamically-nested loop over tensor's dimensions
    ijkl = (int *) cc_malloc(sizeof(int) * rank);
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

        // check if this block is symmetry-allowed
        if (!dpd_contains_totsym_rank[rank / 2 - 1](sym)) {
            goto next_symblock;
        }

        int is_zero = is_symblock_zero(rank, ijkl, dg->qparts, dg->valence);
        if (is_zero) {
            goto next_symblock;
        }

        // create empty block
        sb = symblock_new(rank, ijkl, dg->qparts, dg->valence, dg->order, storage_type, only_unique);

        if (sb == NULL) {
            goto next_symblock;
        }

        // assign fields
        sb->dg_id = dg->dg_id;
        sb->sb_id = isymb;

        // save this block (and count it)
        sbs_temp[isymb] = sb;
        isymb++;

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

    int dims[CC_DIAGRAM_MAX_RANK];
    for (i = 0; i < CC_DIAGRAM_MAX_RANK; i++) {
        dims[i] = n_spinor_blocks;
    }

    // bind non-zero blocks to the diagram
    dg->n_blocks = isymb;
    dg->blocks = (block_t **) cc_malloc(sizeof(block_t *) * isymb);
    for (i = 0; i < dg->n_blocks; i++) {
        dg->blocks[i] = sbs_temp[i];
        dg->blocks[i]->dg_id = dg->dg_id;
        dg->blocks[i]->sb_id = i;

        // add block to the inverted index
        size_t ii = as_linear_index(rank, dims, dg->blocks[i]->spinor_blocks);
        dg->inv_index[ii] = i;
    }

    // cleanup
    cc_free(sbs_temp);
    cc_free(spinor_blocks_list);
    cc_free(ijkl);

    return dg;
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
 * get number of symblock (its index)
 * using inverted index
 *
 * BUG: returns wrong pointer! to the moment only the block_index parameter
 * can be used safely!
 */
block_t *diagram_get_block(diagram_t *dg, int *spinor_blocks_nums, size_t *block_index)
{
    static int dims[CC_DIAGRAM_MAX_RANK];
    static int initialized = 0;

    if (initialized == 0) {
        initialized = 1;
        for (int i = 0; i < CC_DIAGRAM_MAX_RANK; i++) {
            dims[i] = n_spinor_blocks;
        }
    }

    size_t ii = as_linear_index(dg->rank, dims, spinor_blocks_nums);
    *block_index = dg->inv_index[ii];
    block_t *ret = dg->blocks[*block_index];

    for (int i = 0; i < dg->rank; i++) {
        if (ret->spinor_blocks[i] != spinor_blocks_nums[i]) {
            return NULL;
        }
    }

    return ret;
}


/*******************************************************************************
 * diagram_delete
 *
 * destructor -- deallocates all memory associated with this object.
 ******************************************************************************/
void diagram_delete(diagram_t *dg)
{
    for (size_t i = 0; i < dg->n_blocks; i++) {
        symblock_delete(dg->blocks[i]);
    }
    cc_free(dg->blocks);
    cc_free(dg->inv_index);
    cc_free(dg);
}


/*******************************************************************************
 * diagram_copy
 *
 * "Copy constructor".
 * Creates a copy of a diagram 'dg' (the same h/p, valence, order, indices, data
 * BUT the other IDs, names of files with integrals, etc)
 ******************************************************************************/
diagram_t *diagram_copy(diagram_t *dg)
{
    diagram_t *clone;
    char qparts[CC_DIAGRAM_MAX_RANK];
    char valence[CC_DIAGRAM_MAX_RANK];
    char order[CC_DIAGRAM_MAX_RANK];

    diagram_get_quasiparticles(dg, qparts);
    diagram_get_valence(dg, valence);
    diagram_get_order(dg, order);

    clone = diagram_new(dg->name, qparts, valence, order, dg->only_unique);

    // name of the new diagram consist of the old one + "_copy_" + clone's ID
    sprintf(clone->name, "%s_copy_%ld", dg->name, clone->dg_id);

    // copy data: block-by-block
    for (size_t i = 0; i < dg->n_blocks; i++) {
        symblock_copy_data(clone->blocks[i], dg->blocks[i]);
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


/*******************************************************************************
 * diagram_get
 *
 * TODO: error handling (only errquit is now implemented)
 ******************************************************************************/
double complex diagram_get(diagram_t *dg, int *idx)
{
    int rank = dg->rank;
    int sp_blocks[CC_DIAGRAM_MAX_RANK];
    int i, is;

    for (i = 0; i < rank; i++) {
        sp_blocks[i] = spinor_info[idx[i]].blockno;
    }

    // find required block
    for (is = 0; is < dg->n_blocks; is++) {
        // через spinor_block как сигнатуру
        // но это не настолько нужно
        if (intcmp(rank, dg->blocks[is]->spinor_blocks, sp_blocks) == 0) { // block found
            block_t *b = dg->blocks[is];
            if (b->is_unique == 0) {
                int idx_uniq[CC_DIAGRAM_MAX_RANK];
                int uniq_spinor_blocks[CC_DIAGRAM_MAX_RANK];

                transform(b->rank, b->spinor_blocks, uniq_spinor_blocks, b->perm_to_unique, 0);
                transform(b->rank, idx, idx_uniq, b->perm_to_unique, 0);

                size_t unique_index = 0;
                //block_t *uniq_block =
                diagram_get_block(dg, uniq_spinor_blocks, &unique_index);
                block_t *uniq_block = dg->blocks[unique_index];

                return b->sign * symblock_get(uniq_block, idx_uniq);
            }
            else {
                return symblock_get(b, idx);
            }
        }
    }

    return 0.0 + 0.0 * I;
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
    int i, is;

    // extract symmetry of all indices
    for (i = 0; i < rank; i++) {
        spinor_blocks[i] = spinor_info[idx[i]].blockno;
    }

    // find required block
    for (is = 0; is < dg->n_blocks; is++) {
        if (intcmp(rank, dg->blocks[is]->spinor_blocks, spinor_blocks) == 0) { // block found
            symblock_set(dg->blocks[is], val, idx);
            return;
        }
    }
}


/**
 * Prints diagram 'dg' to stdout (both metainfo & all integrals).
 * May cause extremely huge output; only for debug purposes.
 */
void diagram_print(diagram_t *dg)
{
    size_t i;

    printf("=================\n");
    printf("DIAGRAM %s\n", dg->name);
    printf("  id   = %ld\n", dg->dg_id);
    //printf("  next = %p\n", dg->next);
    printf("  rank = %d\n", dg->rank);
    printf("  quasiparticles: ");
    for (i = 0; i < dg->rank; i++) {
        printf("%c", dg->qparts[i]);
    }
    printf("\n");
    printf("  valence: ");
    for (i = 0; i < dg->rank; i++) {
        printf("%1d", dg->valence[i]);
    }
    printf("\n");
    printf("  order: ");
    for (i = 0; i < dg->rank; i++) {
        printf("%1d", dg->order[i]);
    }
    printf("\n");
    printf("number of non-zero blocks = %ld\n", dg->n_blocks);
    for (i = 0; i < dg->n_blocks; i++) {
        if (dg->blocks[i]->is_unique == 0) {
            printf("non unique, restoring...\n");
            restore_block(dg, dg->blocks[i]);
        }
        symblock_print(dg->blocks[i]);
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
    char str_order[CC_DIAGRAM_MAX_RANK + 1];
    size_t n_unique = 0;

    diagram_get_quasiparticles(diag, str_qparts);
    diagram_get_valence(diag, str_valence);
    diagram_get_order(diag, str_order);

    for (size_t i = 0; i < diag->n_blocks; i++) {
        block_t *b = diag->blocks[i];
        if (b->is_unique == 1) {
            n_unique++;
        }
    }

    printf(" diagram %s: %s %s %s %d/%d\n", diag->name, str_qparts, str_valence, str_order, n_unique, diag->n_blocks);
}


/**
 * Sets all matrix elements to zero.
 * NOTE: All diagrams must be allocated and initialized!
 */
diagram_t *diagram_clear(diagram_t *dg)
{
    assert(dg != NULL);

    for (size_t isb = 0; isb < dg->n_blocks; isb++) {
        symblock_clear(dg->blocks[isb]);
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
 * diagram_write
 *
 * Writes diagram to disk (binary format). The name of the file is 'dg->name.dg'
 * The platform-independent I/O interface is used (prefixes 'io_')
 * (see io.h and io.c for details)
 */
void diagram_write(diagram_t *dg, char *file_name)
{
    int f;  // file descriptor
    int i;
    block_t *sb;
    static const int32_t DG_MAGIC = 0x6f6c6579;

    f = io_open(file_name, "w");

    // magic number
    io_write(f, &DG_MAGIC, sizeof(DG_MAGIC));

    // diagram's ID
    io_write(f, &dg->dg_id, sizeof(dg->dg_id));

    // diagram's name (len = CC_DIAGRAM_MAX_NAME)
    io_write(f, dg->name, CC_DIAGRAM_MAX_NAME * sizeof(char));

    // diagram's rank
    io_write_compressed(f, &dg->rank, sizeof(dg->rank));

    // about uniqueness
    io_write_compressed(f, &dg->only_unique, sizeof(dg->only_unique));

    // quasiparticle type list (h/p), valence, order
    io_write(f, dg->qparts, sizeof(int) * dg->rank);
    io_write(f, dg->valence, sizeof(int) * dg->rank);
    io_write(f, dg->order, sizeof(int) * dg->rank);

    // inverted index
    io_write(f, dg->inv_index, sizeof(size_t) * int_pow(n_spinor_blocks, dg->rank));

    // total number of symmetry blocks in this diagram
    io_write(f, &dg->n_blocks, sizeof(dg->n_blocks));

    for (i = 0; i < dg->n_blocks; i++) {
        sb = dg->blocks[i];
        symblock_write(f, sb);
    }

    io_close(f);
}


/**
 * diagram_read
 *
 * Reads diagram from disk (from binary file).
 * The platform-independent I/O interface is used (prefixes 'io_')
 * (see io.h and io.c for details)
 */
diagram_t *diagram_read(char *file_name)
{
    int f;  // file descriptor
    int i;
    int32_t magic;
    diagram_t *dg;
    static const int32_t DG_MAGIC = 0x6f6c6579;

    f = io_open(file_name, "r");

    if (f == -1) {
        printf(" Unable to open diagram file %s\n", file_name);
        return NULL;
    }

    // magic number
    io_read(f, &magic, sizeof(int32_t));
    if (magic != DG_MAGIC) {
        printf("diagram_read: file_name = %s\n", file_name);
        printf("magic expected = %d\n", DG_MAGIC);
        printf("magic found = %d\n", magic);
        printf("wrong magic number! It is not a diagram file!\n");
        return NULL;
    }

    dg = (diagram_t *) cc_malloc(1 * sizeof(diagram_t));

    // diagram's ID
    io_read(f, &dg->dg_id, sizeof(dg->dg_id));
    // ID is new!
    dg->dg_id = diagrams_count++;

    // diagram's name (len = CC_DIAGRAM_MAX_NAME)
    io_read(f, dg->name, CC_DIAGRAM_MAX_NAME * sizeof(char));

    // diagram's rank
    io_read_compressed(f, &dg->rank, sizeof(dg->rank));

    // about uniqueness
    io_read_compressed(f, &dg->only_unique, sizeof(dg->only_unique));

    // quasiparticle type list (h/p), valence, order
    io_read(f, dg->qparts, sizeof(int) * dg->rank);
    io_read(f, dg->valence, sizeof(int) * dg->rank);
    io_read(f, dg->order, sizeof(int) * dg->rank);

    // inverted index
    dg->inv_index = (size_t *) cc_malloc(sizeof(size_t) * int_pow(n_spinor_blocks, dg->rank));
    io_read(f, dg->inv_index, sizeof(size_t) * int_pow(n_spinor_blocks, dg->rank));

    // total number of symmetry blocks in this diagram
    io_read(f, &dg->n_blocks, sizeof(dg->n_blocks));

    // alloc memory for blocks pointers & then read blocks
    dg->blocks = (block_t **) cc_malloc(dg->n_blocks * sizeof(block_t *));
    for (i = 0; i < dg->n_blocks; i++) {
        dg->blocks[i] = symblock_read(f);
        // set new diagram's unique ID
        dg->blocks[i]->dg_id = dg->dg_id;
    }

    io_close(f);

    // try to find in the stack diagram with the same name. if found -- replace it
    // with the new diagram, if not found -- add new diagram to the stack
    i = diagram_stack_find_index(dg->name);
    if (i == -1) {
        diagram_stack_push(dg);
    }
    else {
        diagram_stack_replace(dg->name, dg);
    }

    return dg;
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
