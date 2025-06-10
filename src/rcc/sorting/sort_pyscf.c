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


#include "sort.h"

#include <assert.h>
#include <complex.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sorting_request.h"

#include "io.h"
#include "memory.h"
#include "options.h"
#include "spinors.h"
#include "timer.h"

enum {
    CC_DIRECT,
    CC_EXCHANGE
};


double complex ****allocate_twoel_buffer(size_t dim);

void free_twoel_buffer(size_t dim, double complex ****v_ints);

size_t load_coulomb_block(double complex ****v_ints,
                          int spinor_block_1, int spinor_block_2, int spinor_block_3, int spinor_block_4);

size_t load_gaunt_block(double complex ****v_ints,
                        int spinor_block_1, int spinor_block_2, int spinor_block_3, int spinor_block_4);

size_t load_twoprop_blocks(double complex ****v_ints,
                           int spinor_block_1, int spinor_block_2, int spinor_block_3, int spinor_block_4);

void clear_twoel_buffer(double complex ****v_ints, int spinor_block_1, int spinor_block_2, int spinor_block_3,
                        int spinor_block_4);

void fill_block_twoelec(block_t *block, int sign, double complex ****v_ints, int direct);

void fill_block_twoelec_complex_direct(block_t *block, int sign, double complex ****v_ints);

void fill_block_twoelec_complex_exchange(block_t *block, int sign, double complex ****v_ints);

void fill_block_twoelec_real_direct(block_t *block, int sign, double complex ****v_ints);

void fill_block_twoelec_real_exchange(block_t *block, int sign, double complex ****v_ints);

int is_pppp_diagram(char *holes_particles, char *valence, char *order);

double complex get_eri(int nspinors, double complex *eris, size_t *eri_index, int p, int q, int r, int s);

double complex get_eri_goldstone(int nspinors, double complex *eris, int p, int q, int r, int s);

void fill_block_twoelec_complex_pyscf(block_t *block);

extern double complex *pyscf_unique_eri;
extern size_t *pyscf_eri_index;


void fill_block_one_elec(block_t *block, double complex *ints_matrix, int ignore_diagonal);


void sort_pyscf_one_electron()
{
    int nspinors = get_num_spinors();
    double complex *one_electron_ints = cc_calloc(nspinors * nspinors, sizeof(double complex));

    // only for debug purposes:
    /*for (int i = 0; i < nspinors; i++) {
        for (int j = 0; j < nspinors; j++) {
            if (i != j) {
                one_electron_ints[i * nspinors + j] = 1e-2;
            }
        }
    }*/


    for (int ireq = 0; ireq < n_requests; ireq++) {
        sorting_request_t *req = &sorting_requests[ireq];

        diagram_t *dg = req->dg;
        if (dg->rank != 2) {
            continue;
        }

        int spinor_blocks_nums[2];
        spinor_blocks_nums[0] = 0;
        spinor_blocks_nums[1] = 0;

        block_t *sb = diagram_get_block(dg, spinor_blocks_nums);
        if (sb == NULL) {
            continue;
        }
        if (sb->is_unique == 0) {
            continue;
        }

        // block was found. fill it!
        block_load(sb);

        fill_block_one_elec(sb, one_electron_ints, 1);

        block_store(sb);
    } // end of loop over requests

    cc_free(one_electron_ints);
}


void sort_pyscf_two_electron()
{
    if (num_twoelec_requests(sorting_requests, n_requests) == 0) {
        // nothing to sort
        return;
    }

    for (int ireq = 0; ireq < n_requests; ireq++) {
        sorting_request_t *req = &sorting_requests[ireq];

        diagram_t *dg = req->dg;
        if (dg->rank != 4) { // one-electron diagrams will be sorted later
            continue;
        }

        if (strcmp(req->dg_name, "hhhh") == 0 ||
            strcmp(req->dg_name, "hhpp") == 0 ||
            strcmp(req->dg_name, "pphh") == 0 ||
            strcmp(req->dg_name, "phhp") == 0 ||
            strcmp(req->dg_name, "pphp") == 0 ||
            strcmp(req->dg_name, "phpp") == 0 ||
            strcmp(req->dg_name, "phhh") == 0 ||
            strcmp(req->dg_name, "hhhp") == 0 ||
            strcmp(req->dg_name, "hphh") == 0 ||
            strcmp(req->dg_name, "hphp") == 0 ||
            strcmp(req->dg_name, "pppxpr") == 0) {

            block_t *block = dg->blocks[0];
            block_load(block);

            size_t dims_1 = block->shape[0];
            size_t dims_2 = block->shape[1];
            size_t dims_3 = block->shape[2];
            size_t dims_4 = block->shape[3];
            printf(" reading %s...\n", req->dg_name);

            char file_name[256];
            sprintf(file_name, "../pyscf_%s.dat", req->dg_name);

            int f_hhhh = io_open(file_name, "r");
            io_read(f_hhhh, block->buf, sizeof(double complex) * dims_1 * dims_2 * dims_3 * dims_4);
            io_close(f_hhhh);

            printf(" end %s\n", req->dg_name);

            block_store(block);

            continue;
        }

        if (strcmp(req->dg_name, "ppppr") == 0) {

            printf("begin pppp\n");

            block_t *block = dg->blocks[0];
            size_t dims_1 = block->shape[0];
            size_t dims_2 = block->shape[1];
            size_t dims_3 = block->shape[2];
            size_t dims_4 = block->shape[3];
            size_t num_part = dims_1;

            size_t coef_1 = num_part * num_part * num_part;
            size_t coef_2 = num_part * num_part;
            size_t coef_3 = num_part;

            const size_t buf_size = 100000000;
            double complex *buf = cc_calloc(buf_size, sizeof(double complex));

            int f_pppp = io_open("../pyscf_pppp.dat", "r");

            //size_t index = 0;
            size_t total_size = dims_1 * dims_2 * dims_3 * dims_4;
            size_t num_remain = total_size;
            size_t index = 0;

            while (num_remain > 0) {

                size_t num_read = (num_remain < buf_size) ? num_remain : buf_size;
                io_read(f_pppp, buf, sizeof(double complex) * num_read);
                num_remain -= num_read;

                printf("num_read = %ld\n", num_read);

                //index += num_read;

                for (size_t i = 0; i < num_read; i++) {

                    // linear index -> 4-index
                    size_t linear_index = index;

                    size_t idx_1 = linear_index / coef_1;
                    linear_index = linear_index % coef_1;

                    size_t idx_2 = linear_index / coef_2;
                    linear_index = linear_index % coef_2;

                    size_t idx_3 = linear_index / coef_3;
                    linear_index = linear_index % coef_3;

                    size_t idx_4 = linear_index;

                    //printf("%ld -> %ld %ld %ld %ld\n", index, idx_1, idx_2, idx_3, idx_4);

                    block->buf[coef_1 * idx_1 + coef_2 * idx_3 + coef_3 * idx_2 + idx_4] = conj(buf[i]);

                    // 4-index -> linear index

                    index++;
                }
            }

            io_close(f_pppp);

            cc_free(buf);

            printf("end pppp\n");
            continue;
        }

        /*int spinor_blocks_nums[4];
        spinor_blocks_nums[0] = 0;
        spinor_blocks_nums[1] = 0;
        spinor_blocks_nums[2] = 0;
        spinor_blocks_nums[3] = 0;

        block_t *sb = diagram_get_block(dg, spinor_blocks_nums);
        if (sb == NULL) {
            continue;
        }
        if (sb->is_unique == 0) {
            continue;
        }

        // block was found. fill it!
        block_load(sb);

        fill_block_twoelec_complex_pyscf(sb);

        if (is_pppp_diagram(req->hp, req->valence, req->order) && (arith == CC_ARITH_COMPLEX)) {
            conj_vector(sb->size, sb->buf);
        }

        //fill_block_twoelec(sb, 1, v_ints, CC_DIRECT);
        block_store(sb);*/
    } // end of loop over requests
}


void fill_block_twoelec_complex_pyscf(block_t *block)
{
    assert(block->rank == 4);

    int nspinors = get_num_spinors();

    int dims_1 = block->shape[0];
    int dims_2 = block->shape[1];
    int dims_3 = block->shape[2];
    int dims_4 = block->shape[3];
    int *block_indices_1 = block->indices[0];
    int *block_indices_2 = block->indices[1];
    int *block_indices_3 = block->indices[2];
    int *block_indices_4 = block->indices[3];

    size_t index = 0;
    for (int i1 = 0; i1 < dims_1; i1++) {
        int idx_1 = block_indices_1[i1]; // relative index to spinor indices
        //idx_1 = spinor_index_global2local[idx_1];
        for (int i2 = 0; i2 < dims_2; i2++) {
            int idx_2 = block_indices_2[i2];
            //idx_2 = spinor_index_global2local[idx_2];
            for (int i3 = 0; i3 < dims_3; i3++) {
                int idx_3 = block_indices_3[i3];
                //idx_3 = spinor_index_global2local[idx_3];
                for (int i4 = 0; i4 < dims_4; i4++) {
                    int idx_4 = block_indices_4[i4];
                    //idx_4 = spinor_index_global2local[idx_4];

                    double complex v_1234 = 0.0 + 0.0 * I;
                    if (cc_opts->use_goldstone) {
                        v_1234 = get_eri_goldstone(nspinors, pyscf_unique_eri, idx_1, idx_2, idx_3, idx_4);
                    }
                    else {
                        v_1234 = get_eri(nspinors, pyscf_unique_eri, pyscf_eri_index, idx_1, idx_2, idx_3, idx_4);
                    }

                    block->buf[index] = v_1234;
                    index++;
                }
            }
        }
    }
}


void pyscf_data_free()
{
    if (pyscf_eri_index) {
        cc_free(pyscf_eri_index);
        pyscf_eri_index = NULL;
    }

    if (pyscf_unique_eri) {
        cc_free(pyscf_unique_eri);
        pyscf_unique_eri = NULL;
    }
}
