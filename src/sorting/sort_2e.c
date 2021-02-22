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

#include "sort.h"

#include <assert.h>
#include <complex.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "sorting_request.h"

#include "datamodel.h"
#include "io.h"
#include "memory.h"
#include "options.h"
#include "spinors.h"
#include "timer.h"

enum {
    CC_DIRECT,
    CC_EXCHANGE
};

extern int32_t *symblock_indices;

double complex ****allocate_twoel_buffer(size_t dim);

void free_twoel_buffer(size_t dim, double complex ****v_ints);

int load_coulomb_block(double complex ****v_ints, int spinor_block_1, int spinor_block_2, int spinor_block_3,
                       int spinor_block_4);

void clear_twoel_buffer(double complex ****v_ints, int spinor_block_1, int spinor_block_2, int spinor_block_3,
                        int spinor_block_4);

void fill_block_twoelec(block_t *block, int sign, double complex ****v_ints, int direct);
void fill_block_twoelec_complex_direct(block_t *block, int sign, double complex ****v_ints);
void fill_block_twoelec_complex_exchange(block_t *block, int sign, double complex ****v_ints);
void fill_block_twoelec_real_direct(block_t *block, int sign, double complex ****v_ints);
void fill_block_twoelec_real_exchange(block_t *block, int sign, double complex ****v_ints);

int is_pppp_diagram(char *holes_particles, char *valence, char *order);


void sort_twoel()
{
    int spinor_block_1, spinor_block_2, spinor_block_3, spinor_block_4;
    size_t n_integrals_read_total = 0;
    size_t n_blocks_processed_total = 0;
    sorting_request_t *req;
    double complex ****v_ints = NULL;
    const double BYTES_TO_GB = 1.0 / (1024.0 * 1024.0 * 1024.0);

    if (num_twoelec_requests(sorting_requests, n_requests) == 0) {
        // nothing to sort
        return;
    }

    printf("   sorting two-electron integrals\n");
    printf("     step    #blocks   ints read    time,s  rate,G/s\n");
    double sort_twoel_time_start = abs_time();

    v_ints = allocate_twoel_buffer(get_max_spinor_block_size());

    for (spinor_block_1 = 0; spinor_block_1 < n_spinor_blocks; spinor_block_1++) {
        int n_blocks_processed = 0;
        int n_integrals_read = 0;
        double time_start = abs_time();

        printf("   %3d /%3d ", spinor_block_1, n_spinor_blocks);

        for (spinor_block_2 = 0; spinor_block_2 < n_spinor_blocks; spinor_block_2++) {
            for (spinor_block_3 = 0; spinor_block_3 < n_spinor_blocks; spinor_block_3++) {
                for (spinor_block_4 = 0; spinor_block_4 < n_spinor_blocks; spinor_block_4++) {

                    clear_twoel_buffer(v_ints, spinor_block_1, spinor_block_2, spinor_block_3, spinor_block_4);

                    int n_ints = load_coulomb_block(v_ints, spinor_block_1, spinor_block_2, spinor_block_3,
                                                          spinor_block_4);
                    if (n_ints == 0) {
                        continue;
                    }

                    // direct contribution
                    for (int ireq = 0; ireq < n_requests; ireq++) {
                        req = &sorting_requests[ireq];
                        diagram_t *dg = req->dg;
                        if (dg->rank != 4) { // one-electron diagrams will be sorted later
                            continue;
                        }

                        int spinor_blocks_nums[4];
                        spinor_blocks_nums[0] = spinor_block_1;
                        spinor_blocks_nums[1] = spinor_block_2;
                        spinor_blocks_nums[2] = spinor_block_3;
                        spinor_blocks_nums[3] = spinor_block_4;
                        size_t block_index;
                        block_t *sb = diagram_get_block(dg, spinor_blocks_nums, &block_index);
                        if (sb == NULL) {
                            continue;
                        }
                        if (sb->is_unique == 0) {
                            continue;
                        }

                        // block was found. fill it!
                        symblock_load(sb);
                        fill_block_twoelec(sb, 1, v_ints, CC_DIRECT);
                        symblock_store(sb);
                    } // end of loop over requests

                    // exchange contribution
                    for (int ireq = 0; ireq < n_requests; ireq++) {
                        req = &sorting_requests[ireq];
                        diagram_t *dg = req->dg;
                        if (dg->rank != 4) { // one-electron diagrams will be sorted later
                            continue;
                        }

                        int spinor_blocks_nums[4];
                        spinor_blocks_nums[0] = spinor_block_1;
                        spinor_blocks_nums[1] = spinor_block_2;
                        spinor_blocks_nums[2] = spinor_block_4;
                        spinor_blocks_nums[3] = spinor_block_3;
                        size_t block_index;
                        block_t *sb = diagram_get_block(dg, spinor_blocks_nums, &block_index);
                        if (sb == NULL) {
                            continue;
                        }
                        if (sb->is_unique == 0) {
                            continue;
                        }

                        // block was found. fill it!
                        symblock_load(sb);
                        fill_block_twoelec(sb, -1, v_ints, CC_EXCHANGE);
                        if (is_pppp_diagram(req->hp, req->valence, req->order) && carith) {
                            conj_vector(sb->size, sb->buf);
                        }
                        symblock_store(sb);

                        // update statistics
                        n_integrals_read += n_ints;
                        n_integrals_read_total += n_ints;
                        n_blocks_processed += 1;
                        n_blocks_processed_total += 1;
                    } // end of loop over requests
                }
            }
        }
        // intermediate statistics
        int n_bytes_read = n_integrals_read * (carith ? sizeof(double complex) : sizeof(double));
        double time_elapsed = abs_time() - time_start;
        printf("%8d", n_blocks_processed);
        printf("%12d", n_integrals_read);
        printf("%10.0f", time_elapsed);
        printf("%10.2f", n_bytes_read / time_elapsed * BYTES_TO_GB);
        printf("\n");
    } // end of 4 nested 'for' loops over spinor blocks

    free_twoel_buffer(get_max_spinor_block_size(), v_ints);

    // print statistics
    double sort_twoel_time_elapsed = abs_time() - sort_twoel_time_start;
    size_t n_bytes_read_total = n_integrals_read_total * (carith ? sizeof(double complex) : sizeof(double));
    printf("     total  ");
    printf("%8d", n_blocks_processed_total);
    printf("%12d", n_integrals_read_total);
    printf("%10.0f", sort_twoel_time_elapsed);
    printf("%10.2f", n_bytes_read_total / sort_twoel_time_elapsed * BYTES_TO_GB);
    printf("\n");
}


/**
 * allocates 4-dimensional array for two electron integrals of size
 * dim x dim x dim x dim.
 * returns pointer to the allocated array.
 */
double complex ****allocate_twoel_buffer(size_t dim)
{
    double complex ****buf = NULL;

    buf = (double complex ****) cc_malloc(sizeof(double complex ***) * dim);
    for (size_t i = 0; i < dim; i++) {
        buf[i] = (double complex ***) cc_malloc(sizeof(double complex **) * dim);
        for (size_t j = 0; j < dim; j++) {
            buf[i][j] = (double complex **) cc_malloc(sizeof(double complex *) * dim);
            for (size_t k = 0; k < dim; k++) {
                buf[i][j][k] = (double complex *) cc_malloc(sizeof(double complex) * dim);
            }
        }
    }

    return buf;
}


int num_twoelec_requests(sorting_request_t *requests, int num_requests)
{
    // count numbers of 1-el and 2-el diagrams to be sorted
    int n_requests_2el = 0;

    for (int ireq = 0; ireq < num_requests; ireq++) {
        sorting_request_t *req = requests + ireq;
        if (req->dg->rank == 4) {
            n_requests_2el++;
        }
    }

    return n_requests_2el;
}


/**
 * Returns number of integrals read
 * vint_array = factor1 * vint_array + factor2 * new_integrals
 */
int read_twoel_block_unformatted(char *vint_file_name, double complex ****vint_array, double complex factor1,
                                 double complex factor2, size_t *n_bytes_read)
{
    int32_t nint;
    int32_t iint;
    int i, j, k, l;
    size_t n_integrals_read = 0;
    size_t nbr = 0;
    int fd;
    int nloop = 0;
    int empty_file = 0;
    double dfactor1 = creal(factor1);
    double dfactor2 = creal(factor2);

    static double complex buf_integrals[CC_SORTING_IO_BUF_SIZE];
    static int16_t buf_indices[4 * CC_SORTING_IO_BUF_SIZE];

    fd = io_open(vint_file_name, "r");

    while (1) {
        io_read_compressed(fd, &nint, sizeof(int32_t));
        nloop++;

        if (nint == 0) {
            if (nloop == 1) { empty_file = 1; }
            break;
        }

        io_read_compressed(fd, buf_indices, nint * 4 * sizeof(int16_t));
        if (carith) {
            io_read_compressed(fd, buf_integrals, nint * sizeof(double complex));
        }
        else {
            io_read_compressed(fd, buf_integrals, nint * sizeof(double));
        }
        for (iint = 0; iint < nint; iint++) {
            // absolute spinor indices
            i = buf_indices[iint * 4] - 1; // -1 since numeration from 0 in C and from 1 in F90
            j = buf_indices[iint * 4 + 1] - 1;
            k = buf_indices[iint * 4 + 2] - 1;
            l = buf_indices[iint * 4 + 3] - 1;
            i = spinor_index_global2local[i];
            j = spinor_index_global2local[j];
            k = spinor_index_global2local[k];
            l = spinor_index_global2local[l];

            // here "local" spinor indices are used
            // (inside the spinor block under consideration)
            if (carith) {
                double complex v_ijkl = buf_integrals[iint];
                vint_array[i][j][k][l] = factor1 * vint_array[i][j][k][l] + factor2 * v_ijkl;
            }
            else {
                double v_ijkl = ((double *) buf_integrals)[iint];
                vint_array[i][j][k][l] = dfactor1 * vint_array[i][j][k][l] + dfactor2 * v_ijkl + 0.0 * I;
            }
        }
        n_integrals_read += nint;
        if (arith == CC_ARITH_REAL) {
            nbr += nint * sizeof(double) + nint * 4 * sizeof(int16_t) + sizeof(int32_t);
        }
        else {  // complex arith
            nbr += nint * sizeof(double complex) + nint * 4 * sizeof(int16_t) + sizeof(int32_t);
        }
    }
    io_close(fd);
    *n_bytes_read = nbr;

    if (empty_file) {
        return 0;
    }

    return n_integrals_read;
}


void clear_twoel_buffer(double complex ****v_ints, int spinor_block_1, int spinor_block_2, int spinor_block_3,
                        int spinor_block_4)
{
    // sizes of spinor blocks
    int block1_size = get_spinor_block_size(spinor_block_1);
    int block2_size = get_spinor_block_size(spinor_block_2);
    int block3_size = get_spinor_block_size(spinor_block_3);
    int block4_size = get_spinor_block_size(spinor_block_4);

    // init buffer array with zeroes
    for (int i = 0; i < block1_size; i++) {
        for (int j = 0; j < block2_size; j++) {
            for (int k = 0; k < block3_size; k++) {
                memset(v_ints[i][j][k], 0, block4_size * sizeof(double complex));
            }
        }
    }
}


int load_coulomb_block(double complex ****v_ints, int spinor_block_1, int spinor_block_2, int spinor_block_3,
                       int spinor_block_4)
{
    char vint_file_name[CC_MAX_FILE_NAME_LENGTH];
    int n_integrals_read;
    size_t n_bytes_read;

    sprintf(vint_file_name, "VINT-%d-%d-%d-%d",
            spinor_block_1 + 1, spinor_block_2 + 1, spinor_block_3 + 1, spinor_block_4 + 1);
    if (io_file_exists(vint_file_name) == 0) {
        return 0;
    }

    n_integrals_read = read_twoel_block_unformatted(vint_file_name, v_ints, 0.0 + 0.0 * I,
                                                    1.0 + 0.0 * I, &n_bytes_read);

    return n_integrals_read;
}


void fill_block_twoelec(block_t *block, int sign, double complex ****v_ints, int direct_flag)
{
    assert(block->rank == 4);
    assert(direct_flag == CC_DIRECT || direct_flag == CC_EXCHANGE);

    // four different functions for maximum performance
    if (carith && direct_flag == CC_DIRECT) {
        fill_block_twoelec_complex_direct(block, sign, v_ints);
    }
    else if (carith && direct_flag == CC_EXCHANGE) {
        fill_block_twoelec_complex_exchange(block, sign, v_ints);
    }
    else if ((!carith) && direct_flag == CC_DIRECT) {
        fill_block_twoelec_real_direct(block, sign, v_ints);
    }
    else if ((!carith) && direct_flag == CC_EXCHANGE) {
        fill_block_twoelec_real_exchange(block, sign, v_ints);
    }
}


void fill_block_twoelec_complex_direct(block_t *block, int sign, double complex ****v_ints)
{
    assert(block->rank == 4);

    int dims_1 = block->indices[0][0];
    int dims_2 = block->indices[1][0];
    int dims_3 = block->indices[2][0];
    int dims_4 = block->indices[3][0];
    int *block_indices_1 = block->indices[0] + 1; // +1 to skip first element (length)
    int *block_indices_2 = block->indices[1] + 1;
    int *block_indices_3 = block->indices[2] + 1;
    int *block_indices_4 = block->indices[3] + 1;

    size_t index = 0;
    for (int i1 = 0; i1 < dims_1; i1++) {
        int idx_1 = block_indices_1[i1]; // relative index to spinor indices
        idx_1 = spinor_index_global2local[idx_1];
        for (int i2 = 0; i2 < dims_2; i2++) {
            int idx_2 = block_indices_2[i2];
            idx_2 = spinor_index_global2local[idx_2];
            for (int i3 = 0; i3 < dims_3; i3++) {
                int idx_3 = block_indices_3[i3];
                idx_3 = spinor_index_global2local[idx_3];
                for (int i4 = 0; i4 < dims_4; i4++) {
                    int idx_4 = block_indices_4[i4];
                    idx_4 = spinor_index_global2local[idx_4];

                    block->buf[index] += sign * v_ints[idx_1][idx_2][idx_3][idx_4];
                    index++;
                }
            }
        }
    }
}


void fill_block_twoelec_complex_exchange(block_t *block, int sign, double complex ****v_ints)
{
    assert(block->rank == 4);

    int dims_1 = block->indices[0][0];
    int dims_2 = block->indices[1][0];
    int dims_3 = block->indices[2][0];
    int dims_4 = block->indices[3][0];
    int *block_indices_1 = block->indices[0] + 1; // +1 to skip first element (length)
    int *block_indices_2 = block->indices[1] + 1;
    int *block_indices_3 = block->indices[2] + 1;
    int *block_indices_4 = block->indices[3] + 1;

    size_t index = 0;
    for (int i1 = 0; i1 < dims_1; i1++) {
        int idx_1 = block_indices_1[i1]; // relative index to spinor indices
        idx_1 = spinor_index_global2local[idx_1];
        for (int i2 = 0; i2 < dims_2; i2++) {
            int idx_2 = block_indices_2[i2];
            idx_2 = spinor_index_global2local[idx_2];
            for (int i3 = 0; i3 < dims_3; i3++) {
                int idx_3 = block_indices_3[i3];
                idx_3 = spinor_index_global2local[idx_3];
                for (int i4 = 0; i4 < dims_4; i4++) {
                    int idx_4 = block_indices_4[i4];
                    idx_4 = spinor_index_global2local[idx_4];

                    block->buf[index] += sign * v_ints[idx_1][idx_2][idx_4][idx_3];
                    index++;
                }
            }
        }
    }
}


void fill_block_twoelec_real_direct(block_t *block, int sign, double complex ****v_ints)
{
    assert(block->rank == 4);

    double *dbuf = (double *) block->buf;

    int dims_1 = block->indices[0][0];
    int dims_2 = block->indices[1][0];
    int dims_3 = block->indices[2][0];
    int dims_4 = block->indices[3][0];
    int *block_indices_1 = block->indices[0] + 1; // +1 to skip first element (length)
    int *block_indices_2 = block->indices[1] + 1;
    int *block_indices_3 = block->indices[2] + 1;
    int *block_indices_4 = block->indices[3] + 1;

    size_t index = 0;
    for (int i1 = 0; i1 < dims_1; i1++) {
        int idx_1 = block_indices_1[i1]; // relative index to spinor indices
        idx_1 = spinor_index_global2local[idx_1];
        for (int i2 = 0; i2 < dims_2; i2++) {
            int idx_2 = block_indices_2[i2];
            idx_2 = spinor_index_global2local[idx_2];
            for (int i3 = 0; i3 < dims_3; i3++) {
                int idx_3 = block_indices_3[i3];
                idx_3 = spinor_index_global2local[idx_3];
                for (int i4 = 0; i4 < dims_4; i4++) {
                    int idx_4 = block_indices_4[i4];
                    idx_4 = spinor_index_global2local[idx_4];

                    dbuf[index] += sign * creal(v_ints[idx_1][idx_2][idx_3][idx_4]);
                    index++;
                }
            }
        }
    }
}


void fill_block_twoelec_real_exchange(block_t *block, int sign, double complex ****v_ints)
{
    assert(block->rank == 4);

    double *dbuf = (double *) block->buf;

    int dims_1 = block->indices[0][0];
    int dims_2 = block->indices[1][0];
    int dims_3 = block->indices[2][0];
    int dims_4 = block->indices[3][0];
    int *block_indices_1 = block->indices[0] + 1; // +1 to skip first element (length)
    int *block_indices_2 = block->indices[1] + 1;
    int *block_indices_3 = block->indices[2] + 1;
    int *block_indices_4 = block->indices[3] + 1;

    size_t index = 0;
    for (int i1 = 0; i1 < dims_1; i1++) {
        int idx_1 = block_indices_1[i1]; // relative index to spinor indices
        idx_1 = spinor_index_global2local[idx_1];
        for (int i2 = 0; i2 < dims_2; i2++) {
            int idx_2 = block_indices_2[i2];
            idx_2 = spinor_index_global2local[idx_2];
            for (int i3 = 0; i3 < dims_3; i3++) {
                int idx_3 = block_indices_3[i3];
                idx_3 = spinor_index_global2local[idx_3];
                for (int i4 = 0; i4 < dims_4; i4++) {
                    int idx_4 = block_indices_4[i4];
                    idx_4 = spinor_index_global2local[idx_4];

                    dbuf[index] += sign * creal(v_ints[idx_1][idx_2][idx_4][idx_3]);
                    index++;
                }
            }
        }
    }
}


void free_twoel_buffer(size_t dim, double complex ****v_ints)
{
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                cc_free(v_ints[i][j][k]);
            }
            cc_free(v_ints[i][j]);
        }
        cc_free(v_ints[i]);
    }
    cc_free(v_ints);
}


int is_pppp_diagram(char *holes_particles, char *valence, char *order)
{
    if (strcmp(holes_particles, "pppp") == 0 &&
        strcmp(valence, "0000") == 0 &&
        strcmp(order, "3412") == 0) {
        return 1;
    }
    return 0;
}