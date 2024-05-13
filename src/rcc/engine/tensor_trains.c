/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2024 The EXP-T developers.
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

#include "tensor_trains.h"

#ifdef TENSOR_TRAIN

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "engine.h"
#include "tt.h"
#include "timer.h"


int all_elements_zero(size_t n, void *buf, const double thresh);

double tensor_train_restoration_error(size_t n, void *original, void *restored);

double size_t_to_gygabytes(size_t num_bytes);


void tensor_decomposition_benchmark(char *diagram_name, double svd_threshold)
{
    printf("\n");
    printf(" tensor decomposition benchmark for diagram '%s'\n", diagram_name);

    diagram_t *diagram = diagram_stack_find(diagram_name);
    if (diagram == NULL) {
        printf(" diagram '%s' not found, nothing to be done\n\n", diagram_name);
        return;
    }

    double time_for_decomposition = 0.0;
    double time_for_restoration = 0.0;
    size_t mem_used_linear_arrays = 0;
    size_t mem_used_tensor_trains = 0;
    double max_restoration_error = 0.0;

    for (int iblock = 0; iblock < diagram->n_blocks; iblock++) {
        block_t *block = diagram->blocks[iblock];
        if (!block->is_unique) {
            continue;
        }
        block_load(block);

        if (all_elements_zero(block->size, block->buf, 1e-14)) {
            block_unload(block);
            continue;
        }

        // prepare tensor's shape as an array of type uint64_t
        int rank = diagram->rank;
        __uint64_t shape[rank];
        for (int i = 0; i < rank; i++) {
            shape[i] = block->shape[i];
        }

        if (arith == CC_ARITH_COMPLEX) {
            // decompose tensor
            double t_start = abs_time();
            void *tt_object = tt_z_new(block->buf, shape, rank, svd_threshold);
            time_for_decomposition += abs_time() - t_start;

            // save information about memory used by the tensor train
            __uint64_t *tt_shape = tt_z_shape(tt_object);
            size_t tt_size = 1;
            for (int i = 0; i < rank; i++) {
                tt_size += tt_shape[3 * i] * tt_shape[3 * i + 1] * tt_shape[3 * i + 2];
            }

            mem_used_linear_arrays += block->size;
            mem_used_tensor_trains += tt_size;

            // restore tensor
            t_start = abs_time();
            double complex *restored_tensor = tt_z_to_tensor_into(tt_object);
            time_for_restoration += abs_time() - t_start;

            double restor_error = tensor_train_restoration_error(block->size, block->buf, restored_tensor);
            if (restor_error > max_restoration_error) {
                max_restoration_error = restor_error;
            }

            // clean up
            free(tt_shape);
            free(restored_tensor);
            tt_free(tt_object);
        }
        else {
            // decompose tensor
            double t_start = abs_time();
            void *tt_object = tt_d_new((double *) block->buf, shape, rank, svd_threshold);
            time_for_decomposition += abs_time() - t_start;

            // save information about memory used by the tensor train
            __uint64_t *tt_shape = tt_d_shape(tt_object);
            size_t tt_size = 1;
            for (int i = 0; i < rank; i++) {
                tt_size += tt_shape[3 * i] * tt_shape[3 * i + 1] * tt_shape[3 * i + 2];
            }

            mem_used_linear_arrays += block->size;
            mem_used_tensor_trains += tt_size;

            // restore tensor
            t_start = abs_time();
            double *restored_tensor = tt_d_to_tensor_into(tt_object);
            time_for_restoration += abs_time() - t_start;

            double restor_error = tensor_train_restoration_error(block->size, block->buf, restored_tensor);
            if (restor_error > max_restoration_error) {
                max_restoration_error = restor_error;
            }

            // clean up
            free(tt_shape);
            free(restored_tensor);
            tt_free(tt_object);
        }

        block_unload(block);
    }

    printf("   svd threshold = %.3e\n", svd_threshold);
    printf("   memory used for linear arrays    %ld numbers = %.3f gb\n", mem_used_linear_arrays,
           size_t_to_gygabytes(mem_used_linear_arrays * SIZEOF_WORKING_TYPE));
    printf("   memory used for tensor trains    %ld numbers = %.3f gb\n", mem_used_tensor_trains,
           size_t_to_gygabytes(mem_used_tensor_trains * SIZEOF_WORKING_TYPE));
    printf("   compression rate       = %.3f\n", ((double) mem_used_linear_arrays) / ((double) mem_used_tensor_trains));
    printf("   max restoration error  = %.3e\n", max_restoration_error);
    printf("   time for decomposition = %.3f sec\n", time_for_decomposition);
    printf("   time for restoration   = %.3f sec\n", time_for_restoration);
    printf("\n");
}


double tensor_train_restoration_error(size_t n, void *original, void *restored)
{
    if (arith == CC_ARITH_COMPLEX) {
        double complex *z_original = original;
        double complex *z_restored = restored;

        double max_diff = 0.0;
        for (size_t i = 0; i < n; i++) {
            double diff = cabs(z_original[i] - z_restored[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }

        return max_diff;
    }
    else {
        double *d_original = original;
        double *d_restored = restored;

        double max_diff = 0.0;
        for (size_t i = 0; i < n; i++) {
            double diff = fabs(d_original[i] - d_restored[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }

        return max_diff;
    }
}


double size_t_to_gygabytes(size_t num_bytes)
{
    return ((double) num_bytes) / (1024.0 * 1024.0 * 1024.0);
}


#endif // TENSOR_TRAIN

