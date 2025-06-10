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

/**
 * implementation of a tensor transposition operation
 * (reordering of dimensions of a tensor).
 * implemented as a set of "template" functions (C++ style).
 * This file is to be #included from reorder.c.
 *
 * For more details on the "template technique" employed, see
 * arnold.uthar.net/index.php?n=Work.TemplatesC
 */


#ifdef TYPENAME

#include "tensor.h"
#include "comdef.h"

#include "c_templates.h"


#ifndef COMPILER_CLANG
#include "omp.h"
#else

void omp_set_num_threads(int num_threads);

void openblas_set_num_threads(int num_threads);

#endif // ifndef COMPILER_CLANG

void reverse_perm(int n, const int *direct_perm, int *inv_perm);

static void TEMPLATE(tensor_transpose_rank4, TYPENAME)(const TYPENAME *tensor, const int *shape, const int *new_shape, const int *perm, TYPENAME *transposed_tensor);

static void TEMPLATE(tensor_transpose_rank6, TYPENAME)(const TYPENAME *tensor, const int *shape, const int *new_shape, const int *perm, TYPENAME *transposed_tensor);

static void TEMPLATE(tensor_transpose_general_rank, TYPENAME)(int rank, const TYPENAME *tensor, const int *shape, const int *new_shape, const int *perm, TYPENAME *transposed_tensor);


void TEMPLATE(tensor_transpose, TYPENAME)(int rank, const TYPENAME *tensor, const int *shape, const int *new_shape, const int *perm,
    TYPENAME *transposed_tensor, int nthreads)
{
    if (nthreads > 1) {
        omp_set_num_threads(nthreads);
    }
    else {
        omp_set_num_threads(1);
    }

    if (rank == 4) {
        TEMPLATE(tensor_transpose_rank4, TYPENAME)(tensor, shape, new_shape, perm, transposed_tensor);
    }
    else if (rank == 6) {
        TEMPLATE(tensor_transpose_rank6, TYPENAME)(tensor, shape, new_shape, perm, transposed_tensor);
    }
    else {
        TEMPLATE(tensor_transpose_general_rank, TYPENAME)(rank, tensor, shape, new_shape, perm, transposed_tensor);
    }
}


static void TEMPLATE(tensor_transpose_rank4, TYPENAME)(const TYPENAME *tensor, const int *shape, const int *new_shape, const int *perm, TYPENAME *transposed_tensor)
{
    // prepare coefficients for recalculation: compound index <-> linear index
    size_t coef2[CC_DIAGRAM_MAX_RANK];
    coef2[0] = new_shape[1] * new_shape[2] * new_shape[3];
    coef2[1] = new_shape[2] * new_shape[3];
    coef2[2] = new_shape[3];
    coef2[3] = 1;

    int dim0 = shape[0];
    int dim1 = shape[1];
    int dim2 = shape[2];
    int dim3 = shape[3];
    int to[4];
    reverse_perm(4, perm, to);

#pragma omp parallel for schedule(static) shared(coef2, tensor)
    for (int i0 = 0; i0 < dim0; i0++) {
        size_t index = i0 * dim1 * dim2 * dim3;
        size_t index2 = i0 * coef2[to[0]];

        for (int i1 = 0; i1 < dim1; i1++) {
            size_t index2_1 = index2 + i1 * coef2[to[1]];

            for (int i2 = 0; i2 < dim2; i2++) {
                size_t index2_2 = index2_1 + i2 * coef2[to[2]];

                for (int i3 = 0; i3 < dim3; i3++) {
                    size_t index2_3 = index2_2 + i3 * coef2[to[3]];

                    // copy matrix element
                    transposed_tensor[index2_3] = tensor[index];
                    index++;
                }
            }
        }
    }
}


static void TEMPLATE(tensor_transpose_rank6, TYPENAME)(const TYPENAME *tensor, const int *shape, const int *new_shape, const int *perm, TYPENAME *transposed_tensor)
{
    // prepare coefficients for recalculation: compound index <-> linear index
    size_t coef2[CC_DIAGRAM_MAX_RANK];
    coef2[0] = new_shape[1] * new_shape[2] * new_shape[3] * new_shape[4] * new_shape[5];
    coef2[1] = new_shape[2] * new_shape[3] * new_shape[4] * new_shape[5];
    coef2[2] = new_shape[3] * new_shape[4] * new_shape[5];
    coef2[3] = new_shape[4] * new_shape[5];
    coef2[4] = new_shape[5];
    coef2[5] = 1;

    int dim0 = shape[0];
    int dim1 = shape[1];
    int dim2 = shape[2];
    int dim3 = shape[3];
    int dim4 = shape[4];
    int dim5 = shape[5];
    int to[6];  // на каком месте окажется данное измерение?
    reverse_perm(6, perm, to);

#pragma omp parallel for schedule(static) shared(coef2, tensor)
    for (int i0 = 0; i0 < dim0; i0++) {
        size_t index = i0 * dim1 * dim2 * dim3 * dim4 * dim5;
        size_t index2 = i0 * coef2[to[0]];

        for (int i1 = 0; i1 < dim1; i1++) {
            size_t index2_1 = index2 + i1 * coef2[to[1]];

            for (int i2 = 0; i2 < dim2; i2++) {
                size_t index2_2 = index2_1 + i2 * coef2[to[2]];

                for (int i3 = 0; i3 < dim3; i3++) {
                    size_t index2_3 = index2_2 + i3 * coef2[to[3]];

                    for (int i4 = 0; i4 < dim4; i4++) {
                        size_t index2_4 = index2_3 + i4 * coef2[to[4]];

                        for (int i5 = 0; i5 < dim5; i5++) {
                            size_t index2_5 = index2_4 + i5 * coef2[to[5]];

                            // copy matrix element
                            transposed_tensor[index2_5] = tensor[index];
                            index++;
                        }
                    }
                }
            }
        }
    }
}


static void TEMPLATE(tensor_transpose_general_rank, TYPENAME)(int rank, const TYPENAME *tensor, const int *shape, const int *new_shape, const int *perm,
                                          TYPENAME *transposed_tensor)
{
    int idx[CC_DIAGRAM_MAX_RANK];
    int idx2[CC_DIAGRAM_MAX_RANK];
    size_t coef1[CC_DIAGRAM_MAX_RANK];
    size_t coef2[CC_DIAGRAM_MAX_RANK];

    // prepare coefficients for recalculation: compound index <-> linear index
    for (int i = 0; i < rank; i++) {
        coef1[i] = 1;
        coef2[i] = 1;
    }
    for (int i = 0; i < rank - 1; i++) {
        coef1[i] = 1;
        coef2[i] = 1;
        for (int j = i + 1; j < rank; j++) {
            coef1[i] *= shape[j];
            coef2[i] *= new_shape[j];
        }
    }

    // general case
    size_t num_elements = tensor_num_elements(rank, shape);
    for (size_t i = 0; i < num_elements; i++) {
        // linear index 'i' -> compound index 'idx'
        size_t index = i;
        for (int idim = 0; idim < rank; idim++) {
            idx[idim] = index / coef1[idim];
            index = index % coef1[idim];
        }

        // reorder compound index
        for (int idim = 0; idim < rank; idim++) {
            idx2[idim] = idx[perm[idim]];
        }

        // compound index 'idx2' -> linear index 'index2'
        size_t index2;
        index2 = 0;
        for (int idim = 0; idim < rank; idim++) {
            index2 += coef2[idim] * idx2[idim];
        }

        // copy matrix element
        transposed_tensor[index2] = tensor[i];
    }
}

#endif /* TYPENAME defined */
