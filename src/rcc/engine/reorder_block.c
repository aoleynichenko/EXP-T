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

/**
 * performs reordering of the symmetry block.
 * implemented as a "template" function (C++ -style).
 * This file is to be #included from reorder.c.
 *
 * For more details on the "template technique" employed, see
 * arnold.uthar.net/index.php?n=Work.TemplatesC
 *
 * 2018-2021 Alexander Oleynichenko
 */

#ifdef TYPENAME

#include "templates.h"


void TEMPLATE(reorder_block, TYPENAME)(block_t *src, block_t *tgt, int *perm)
{
    int dims1[CC_DIAGRAM_MAX_RANK];
    int dims2[CC_DIAGRAM_MAX_RANK];
    int idx[CC_DIAGRAM_MAX_RANK];
    int idx2[CC_DIAGRAM_MAX_RANK];
    size_t i, j, idim;
    size_t coef1[CC_DIAGRAM_MAX_RANK];
    size_t coef2[CC_DIAGRAM_MAX_RANK];
    bool_t cuda_enabled = cc_opts->cuda_enabled;

    block_load(src);
    block_load(tgt);

    int rank = src->rank;
    TYPENAME *v1 = (TYPENAME *) src->buf;
    TYPENAME *v2 = (TYPENAME *) tgt->buf;

    block_get_dims(src, dims1);
    block_get_dims(tgt, dims2);

    // prepare coefficients for recalculation: compound index <-> linear index
    for (i = 0; i < rank; i++) {
        coef1[i] = 1;
        coef2[i] = 1;
    }
    for (i = 0; i < rank - 1; i++) {
        coef1[i] = 1;
        coef2[i] = 1;
        for (j = i + 1; j < rank; j++) {
            coef1[i] *= dims1[j];
            coef2[i] *= dims2[j];
        }
    }

    // for a while: CUDA only for complex arith
#if defined(CUDA_FOUND)
    if (cuda_enabled && rank == 4) {
        int err = reorder_cuda(carith, rank, perm, dims1, dims2, (double complex *)v1, (double complex *)v2);
        if (err > 0) {
            errquit("Error in reorder (CUDA error)");
        }
        goto end_reorder;
    }
#endif /* CUDA_FOUND */

    // special case: rank 4 tensors
    // unrolled inner loops
    if (cc_opts->nthreads > 1) {
        omp_set_num_threads(cc_opts->nthreads);
    }
    else {
        omp_set_num_threads(1);
    }
    if (rank == 4) {
        int dim0 = dims1[0];
        int dim1 = dims1[1];
        int dim2 = dims1[2];
        int dim3 = dims1[3];
        int to[4];
        reverse_perm(4, perm, to);

#pragma omp parallel for schedule(static) shared(coef2, v1)
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
                        v2[index2_3] = v1[index];
                        index++;
                    }
                }
            }
        }

        goto end_reorder;
    }

    // special case: rank 6 tensors
    // unrolled inner loops
    if (rank == 6) {
        int dim0 = dims1[0];
        int dim1 = dims1[1];
        int dim2 = dims1[2];
        int dim3 = dims1[3];
        int dim4 = dims1[4];
        int dim5 = dims1[5];
        int to[6];  // на каком месте окажется данное измерение?
        reverse_perm(6, perm, to);

#pragma omp parallel for schedule(static) shared(coef2, v1)
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
                                v2[index2_5] = v1[index];
                                index++;
                            }
                        }
                    }
                }
            }
        }

        goto end_reorder;
    }

    // general case
    for (i = 0; i < src->size; i++) {
        // linear index 'i' -> compound index 'idx'
        size_t index = i;
        for (idim = 0; idim < rank; idim++) {
            idx[idim] = index / coef1[idim];
            index = index % coef1[idim];
        }
        // reorder compound index
        for (idim = 0; idim < rank; idim++) {
            idx2[idim] = idx[perm[idim]];
        }
        // compound index 'idx2' -> linear index 'index2'
        size_t index2;
        index2 = 0;
        for (idim = 0; idim < rank; idim++) {
            index2 += coef2[idim] * idx2[idim];
        }
        // copy matrix element
        v2[index2] = v1[i];
    }

    end_reorder:
    block_unload(src);
    block_store(tgt);
}


#endif /* TYPENAME defined */
