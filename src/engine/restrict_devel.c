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
 * restrict.c
 * ==========
 *
 * Subroutine for restriction of triples-like operators.
 * Idea: we can suppose that only a subset of triples (for example) amplitudes
 * contributes significantly to the answers (energies etc). The simplest way to
 * introduce such a set is just to set up upper and lower limits on one-electron
 * energies (e_min and e_max) and suppose that only triples amplitudes t_ijkabc
 * corresponding to spinors with energies in range [e_min,e_max] are non-zero:
 *   t_ijkabc != 0 only if e_i, e_j, ... e_c are in range [e_min,e_max]
 * Let this scheme be named "restriction of triples".
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <assert.h>
#include <complex.h>
#include <string.h>

#include "platform.h"
#include "datamodel.h"
#include "error.h"
#include "memory.h"
#include "options.h"
#include "spinors.h"
#include "utils.h"


void restrict_triples(char *dg_name, double erange_1, double erange_2)
{
    diagram_t *dg;
    int dims[CC_DIAGRAM_MAX_RANK];
    int rank;
    double complex *buf;

    assert(carith);

    if (cc_opts->restrict_triples == 0) {
        return;
    }

    dg = diagram_stack_find(dg_name);
    if (dg == NULL) {
        errquit("restrict_triples(): diagram '%s' not found", dg_name);
    }
    rank = dg->rank;

    // extract one-electron energies
    double *eps = (double *) cc_malloc(nspinors * sizeof(double));
    for (int i = 0; i < nspinors; i++) {
        eps[i] = spinor_info[i].eps;
    }

    for (size_t isb = 0; isb < dg->n_blocks; isb++) {
        block_t *block = dg->blocks[isb];

        // calculate block's dimensions
        for (int i = 0; i < rank; i++) {
            dims[i] = block->indices[i][0];
        }

        symblock_load(block);

        buf = block->buf;

        // NOTE: only for 6-dim tensors

        int dims_0 = block->indices[0][0];
        int dims_1 = block->indices[1][0];
        int dims_2 = block->indices[2][0];
        int dims_3 = block->indices[3][0];
        int dims_4 = block->indices[4][0];
        int dims_5 = block->indices[5][0];
        int coef_0 = dims_1 * dims_2 * dims_3 * dims_4 * dims_5;
        int coef_1 = dims_2 * dims_3 * dims_4 * dims_5;
        int coef_2 = dims_3 * dims_4 * dims_5;
        int coef_3 = dims_4 * dims_5;
        int coef_4 = dims_5;
        //int coef_5 = 1;
        int *block_indices_0 = block->indices[0] + 1; // +1 to skip first element (length)
        int *block_indices_1 = block->indices[1] + 1;
        int *block_indices_2 = block->indices[2] + 1;
        int *block_indices_3 = block->indices[3] + 1;
        int *block_indices_4 = block->indices[4] + 1;
        int *block_indices_5 = block->indices[5] + 1;

        size_t index = 0;

        for (int i0 = 0; i0 < dims_0; i0++) {
            int idx_0 = block_indices_0[i0]; // relative index to spinor indices
            double e0 = eps[idx_0];
            for (int i1 = 0; i1 < dims_1; i1++) {
                int idx_1 = block_indices_1[i1];
                double e1 = eps[idx_1];
                for (int i2 = 0; i2 < dims_2; i2++) {
                    int idx_2 = block_indices_2[i2];
                    double e2 = eps[idx_2];
                    for (int i3 = 0; i3 < dims_3; i3++) {
                        int idx_3 = block_indices_3[i3];
                        double e3 = eps[idx_3];
                        for (int i4 = 0; i4 < dims_4; i4++) {
                            int idx_4 = block_indices_4[i4];
                            double e4 = eps[idx_4];
                            for (int i5 = 0; i5 < dims_5; i5++) {
                                int idx_5 = block_indices_5[i5];
                                double e5 = eps[idx_5];

                                if (!(in_range(e0, erange_1, erange_2) &&
                                      in_range(e1, erange_1, erange_2) &&
                                      in_range(e2, erange_1, erange_2) &&
                                      in_range(e3, erange_1, erange_2) &&
                                      in_range(e4, erange_1, erange_2) &&
                                      in_range(e5, erange_1, erange_2))) {
                                    buf[index] = 0.0 + 0.0 * I;
                                }

                                index++;
                            }
                        }
                    }
                }
            }
        }

        symblock_store(block);
    }

    cc_free(eps);
}
