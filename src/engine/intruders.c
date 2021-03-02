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
 * intruders.c
 *
 * Analysis of possible intruder states.
 *
 * 2020-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "engine.h"
#include "datamodel.h"
#include "error.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "timer.h"

typedef struct {
    int idx[CC_DIAGRAM_MAX_RANK];
    double denom;
} indexed_denom_t;


void try_denom(indexed_denom_t *denoms, int ndenoms, int dim, int *idx, double d)
{
    for (int i = 0; i < ndenoms; i++) {
        if (d > denoms[i].denom) {
            denoms[i].denom = d;
            for (int j = 0; j < dim; j++) {
                denoms[i].idx[j] = idx[j];
            }
            return;
        }
            // omit repeating max values
        else if (fabs(d - denoms[i].denom) < 1e-6) {
            return;
        }
    }
}


/**
 * Denominator-based prediction of possible intruder states.
 *
 * @param diagram_name symbolic name of the diagram containing
 *      some cluster operator
 * @param nmax print information about only 'nmax' intruders
 *      with largest energy denominators
 */
void predict_intruders(char *diagram_name, int nmax)
{
    size_t i, j;
    double complex *buf;
    double denom;
    int idx[CC_DIAGRAM_MAX_RANK];
    int dims[CC_DIAGRAM_MAX_RANK];
    indexed_denom_t *max_denoms;
    int occ, active, repno;
    double spinor_e;
    int rk;

    timer_new_entry("intruders", "Prediction of intruder states");
    timer_start("intruders");

    diagram_t *dg = diagram_stack_find(diagram_name);
    rk = rank(diagram_name);
    if (dg == NULL) {
        errquit("predict_intruders(): diagram '%s' not found", diagram_name);
    }

    int h = cc_opts->curr_sector_h;
    int p = cc_opts->curr_sector_p;

    max_denoms = (indexed_denom_t *) cc_malloc(sizeof(indexed_denom_t) * nmax);
    memset(max_denoms, 0, sizeof(indexed_denom_t) * nmax);

    // extract one-electron energies
    double *eps = (double *) cc_malloc(NSPINORS * sizeof(double));
    for (i = 0; i < NSPINORS; i++) {
        eps[i] = spinor_info[i].eps;
    }

    for (size_t isb = 0; isb < dg->n_blocks; isb++) {
        block_t *block = dg->blocks[isb];
        size_t size = block->size;
        symblock_get_dims(block, dims);
        symblock_load(block);
        buf = block->buf;
        if (block->is_unique == 0) {
            continue;
        }

        // special case: 2-dim tensors
        // hand-written loop unrolling
        if (rk == 2) {
            int dims_0 = block->indices[0][0];
            int dims_1 = block->indices[1][0];
            int coef_0 = dims_1;
            int *block_indices_0 = block->indices[0] + 1; // +1 to skip first element (length)
            int *block_indices_1 = block->indices[1] + 1;

            // loop over matrix elements
            for (i = 0; i < size; i++) {
                // calculate compound index
                int index = i;
                int idx_0 = index / coef_0;
                index = index % coef_0;
                int idx_1 = index;    // since coef_1 == 1
                // relative index to spinor indices
                idx_0 = block_indices_0[idx_0];
                idx_1 = block_indices_1[idx_1];
                // calculate denominator and divide matrix element by it
                denom = eps[idx_0] - eps[idx_1];
                // cancel all-closed terms
                if (is_active(idx_0) && is_active(idx_1)) {
                    continue;
                }
                idx[0] = idx_0;
                idx[1] = idx_1;
                try_denom(max_denoms, nmax, 2, idx, denom);
            }
        }
            // special case: 4-dim tensors
            // hand-written loop unrolling
        else if (rk == 4) {
            int dims_0 = block->indices[0][0];
            int dims_1 = block->indices[1][0];
            int dims_2 = block->indices[2][0];
            int dims_3 = block->indices[3][0];
            int coef_0 = dims_1 * dims_2 * dims_3;
            int coef_1 = dims_2 * dims_3;
            int coef_2 = dims_3;
            //int coef_3 = 1;
            int *block_indices_0 = block->indices[0] + 1; // +1 to skip first element (length)
            int *block_indices_1 = block->indices[1] + 1;
            int *block_indices_2 = block->indices[2] + 1;
            int *block_indices_3 = block->indices[3] + 1;

            // loop over matrix elements
            for (i = 0; i < size; i++) {
                // calculate compound index
                int index = i;
                int idx_0 = index / coef_0;
                index = index % coef_0;
                int idx_1 = index / coef_1;
                index = index % coef_1;
                int idx_2 = index / coef_2;
                index = index % coef_2;
                int idx_3 = index;    // since coef_3 == 1

                // relative index to spinor indices
                idx_0 = block_indices_0[idx_0];
                idx_1 = block_indices_1[idx_1];
                idx_2 = block_indices_2[idx_2];
                idx_3 = block_indices_3[idx_3];

                // calculate denominator and divide matrix element by it
                denom = eps[idx_0] + eps[idx_1] - eps[idx_2] - eps[idx_3];
                // cancel all-closed terms
                if (is_active(idx_0) && is_active(idx_1) &&
                    is_active(idx_2) && is_active(idx_3)) {
                    continue;
                }
                // exclude epv terms
                if (idx_0 == idx_1 || idx_2 == idx_3) {
                    continue;
                }
                idx[0] = idx_0;
                idx[1] = idx_1;
                idx[2] = idx_2;
                idx[3] = idx_3;
                try_denom(max_denoms, nmax, 4, idx, denom);
            }
        }
            // special case: 6-dim tensors
            // hand-written loop unrolling
        else if (rk == 6) {
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
                double denom_0 = eps[idx_0];
                for (int i1 = 0; i1 < dims_1; i1++) {
                    int idx_1 = block_indices_1[i1];
                    double denom_1 = denom_0 + eps[idx_1];
                    for (int i2 = 0; i2 < dims_2; i2++) {
                        int idx_2 = block_indices_2[i2];
                        double denom_2 = denom_1 + eps[idx_2];
                        for (int i3 = 0; i3 < dims_3; i3++) {
                            int idx_3 = block_indices_3[i3];
                            double denom_3 = denom_2 - eps[idx_3];
                            for (int i4 = 0; i4 < dims_4; i4++) {
                                int idx_4 = block_indices_4[i4];
                                double denom_4 = denom_3 - eps[idx_4];
                                for (int i5 = 0; i5 < dims_5; i5++) {
                                    int idx_5 = block_indices_5[i5];
                                    double denom = denom_4 - eps[idx_5];
                                    // exclude all-closed terms
                                    if (is_active(idx_0) && is_active(idx_1) &&
                                        is_active(idx_2) && is_active(idx_3) &&
                                        is_active(idx_4) && is_active(idx_5)) {
                                        continue;
                                    }
                                    // exclude epv terms
                                    if (idx_0 == idx_1 || idx_1 == idx_2 || idx_0 == idx_2 ||
                                        idx_3 == idx_4 || idx_4 == idx_5 || idx_3 == idx_5) {
                                        continue;
                                    }
                                    idx[0] = idx_0;
                                    idx[1] = idx_1;
                                    idx[2] = idx_2;
                                    idx[3] = idx_3;
                                    idx[4] = idx_4;
                                    idx[5] = idx_5;
                                    try_denom(max_denoms, nmax, 6, idx, denom);
                                }
                            }
                        }
                    }
                }
            }
        }
            // general case: n-dim tensor
            // @todo generalize for the case of real arith!
            // @todo write code !
        else {
            assert(carith);
            // loop over matrix elements
            for (i = 0; i < size; i++) {
                as_compound_index(rk, dims, i /* linear index */, idx);
                // relative indices -> spinor indices
                for (j = 0; j < rk; j++) {
                    idx[j] = block->indices[j][idx[j] + 1];
                }
                denom = 0.0;
                for (j = 0; j < rk / 2; j++) {
                    denom += eps[idx[j]];
                }
                for (j = rk / 2; j < rk; j++) {
                    denom -= eps[idx[j]];
                }
                // buf[i] bla bla bla ...
                // todo: implement for the general case!
            }
        }
        symblock_store(block);
    }

    // print max denoms

    int n_positive_denoms = 0;
    for (int iden = 0; iden < nmax; iden++) {
        if (max_denoms[iden].denom > 0) {
            n_positive_denoms++;
        }
    }

    if (n_positive_denoms > 0) {
        printf("\n Warning: intruder states are possible!\n");
    }

    for (int iden = 0; iden < nmax && n_positive_denoms > 0; iden++) {
        if (max_denoms[iden].denom <= 0) {
            continue;
        }

        printf(" [%d] Denominator [ ", iden + 1);
        for (i = 0; i < rk / 2; i++) {
            printf("%d ", max_denoms[iden].idx[i] + 1);
        }
        printf("-> ");
        for (i = rk / 2; i < rk; i++) {
            printf("%d ", max_denoms[iden].idx[i] + 1);
        }
        printf("] = ");
        printf("%.4E\n", max_denoms[iden].denom);

        // print more information about these spinors in order to simplify analysis
        for (i = 0; i < rk; i++) {
            get_spinor_info(max_denoms[iden].idx[i], &repno, &occ, &active, &spinor_e);
            printf("     [%4d] eps=%16.8f rep=%s %s %s\n",
                   max_denoms[iden].idx[i] + 1, spinor_e, get_irrep_name(repno),
                   (active == 1) ? "active" : "inactive",
                   (occ == 1) ? "occ" : "virt");
        }
    }

    if (n_positive_denoms > 0) {
        printf("\n");
    }

    cc_free(max_denoms);
    cc_free(eps);

    timer_stop("intruders");
}

