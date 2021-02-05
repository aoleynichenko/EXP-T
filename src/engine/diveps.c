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
 * diveps.c
 *
 * Division of the diagram by the energy denominators (DIVide by EPSilon).
 *  + basic denominators;
 *  - intermediate Hamiltonians;
 *  + Andrei Zaitsevskii's dynamic shifts.
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "engine.h"
#include "datamodel.h"
#include "error.h"
#include "options.h"
#include "spinors.h"
#include "timer.h"

const double ZERO_THRESH = 1e-14;

void divide_with_shift(double complex *val, double denom, int shift_type, double shift, int npower);


/*******************************************************************************
 * diveps
 *
 * Division of the diagram by the energy denominators
 ******************************************************************************/
void diveps(char *name)
{
    diagram_t *dg;

    timer_new_entry("diveps", "Energy denominators (diveps)");
    timer_start("diveps");

    dg = diagram_stack_find(name);
    if (dg == NULL) {
        errquit("diveps(): diagram '%s' not found", name);
    }

    diagram_diveps(dg);

    timer_stop("diveps");
}


/*******************************************************************************
 * diagram_diveps
 *
 * Division of the diagram by the energy denominators:
 * V[in|out] = V[in|out] / D_K
 * D_K = sum eps[in] - sum eps[out]
 ******************************************************************************/
diagram_t *diagram_diveps(diagram_t *dg)
{
    size_t isb;
    size_t i, j;
    block_t *block;
    double complex *buf;
    size_t size;
    double denom;
    int idx[CC_DIAGRAM_MAX_RANK];
    int dims[CC_DIAGRAM_MAX_RANK];
    double *eps;
    int do_orbshift = 0;

    int h = cc_opts->curr_sector_h;
    int p = cc_opts->curr_sector_p;
    int rank = dg->rank;
    int nparts = rank / 2;
    cc_shifttype_t shift_type = cc_opts->shifts[h][p].type;
    int npower = cc_opts->shifts[h][p].power;
    double shift = cc_opts->shifts[h][p].shifts[nparts - 1];

    // orbital shifts: its own options
    if (cc_opts->do_orbshift && (h != 0 || p != 0)) {
        do_orbshift = 1;
        shift = cc_opts->orbshift;
        npower = cc_opts->orbshift_power;
    }
    else if (cc_opts->do_orbshift_0h0p == 1 && (h == 0 && p == 0)) {
        do_orbshift = 1;
        shift = cc_opts->orbshift_0h0p[0];
        npower = cc_opts->orbshift_0h0p_power;
    }
    else if (cc_opts->do_orbshift_0h0p == 2 && (h == 0 && p == 0)) {
        do_orbshift = 2;
        shift = cc_opts->orbshift_0h0p[0];
        npower = cc_opts->orbshift_0h0p_power;
    }

    // extract one-electron energies
    eps = (double *) cc_malloc(nspinors * sizeof(double));
    for (i = 0; i < nspinors; i++) {
        eps[i] = spinor_info[i].eps;
    }

    for (isb = 0; isb < dg->n_blocks; isb++) {
        block = dg->blocks[isb];
        size = block->size;
        symblock_get_dims(block, dims);
        symblock_load(block);
        buf = block->buf;
        if (block->is_unique == 0) {
            continue;
        }

        // special case: 2-dim tensors
        // hand-written loop unrolling
        if (rank == 2) {
            int dims_0 = block->indices[0][0];
            int dims_1 = block->indices[1][0];
            int coef_0 = dims_1;
            int *block_indices_0 = block->indices[0] + 1; // +1 to skip first element (length)
            int *block_indices_1 = block->indices[1] + 1;

            // loop over matrix elements
            for (i = 0; i < size; i++) {
                // if matrix element is zero we can skip it
                double complex t = carith ? buf[i] : ((double *) buf)[i] + 0.0 * I;
                if (cabs(t) < ZERO_THRESH) {
                    continue;
                }
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
                if (shift_type == CC_SHIFT_NONE) {
                    t /= denom;
                }
                else if (do_orbshift == 1) {
                    // count number of valence indices
                    int nvalence = is_active(idx_0) + is_active(idx_1);
                    double actual_shift = nvalence * (0.5 * shift); // n*s/2
                    divide_with_shift(&t, denom, shift_type, actual_shift, npower);
                }
                else if (do_orbshift == 2) {
                    // count number of valence indices
                    int nvalence = is_active(idx_0) + is_active(idx_1);
                    double actual_shift = nvalence > 0 ? cc_opts->orbshift_0h0p[nvalence - 1] : 0.0;
                    divide_with_shift(&t, denom, shift_type, actual_shift, npower);
                }
                else {
                    divide_with_shift(&t, denom, shift_type, shift, npower);
                }
                if (carith) {
                    buf[i] = t;
                }
                else {
                    ((double *) buf)[i] = creal(t);
                }
            }
        }
            // special case: 4-dim tensors
            // hand-written loop unrolling
        else if (rank == 4) {
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
                // if matrix element is zero we can skip it
                double complex t = carith ? buf[i] : ((double *) buf)[i] + 0.0 * I;
                if (cabs(t) < ZERO_THRESH) {
                    continue;
                }
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
                if (shift_type == CC_SHIFT_NONE) {
                    t /= denom;
                }
                else if (do_orbshift == 1) {
                    // count number of valence indices
                    int nvalence = is_active(idx_0) + is_active(idx_1) +
                                   is_active(idx_2) + is_active(idx_3);
                    double actual_shift = nvalence * (0.5 * shift); // n*s/2
                    divide_with_shift(&t, denom, shift_type, actual_shift, npower);
                }
                else if (do_orbshift == 2) {
                    // count number of valence indices
                    int nvalence = is_active(idx_0) + is_active(idx_1) +
                                   is_active(idx_2) + is_active(idx_3);
                    double actual_shift = nvalence > 0 ? cc_opts->orbshift_0h0p[2 + nvalence - 1] : 0.0;
                    divide_with_shift(&t, denom, shift_type, actual_shift, npower);
                }
                else {
                    divide_with_shift(&t, denom, shift_type, shift, npower);
                }
                if (carith) {
                    buf[i] = t;
                }
                else {
                    ((double *) buf)[i] = creal(t);
                }
            }
        }
            // special case: 6-dim tensors
            // hand-written loop unrolling
        else if (rank == 6) {
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

                                    double complex t = carith ?
                                                       buf[index] : ((double *) buf)[index] + 0.0 * I;
                                    if (cabs(t) < ZERO_THRESH) {
                                        index++;
                                        continue;
                                    }

                                    // calculate denominator and divide matrix element by it
                                    //denom = eps[idx_0] + eps[idx_1] + eps[idx_2] - eps[idx_3] - eps[idx_4] - eps[idx_5];
                                    if (shift_type == CC_SHIFT_NONE) {
                                        t /= denom;
                                    }
                                    else if (do_orbshift) {
                                        // count number of valence indices
                                        int nvalence = is_active(idx_0) + is_active(idx_1) +
                                                       is_active(idx_2) + is_active(idx_3) +
                                                       is_active(idx_4) + is_active(idx_5);
                                        double actual_shift = nvalence * (0.5 * shift); // n*s/2
                                        divide_with_shift(&t, denom, shift_type, actual_shift, npower);
                                    }
                                    else {
                                        divide_with_shift(&t, denom, shift_type, shift, npower);
                                    }

                                    if (carith) {
                                        buf[index] = t;
                                    }
                                    else {
                                        ((double *) buf)[index] = creal(t);
                                    }
                                    index++;
                                }
                            }
                        }
                    }
                }
            }
        }
            // general case: n-dim tensor
            // @todo generalize for the case of real arith!
        else {
            assert(carith);
            // loop over matrix elements
            for (i = 0; i < size; i++) {
                as_compound_index(rank, dims, i /* linear index */, idx);
                // relative indices -> spinor indices
                for (j = 0; j < rank; j++) {
                    idx[j] = block->indices[j][idx[j] + 1];
                }
                denom = 0.0;
                for (j = 0; j < rank / 2; j++) {
                    denom += eps[idx[j]];
                }
                for (j = rank / 2; j < rank; j++) {
                    denom -= eps[idx[j]];
                }
                // divide by denominator
                if (shift_type == CC_SHIFT_NONE) {
                    buf[i] /= denom;
                }
                else if (do_orbshift) {
                    // count number of valence indices
                    int nvalence = 0;
                    for (j = 0; j < rank; j++) {
                        if (is_active(idx[j])) {
                            nvalence++;
                        }
                    }
                    shift = nvalence * (0.5 * shift); // n*s/2
                    divide_with_shift(&buf[i], denom, shift_type, shift, npower);
                }
                else {
                    divide_with_shift(&buf[i], denom, shift_type, shift, npower);
                }
            }
        }
        symblock_store(block);
    }

    cc_free(eps);

    return dg;
}


/*******************************************************************************
 * divide_with_shift
 *
 * Divides cluster amplitude 'val' by energy shifted energy denominator.
 * Basic papers:
 * (1) A. Zaitsevskii, N. S. Mosyagin, A. V. Stolyarov, E. Eliav.
 *     Phys. Rev. A, V. 96, P. 022516 (2017)
 * (2) A. Zaitsevskii, E. Eliav. Int. J. Quantum Chem. V. 118, P. e25772 (2018)
 ******************************************************************************/
void divide_with_shift(double complex *val, double denom, int shift_type, double shift, int npower)
{
    double complex cdenom;
    double top, bot;

    switch (shift_type) {
        case CC_SHIFT_NONE:
            *val = *val / denom;
            break;
        case CC_SHIFT_REAL:
            // real shift
            // D' = D + S*(S/(D+S))^m
            denom = denom + shift * pow(shift / (denom + shift), npower);
            *val = *val / denom;
            break;
        case CC_SHIFT_REALIMAG:
            // real simulation of imaginary shift
            // D' = D + S^2/D * (S^2/(D^2+S^2))^m
            denom = denom + (shift) * (shift) / denom *
                            pow(shift * shift / (denom * denom + shift * shift), npower);
            *val = *val / denom;
            break;
        case CC_SHIFT_IMAG:
            // imaginary shift
            // D' = D + iS*(|S|/|D+iS|)^m
            cdenom = denom + I * shift * pow(fabs(shift) / cabs(denom + I * shift), npower);
            *val = *val / cdenom;
            break;
        case CC_SHIFT_TAYLOR:
            // "taylor" shift
            // D' = (D+S)*(1-S/(D+S))/(1-(S/(D+S))^(m+1))
            top = 1.0 - shift / (denom + shift);
            bot = 1.0 - pow(shift / (denom + shift), npower + 1);
            denom = (denom + shift) * top / bot;
            *val = *val / denom;
            break;
        default:
            errquit("in divide_with_shift(): unknown type of shift");
            break;
    }
}
