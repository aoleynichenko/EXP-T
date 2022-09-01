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

/*
 * Division of the diagram by the energy denominators (DIVide by EPSilon).
 *  + basic denominators;
 *  + "intermediate Hamiltonians";
 *  + Andrei Zaitsevskii's dynamic shifts.
 *
 * 2018-2021 Alexander Oleynichenko
 */

#include <math.h>

#include "engine.h"
#include "datamodel.h"
#include "error.h"
#include "options.h"
#include "spinors.h"
#include "timer.h"

static int *curr_valence;

static const double ZERO_THRESH = 1e-14;

void diveps_block_rank_2(block_t *block);

void diveps_block_rank_4(block_t *block);

void diveps_block_rank_6(block_t *block);

void diveps_block_general(block_t *block);

int get_attenuation_parameter(int sect_h, int sect_p);

cc_shifttype_t get_shift_type(int sect_h, int sect_p);

double get_shift_value(int sect_h, int sect_p, int rank, int *indices);

void divide_with_shift(double complex *val, double denom, int shift_type, double shift, int npower);

int count_active_indices(int rank, int *indices);

int count_active_intermediate_indices(int rank, int *indices);

int count_active_interm_quasipart_annihilation_indices(int sect_h, int sect_p, int rank, int *indices);

/**
 * Division of the diagram by the energy denominators
 */
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


/**
 * Division of the diagram by the energy denominators:
 * V[in|out] = V[in|out] / D_K
 * D_K = sum eps[in] - sum eps[out]
 */
diagram_t *diagram_diveps(diagram_t *dg)
{
    int rank = dg->rank;
    curr_valence = dg->valence;

    for (size_t isb = 0; isb < dg->n_blocks; isb++) {
        block_t *block = dg->blocks[isb];
        block_load(block);
        if (block->is_unique == 0) {
            continue;
        }

        if (rank == 2) {
            diveps_block_rank_2(block);
        }
        else if (rank == 4) {
            diveps_block_rank_4(block);
        }
        else if (rank == 6) {
            diveps_block_rank_6(block);
        }
        else {
            diveps_block_general(block);
        }

        block_store(block);
    }

    return dg;
}


void diveps_block_rank_2(block_t *block)
{
    int nspinors = get_num_spinors();
    int sect_h = cc_opts->curr_sector_h;
    int sect_p = cc_opts->curr_sector_p;
    int shift_type = get_shift_type(sect_h, sect_p);
    int npower = get_attenuation_parameter(sect_h, sect_p);

    int dims_0 = block->indices[0][0];
    int dims_1 = block->indices[1][0];
    int coef_0 = dims_1;
    int *block_indices_0 = block->indices[0] + 1; // +1 to skip first element (length)
    int *block_indices_1 = block->indices[1] + 1;

    double *eps = (double *) cc_malloc(nspinors * sizeof(double));
    get_spinor_energies(nspinors, eps);

    // loop over matrix elements
    for (size_t i = 0; i < block->size; i++) {
        // if matrix element is zero we can skip it
        double complex t = (arith == CC_ARITH_COMPLEX) ? block->buf[i] : ((double *) block->buf)[i] + 0.0 * I;
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
        double denom = eps[idx_0] - eps[idx_1];

        extern int diveps_lambda;
        if (diveps_lambda) {
            denom *= -1;
        }

        if (shift_type == CC_SHIFT_NONE) {
            t /= denom;
        }
        else {
            int indices[2];
            indices[0] = idx_0;
            indices[1] = idx_1;
            double shift = get_shift_value(sect_h, sect_p, 2, indices);
            divide_with_shift(&t, denom, shift_type, shift, npower);
        }

        if (arith == CC_ARITH_COMPLEX) {
            block->buf[i] = t;
        }
        else {
            ((double *) block->buf)[i] = creal(t);
        }
    }

    cc_free(eps);
}


void diveps_block_rank_4(block_t *block)
{
    int nspinors = get_num_spinors();
    int sect_h = cc_opts->curr_sector_h;
    int sect_p = cc_opts->curr_sector_p;
    int shift_type = get_shift_type(sect_h, sect_p);
    int npower = get_attenuation_parameter(sect_h, sect_p);

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

    double *eps = (double *) cc_malloc(nspinors * sizeof(double));
    get_spinor_energies(nspinors, eps);

    // loop over matrix elements
    for (size_t i = 0; i < block->size; i++) {
        // if matrix element is zero we can skip it
        double complex t = (arith == CC_ARITH_COMPLEX) ? block->buf[i] : ((double *) block->buf)[i] + 0.0 * I;
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
        double denom = eps[idx_0] + eps[idx_1] - eps[idx_2] - eps[idx_3];


        extern int diveps_lambda;
        if (diveps_lambda) {
            denom *= -1;
        }

        if (shift_type == CC_SHIFT_NONE) {
            t /= denom;
        }
        else {
            int indices[4];
            indices[0] = idx_0;
            indices[1] = idx_1;
            indices[2] = idx_2;
            indices[3] = idx_3;
            double shift = get_shift_value(sect_h, sect_p, 4, indices);
            /*if (fabs(shift) > 1e-3) {
                printf("denom %d-%d-%d-%d = %f\n", idx_0, idx_1, idx_2, idx_3, denom);
                printf("shift = %f\n", shift);
            }*/
            divide_with_shift(&t, denom, shift_type, shift, npower);
        }

        if (arith == CC_ARITH_COMPLEX) {
            block->buf[i] = t;
        }
        else {
            ((double *) block->buf)[i] = creal(t);
        }
    }

    cc_free(eps);
}


void diveps_block_rank_6(block_t *block)
{
    int nspinors = get_num_spinors();
    int sect_h = cc_opts->curr_sector_h;
    int sect_p = cc_opts->curr_sector_p;
    int shift_type = get_shift_type(sect_h, sect_p);
    int npower = get_attenuation_parameter(sect_h, sect_p);

    int dims_1 = block->indices[0][0];
    int dims_2 = block->indices[1][0];
    int dims_3 = block->indices[2][0];
    int dims_4 = block->indices[3][0];
    int dims_5 = block->indices[4][0];
    int dims_6 = block->indices[5][0];

    int *block_indices_1 = block->indices[0] + 1; // +1 to skip first element (length)
    int *block_indices_2 = block->indices[1] + 1;
    int *block_indices_3 = block->indices[2] + 1;
    int *block_indices_4 = block->indices[3] + 1;
    int *block_indices_5 = block->indices[4] + 1;
    int *block_indices_6 = block->indices[5] + 1;

    double *eps = (double *) cc_malloc(nspinors * sizeof(double));
    get_spinor_energies(nspinors, eps);

    size_t index = 0;
    for (int i1 = 0; i1 < dims_1; i1++) {
        int idx_1 = block_indices_1[i1]; // relative index to spinor indices
        double denom_1 = eps[idx_1];
        for (int i2 = 0; i2 < dims_2; i2++) {
            int idx_2 = block_indices_2[i2];
            double denom_2 = denom_1 + eps[idx_2];
            for (int i3 = 0; i3 < dims_3; i3++) {
                int idx_3 = block_indices_3[i3];
                double denom_3 = denom_2 + eps[idx_3];
                for (int i4 = 0; i4 < dims_4; i4++) {
                    int idx_4 = block_indices_4[i4];
                    double denom_4 = denom_3 - eps[idx_4];
                    for (int i5 = 0; i5 < dims_5; i5++) {
                        int idx_5 = block_indices_5[i5];
                        double denom_5 = denom_4 - eps[idx_5];
                        for (int i6 = 0; i6 < dims_6; i6++) {
                            int idx_6 = block_indices_6[i6];
                            double denom = denom_5 - eps[idx_6];

                            double complex t = (arith == CC_ARITH_COMPLEX) ?
                                               block->buf[index] : ((double *) block->buf)[index] + 0.0 * I;
                            if (cabs(t) < ZERO_THRESH) {
                                index++;
                                continue;
                            }

                            // calculate denominator and divide matrix element by it
                            //denom = eps[idx_1] + eps[idx_2] + eps[idx_3] - eps[idx_4] - eps[idx_5] - eps[idx_6];
                            if (shift_type == CC_SHIFT_NONE) {
                                t /= denom;
                            }
                            else {
                                int indices[6];
                                indices[0] = idx_1;
                                indices[1] = idx_2;
                                indices[2] = idx_3;
                                indices[3] = idx_4;
                                indices[4] = idx_5;
                                indices[5] = idx_6;
                                double shift = get_shift_value(sect_h, sect_p, 6, indices);
                                divide_with_shift(&t, denom, shift_type, shift, npower);
                            }

                            if (arith == CC_ARITH_COMPLEX) {
                                block->buf[index] = t;
                            }
                            else {
                                ((double *) block->buf)[index] = creal(t);
                            }
                            index++;
                        }
                    }
                }
            }
        }
    }

    cc_free(eps);
}


// general case: n-dim tensor
void diveps_block_general(block_t *block)
{
    int dims[CC_DIAGRAM_MAX_RANK];
    int idx[CC_DIAGRAM_MAX_RANK];

    int rank = block->rank;
    int nspinors = get_num_spinors();
    int sect_h = cc_opts->curr_sector_h;
    int sect_p = cc_opts->curr_sector_p;
    int shift_type = get_shift_type(sect_h, sect_p);
    int npower = get_attenuation_parameter(sect_h, sect_p);

    double *eps = (double *) cc_malloc(get_num_spinors() * sizeof(double));
    get_spinor_energies(nspinors, eps);
    block_get_dims(block, dims);

    // loop over matrix elements
    for (size_t i = 0; i < block->size; i++) {
        as_compound_index(rank, dims, i /* linear index */, idx);
        // relative indices -> spinor indices
        for (int j = 0; j < rank; j++) {
            idx[j] = block->indices[j][idx[j] + 1];
        }
        double denom = 0.0;
        for (int j = 0; j < rank / 2; j++) {
            denom += eps[idx[j]];
        }
        for (int j = rank / 2; j < rank; j++) {
            denom -= eps[idx[j]];
        }

        // get amplitude
        double complex t = (arith == CC_ARITH_COMPLEX) ?
                           block->buf[i] : ((double *) block->buf)[i] + 0.0 * I;
        if (cabs(t) < ZERO_THRESH) {
            continue;
        }

        // divide by denominator
        if (shift_type == CC_SHIFT_NONE) {
            t /= denom;
        }
        else {
            double shift = get_shift_value(sect_h, sect_p, rank, idx);
            divide_with_shift(&t, denom, shift_type, shift, npower);
        }

        // update amplitude
        if (arith == CC_ARITH_COMPLEX) {
            block->buf[i] = t;
        }
        else {
            ((double *) block->buf)[i] = creal(t);
        }
    }

    cc_free(eps);
}


cc_shifttype_t get_shift_type(int sect_h, int sect_p)
{
    if (cc_opts->do_intham_imms) {
        return cc_opts->intham_imms_opts.shift_type;
    }
    else {
        return cc_opts->shifts[sect_h][sect_p].type;
    }
}


int get_attenuation_parameter(int sect_h, int sect_p)
{
    if (cc_opts->do_intham_imms) {
        return cc_opts->intham_imms_opts.npower;
    }
    else {
        return cc_opts->shifts[sect_h][sect_p].power;
    }
}


/*
 * Encapsulates some specific types of shift applied to the given amplitude.
 * Shifts can depend on:
 * - sector
 * - cluster operator rank
 * - number of annihilation lines of active quasiparticles
 * - number of annihilation lines of active intermediate quasiparticles
 */
double get_shift_value(int sect_h, int sect_p, int rank, int *indices)
{
    /*
     * "IH-1": intermediate-Hamiltonian-like shifting technique
     */
    if (cc_opts->do_intham_imms) {
        if (cc_opts->intham_imms_opts.sectors[sect_h][sect_p] == 1) {
            double shift = intham_imms_get_shift_value(sect_h, sect_p, rank, indices, curr_valence);
            return shift;
        }
        else {
            return 0.0;
        }
    }

    /*
     * default: constant shifts, which don't depend on the amplitude or its indices
     */
    double shift = cc_opts->shifts[sect_h][sect_p].shifts[rank / 2 - 1];
    return shift;
}


/**
 * Divides cluster amplitude 'val' by energy shifted energy denominator.
 * Basic papers:
 * (1) A. Zaitsevskii, N. S. Mosyagin, A. V. Stolyarov, E. Eliav.
 *     Phys. Rev. A, V. 96, P. 022516 (2017)
 * (2) A. Zaitsevskii, E. Eliav. Int. J. Quantum Chem. V. 118, P. e25772 (2018)
 */
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

