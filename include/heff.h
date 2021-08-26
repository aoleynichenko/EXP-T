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
 * heff.h
 * ======
 *
 * Tools for construction and diagonalization of the effective Hamiltonian
 * matrix and the subsequent analysis of its eigenvalues and eigenvectors.
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#ifndef CC_HEFF_H_INCLUDED
#define CC_HEFF_H_INCLUDED

#include <complex.h>
#include <stdarg.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

#include "engine.h"
#include "options.h"
#include "../src/heff/slater_rules.h"

typedef struct {
    double complex eigval;
    int repno;
} eigval_t;

void heff_analysis(int sect_h, int sect_p, ...);

int eigval_cmp(const void *aa, const void *bb);

void density_matrix(int sect_h, int sect_p, int rep1, int state1, int rep2, int state2);

void dipole_length_tdms(int sect_h, int sect_p);

void model_space_property(cc_ms_prop_query_t *prop_query);

void model_space_properties_and_natural_orbitals();

void print_model_vector(FILE *f_out, int sect_h, int sect_p,
                        int rep_no, char *rep_name, int state_no, double complex energy,
                        size_t len, double complex *coeffs, slater_det_t *det_list, double coef_thresh);

void construct_diagonalize_heff_silent(int sect_h, int sect_p, ...);

void construct_heff(int sect_h, int sect_p, slater_det_t **det_list, size_t *block_sizes, double complex **heff_blocks,
                    int n_diagrams, char **diagram_names);
void diagonalize_heff(int sect_h, int sect_p, size_t *block_dims, double complex **heff,
                      double complex **eigvalues, double complex **coef_left, double complex **coef_right);
double get_lowest_eigenvalue(size_t *block_dims, double complex **eigvalues);

void write_formatted_heff_0h0p(double total_energy);

#endif /* CC_HEFF_H_INCLUDED */
