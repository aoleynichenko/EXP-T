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
#include <stddef.h>
#include <stdint.h>

#include "engine.h"
#include "options.h"
#include "../src/heff/slater.h"

typedef struct {
    double complex eigval;
    int repno;
} eigval_t;

void diag_heff(int sect_h, int sect_p, ...);

int eigval_cmp(const void *aa, const void *bb);

void eigenvalues_table(size_t n_eigenvalues, eigval_t *eigenvalues, double degen_thresh, double max_energy);

void density_matrix(int sect_h, int sect_p, int rep1, int state1, int rep2, int state2);

void dipole_length_tdms(int sect_h, int sect_p);

void model_space_property(cc_ms_prop_query_t *prop_query);

#endif /* CC_HEFF_H_INCLUDED */
