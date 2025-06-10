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

#ifndef CC_EIGENVALUES_H_INCLUDED
#define CC_EIGENVALUES_H_INCLUDED

#include <complex.h>
#include <stdlib.h>

#include "slater_det.h"

double get_lowest_eigenvalue(const size_t *block_dims, double complex **eigvalues);

void print_eigenvalues_table(
        int sect_h,
        int sect_p,
        size_t *block_dims,
        slater_det_t **det_list,
        double complex **eigvalues,
        double degen_thresh,
        double complex **coef_right
);

size_t first_nonzero_irrep(const size_t *block_dims);

int get_nroots_for_irrep(char *irrep_name);

void get_nroots(size_t *block_dims, double complex **eigvalues, size_t *nroots);

size_t get_nroots_from_cutoff(size_t irrep_dim, double complex *eigvalues, double emin, double cutoff);

#endif // CC_EIGENVALUES_H_INCLUDED
