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

#ifndef CC_NATURAL_SPINORS_H_INCLUDED
#define CC_NATURAL_SPINORS_H_INCLUDED

#include <stdio.h>
#include <complex.h>

void calculate_natural_spinors_and_occ_numbers(char *path_to_dm_file, char *path_to_nat_spinors_file);

void calculate_natural_transition_spinors_and_singular_values(
    char *path_to_tran_dm_file,
    char *path_to_nat_spinors_from_file,
    char *path_to_nat_spinors_to_file
);

void print_natural_spinor(FILE *out, double complex *nat_spinor_coefs, double thresh);

void write_natural_spinors_formatted(
    char *natorb_file_name,
    double *occ_numbers, double complex *nat_orbs,
    double occ_thresh
);

#endif // CC_NATURAL_SPINORS_H_INCLUDED
