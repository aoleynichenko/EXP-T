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

/*
 * Operations with formatted files containing effective Hamiltonians.
 * The format of formatted HEFFF files was primarily introduced by A. Zaitsevskii.
 */

#ifndef CC_FORMATTED_HEFF_FILE_H_INCLUDED
#define CC_FORMATTED_HEFF_FILE_H_INCLUDED

#include <complex.h>
#include <stdio.h>

FILE *open_formatted_heff_file(int sect_h, int sect_p, char *label);

void close_formatted_heff_file(FILE *hefff);

void formatted_heff_file_write_block(FILE *hefff, int carith, int rep_no, size_t dim, double complex *heff);

size_t first_nonzero_irrep(size_t *block_dims);

#endif // CC_FORMATTED_HEFF_FILE_H_INCLUDED
