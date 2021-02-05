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

/**
* Utility functions for FSCC models.
*
* 2020-2021 Alexander Oleynichenko
*/

#ifndef CCUTILS_H_INCLUDED
#define CCUTILS_H_INCLUDED

void flush_block_sizes(char *diag_name, char *file_name);

void damping(int h, int p, char *old_ampl, char *new_ampl, int iter);

void print_max_ampl(char *label, char *dg_name);

double t1_diagnostic(char *dg_name, int nelec);

int check_divergence(char *dg_name);

void print_ampl_norm(char *label, char *dg_name);

#endif /* CCUTILS_H_INCLUDED */
