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

#ifndef EXPT_UTILS_H
#define EXPT_UTILS_H

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

double *new_double_array(int n, double *values);

double **new_2d_array(int n, int m);

void delete_2d_array(double **A, int n, int m);

double ***new_3d_array(int n, int m, int k);

void delete_3d_array(double ***A, int n, int m, int k);

void rescale_array(int len, double *array, double factor);

double *new_zero_matrix(int n);

double *new_tridiagonal_matrix(int n, double diag_values, double offdiag_values);

void print_matrix(int n, double *A, char *comment);

#endif //EXPT_UTILS_H
