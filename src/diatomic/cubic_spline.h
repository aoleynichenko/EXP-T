/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2024 The EXP-T developers.
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

#ifndef EXPT_CUBIC_SPLINE_H
#define EXPT_CUBIC_SPLINE_H

typedef struct {
    int n;
    double *x;
    double *y;
    double *y2;
} cubic_spline_t;


cubic_spline_t *construct_cubic_spline(int n, double *x, double *y);

double evaluate_spline(cubic_spline_t *spline, double x);

void find_spline_minimum(cubic_spline_t *spline, double *xmin, double *ymin, int print_level);

double spline_first_derivative(cubic_spline_t *spline, double x);

double spline_second_derivative(cubic_spline_t *spline, double x);

double spline_third_derivative(cubic_spline_t *spline, double x);

void delete_spline(cubic_spline_t *spline);

#endif //EXPT_CUBIC_SPLINE_H
