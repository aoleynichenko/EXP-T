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

#include "linalg.h"

#include <assert.h>
#include <ctype.h>

double my_ddot(size_t n, const double *x, const double *y);

double complex my_zdot(char *conjx, char *conjy, size_t n, const double complex *zx, const double complex *zy);


/**
 * Scalar product of two vectors:
 *   xdot = (x^conjx, y^conjy)
 * 'X' prefix means that this subroutine is designed to work both with double
 * and complex numbers.
 *
 * @param type   data type used (CC_DOUBLE or CC_DOUBLE COMPLEX)
 * @param conjx  if == "C", the 'x' vector will be complex conjugated
 * @param conjy  if == "C", the 'y' vector will be complex conjugated
 * @param n      number of elements in input vector(s)
 * @param x      real or complex array, dimension == n
 * @param y      real or complex array, dimension == n
 * @return scalar product
 */

double complex xdot(data_type_t type, char *conjx, char *conjy, size_t n, const void *x, const void *y)
{
    assert(type == CC_DOUBLE || type == CC_DOUBLE_COMPLEX);

    if (type == CC_DOUBLE) {
        return my_ddot(n, x, y) + 0.0 * I;
    }
    else { // CC_DOUBLE_COMPLEX
        return my_zdot(conjx, conjy, n, x, y);
    }
}


double my_ddot(size_t n, const double *x, const double *y)
{
    double prod = 0.0;

    for (size_t i = 0; i < n; i++) {
        prod += x[i] * y[i];
    }

    return prod;
}


double complex my_zdot(char *conjx, char *conjy, size_t n, const double complex *zx, const double complex *zy)
{
    int conjugate_x = tolower(conjx[0]) == 'c';
    int conjugate_y = tolower(conjy[0]) == 'c';

    double complex prod = 0.0 + 0.0 * I;

    if (!conjugate_x && !conjugate_y) {
        for (size_t i = 0; i < n; i++) {
            prod += zx[i] * zy[i];
        }
    }
    else if (conjugate_x && !conjugate_y) {
        for (size_t i = 0; i < n; i++) {
            prod += conj(zx[i]) * zy[i];
        }
    }
    else if (!conjugate_x && conjugate_y) {
        for (size_t i = 0; i < n; i++) {
            prod += zx[i] * conj(zy[i]);
        }
    }
    else {  // 'C' && 'C'
        for (size_t i = 0; i < n; i++) {
            prod += conj(zx[i]) * conj(zy[i]);
        }
    }

    return prod;
}