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

#include "linalg.h"

#include <assert.h>
#include <math.h>

size_t idadiffmax(size_t n, const double *a, const double *b, double *diffmax);

size_t izadiffmax(size_t n, const double complex *a, const double complex *b, double *diffmax);


/**
 * Finds the index of the pair of elements with maximum absolute value
 * of difference between them.
 *
 * @param type CC_DOUBLE or CC_DOUBLE_COMPLEX
 * @param n    specifies the number of elements in vectors
 * @param a    vector 1
 * @param b    vector 2
 */
size_t ixadiffmax(data_type_t type, size_t n, const void *a, const void *b, double *diffmax)
{
    assert(type == CC_DOUBLE || type == CC_DOUBLE_COMPLEX);

    if (type == CC_DOUBLE) {
        return idadiffmax(n, (const double *) a, (const double *) b, diffmax);
    }
    else { // type == CC_DOUBLE_COMPLEX
        return izadiffmax(n, (const double complex *) a, (const double complex *) b, diffmax);
    }
}


/**
 * version for real numbers (double)
 */
size_t idadiffmax(size_t n, const double *a, const double *b, double *diffmax)
{
    size_t max_idx = 0;
    double max_val = 0.0;

    for (size_t i = 0; i < n; i++) {
        double absdiff = fabs(a[i] - b[i]);
        if (absdiff > max_val) {
            max_idx = i;
            max_val = absdiff;
        }
    }

    *diffmax = max_val;
    return max_idx;
}


/**
 * version for complex numbers (double complex)
 */
size_t izadiffmax(size_t n, const double complex *a, const double complex *b, double *diffmax)
{
    size_t max_idx = 0;
    double max_val = 0.0;

    for (size_t i = 0; i < n; i++) {
        double absdiff = cabs(a[i] - b[i]);
        if (absdiff > max_val) {
            max_idx = i;
            max_val = absdiff;
        }
    }

    *diffmax = max_val;
    return max_idx;
}
