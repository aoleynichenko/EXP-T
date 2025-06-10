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

size_t my_idamax(size_t n, const double *a, double *max_val);

size_t my_izamax(size_t n, const double complex *a, double *max_val);


/**
 * Finds the index of the element with maximum absolute value.
 *
 * @param type CC_DOUBLE or CC_DOUBLE_COMPLEX
 * @param n    specifies the number of elements in vector
 * @param a    vector
 * @param max_val max absolute value of the element found
 * @return index, max_val
 */
size_t ixamax(data_type_t type, size_t n, const void *a, double *max_val)
{
    assert(type == CC_DOUBLE || type == CC_DOUBLE_COMPLEX);

    if (type == CC_DOUBLE) {
        return my_idamax(n, (const double *) a, max_val);
    }
    else { // type == CC_DOUBLE_COMPLEX
        return my_izamax(n, (const double complex *) a, max_val);
    }
}


/**
 * version for real numbers
 */
size_t my_idamax(size_t n, const double *a, double *max_val)
{
    double mx = 0.0;
    size_t max_idx = 0;

    for (size_t i = 0; i < n; i++) {
        if (fabs(a[i]) > mx) {
            mx = fabs(a[i]);
            max_idx = i;
        }
    }

    *max_val = mx;
    return max_idx;
}


/**
 * version for complex numbers
 */
size_t my_izamax(size_t n, const double complex *a, double *max_val)
{
    double mx = 0.0;
    size_t max_idx = 0;

    for (size_t i = 0; i < n; i++) {
        if (cabs(a[i]) > mx) {
            mx = cabs(a[i]);
            max_idx = i;
        }
    }

    *max_val = mx;
    return max_idx;
}
