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

#include "linalg.h"

#include <assert.h>
#include <string.h>

#include "memory.h"


/**
 * Allocates space for a "n x n" square matrix; sets all diagonal elements
 * matrix to 1.0. Type of the matrix is given by the third argument.
 *
 * @param type   CC_DOUBLE or CC_COMPLEX
 * @param n      square matrix dimension
 * @return pointer to the allocated matrix
 *
 * @note matrix must be manually deallocated by the user
 */
void *x_identity(data_type_t type, size_t n)
{
    assert(type == CC_DOUBLE || type == CC_DOUBLE_COMPLEX);

    if (type == CC_DOUBLE) {
        return d_identity(n);
    }
    else {  // type == CC_DOUBLE_COMPLEX)
        return z_identity(n);
    }
}


double *d_identity(size_t n)
{
    double *identity = (double *) x_zeros(CC_DOUBLE, n, n);

    for (size_t i = 0; i < n; i++) {
        identity[n * i + i] = 1.0;
    }

    return identity;
}


double complex *z_identity(size_t n)
{
    double complex *identity = (double complex *) x_zeros(CC_DOUBLE_COMPLEX, n, n);

    for (size_t i = 0; i < n; i++) {
        identity[n * i + i] = 1.0 + 0.0*I;
    }

    return identity;
}
