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
#include <string.h>

#include "memory.h"


/**
 * Allocates space for a matrix of a given size (n x m); sets all matrix
 * elements to zero. Type of the matrix is given by the third argument.
 *
 * @param type CC_DOUBLE or CC_COMPLEX
 * @param n number of rows
 * @param m number of columns
 * @return pointer to the allocated matrix
 *
 * @note matrix must be manually deallocated by the user
 */
void *x_zeros(data_type_t type, size_t n, size_t m)
{
    assert(type == CC_DOUBLE || type == CC_DOUBLE_COMPLEX);

    void *p = cc_malloc(n * m * type);
    memset(p, 0, n * m * type);

    return p;
}


double *d_zeros(size_t n, size_t m)
{
    return (double *) x_zeros(CC_DOUBLE, n, m);
}


double complex *z_zeros(size_t n, size_t m)
{
    return (double complex *) x_zeros(CC_DOUBLE_COMPLEX, n, m);
}
