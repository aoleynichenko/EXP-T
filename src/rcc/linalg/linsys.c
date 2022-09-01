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

#include <stdlib.h>
#include <string.h>

#include "cblas_lapacke.h"


/**
 * Solves system of linear equations:
 * Ax = b
 * On exit, A and b arrays are preserved.
 *
 * Returns:
 * = 0:  successful exit
 * < 0:  if INFO = -i, the i-th argument had an illegal value
 * > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
 *       has been completed, but the factor U is exactly
 *       singular, so the solution could not be computed.
 */
int linsys(int n, double *A, double *b, double *x)
{
    int *ipiv = (int *) malloc(sizeof(int) * n);
    double *A_copy = (double *) malloc(sizeof(double) * n * n);

    memcpy(A_copy, A, sizeof(double) * n * n);
    memcpy(x, b, sizeof(double) * n);

    int info = LAPACKE_dgesv(CblasRowMajor, n, 1, A_copy, n, ipiv, x, 1);

    free(A_copy);
    free(ipiv);

    return info;
}
