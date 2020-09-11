/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2020 The EXP-T developers.
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

#include <stdio.h>
#include <string.h>

#include "memory.h"


/**
 * calculates inverse matrix
 * @param n matrix dimensino
 * @param A square n x n matrix; will not be destroyed
 * @param Ainv inverse matrix Ainv = A^{-1}
 */
int inv(size_t n, double complex *A, double complex *Ainv)
{
    int *ipiv;
    int ret;

    ipiv = (int *) cc_malloc(sizeof(int) * (n+1));

    memmove(Ainv, A, sizeof(double complex) * n * n);

#ifdef BLAS_MKL
    ret =  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, (MKL_Complex16 *) Ainv, n, ipiv);
#else
    ret =  LAPACKE_zgetrf(CblasRowMajor, n, n, Ainv, n, ipiv);
#endif

    if (ret !=0) {
        return ret;
    }

#ifdef BLAS_MKL
    ret = LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, (MKL_Complex16 *) Ainv, n, ipiv);
#else
    ret = LAPACKE_zgetri(CblasRowMajor, n, Ainv, n, ipiv);
#endif

    return ret;
}
