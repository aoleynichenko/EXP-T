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

void my_daxpy(size_t n, double alpha, const double *x, double *y);

void my_zaxpy(size_t n, double alpha,
              const double complex *x, double complex *y);

void my_daxpby(size_t n, double alpha, const double *x,
               double beta, const double *y, double *z);

void my_zaxpby(size_t n, double complex alpha, const double complex *x,
               double complex beta, const double complex *y, double complex *z);


/**
 * Computes a constant times a vector plus a vector:
 *   y := alpha*x + y
 *
 * @param type  CC_DOUBLE or CC_DOUBLE_COMPLEX
 * @param n     number of elements in the vectors
 * @param alpha scaling factor for the values in X
 * @param x     input vector X
 * @param y     input and output vector Y
 */
void xaxpy(data_type_t type, size_t n,
           double alpha, const void *x, void *y)
{
    assert(type == CC_DOUBLE || type == CC_DOUBLE_COMPLEX);

    if (type == CC_DOUBLE) {
        my_daxpy(n, alpha, (double *) x, (double *) y);
    }
    else { // CC_DOUBLE_COMPLEX
        my_zaxpy(n, alpha, (double complex *) x, (double complex *) y);
    }
}


void my_daxpy(size_t n, double alpha, const double *x, double *y)
{
    for (size_t i = 0; i < n; i++) {
        y[i] = alpha * x[i] + y[i];
    }
    // void cblas_daxpy(const int N, const double alpha, const double *X,
    //                 const int incX, double *Y, const int incY);
    //cblas_daxpy(n, alpha, x, 1, y, 1);
}


void my_zaxpy(size_t n, double alpha,
              const double complex *x, double complex *y)
{
    for (size_t i = 0; i < n; i++) {
        y[i] = alpha * x[i] + y[i];
    }
    //double complex a = alpha + 0.0*I;
    // void cblas_zaxpy(const int __N, const void *__alpha, const void *__X, const int __incX, void *__Y, const int __incY);
    //cblas_zaxpy(n, &a, x, 1, y, 1);
}


/**
 * Computes a constant times a vector plus a constant times a vector
 * z := alpha * x + beta * y
 *
 * @param type  CC_DOUBLE or CC_DOUBLE_COMPLEX
 * @param n     number of elements in the vectors
 * @param alpha scaling factor for the values in X
 * @param x     input vector X
 * @param beta  scaling factor for the values in Y
 * @param y     input vector Y
 * @param z     output vector
 */
void xaxpby(data_type_t type, size_t n, const void *alpha,
            const void *x, const void *beta, const void *y, void *z)
{
    assert(type == CC_DOUBLE || type == CC_DOUBLE_COMPLEX);

    if (type == CC_DOUBLE) {
        double a = *((double *) alpha);
        double b = *((double *) beta);
        my_daxpby(n, a, (double *) x, b, (double *) y, (double *) z);
    }
    else { // CC_DOUBLE_COMPLEX
        double complex a = *((double complex *) alpha);
        double complex b = *((double complex *) beta);
        my_zaxpby(n, a, (double complex *) x, b, (double complex *) y, (double complex *) z);
    }
}


void my_daxpby(size_t n, double alpha, const double *x,
               double beta, const double *y, double *z)
{
    for (size_t i = 0; i < n; i++) {
        z[i] = alpha * x[i] + beta * y[i];
    }
}


void my_zaxpby(size_t n, double complex alpha, const double complex *x,
               double complex beta, const double complex *y, double complex *z)
{
    for (size_t i = 0; i < n; i++) {
        z[i] = alpha * x[i] + beta * y[i];
    }
}
