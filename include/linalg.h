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

/*******************************************************************************
* linalg.h
* ========
*
* Interface to linear algebra packages (BLAS/LAPACK) and some related tools for
* manipulating with linear algebra objects (vectors, matrices, etc).
* See linalg.c for the implementation.
*
* 2018-2021 Alexander Oleynichenko
******************************************************************************/

#ifndef CC_LINALG_H_INCLUDED
#define CC_LINALG_H_INCLUDED

#include <complex.h>
#include <stddef.h>

#include "platform.h"

typedef enum {
    CC_FLOAT = sizeof(float),
    CC_DOUBLE = sizeof(double),
    CC_REAL = CC_DOUBLE,
    CC_FLOAT_COMPLEX = sizeof(float complex),
    CC_DOUBLE_COMPLEX = sizeof(double complex),
    CC_COMPLEX = CC_DOUBLE_COMPLEX
} data_type_t;

// absolute value
double xabs(data_type_t type, void *x);

// transposition/copying/conjugation of matrices
void xmatcopy(data_type_t type, char op, size_t rows, size_t cols, void *A, void *B);
void xtranspose(data_type_t type, size_t rows, size_t cols, void *A);

// a constant times a vector plus a vector
void xaxpy(data_type_t type, size_t n,
           double alpha, const void *x, void *y);

// linear combination of two vectors: z = alpha*x+beta*y
void xaxpby(data_type_t type, size_t n, const void *alpha,
            const void *x, const void *beta, const void *y, void *z);

size_t ixadiffmax(data_type_t type, size_t n, const void *a, const void *b, double *diffmax);

size_t ixamax(data_type_t type, size_t n, const void *a, double *max_val);

// scalar product
double complex xdot(data_type_t type, char *conjx, char *conjy, size_t n, const void *x, const void *y);

// wrapper for matrix multiplications
void xgemm(data_type_t data_type, char *trans_a, char *trans_b,
           int m, int n, int k, void *alpha, void *A, int lda, void *B, int ldb,
           void *beta, void *C, int ldc);

// print matrix
void xprimat(data_type_t type, const void *A, int n, int m, const char *comment);

// return a new matrix of given type, filled with zeros
void *xzeros(data_type_t type, size_t n, size_t m);

double *dzeros(size_t n, size_t m);

double complex *zzeros(size_t n, size_t m);

// solves eigenvalue problem for matrix A
int eig(int n, double complex *A, double complex *ev, double complex *vl, double complex *vr);

// singular value decomposition (SVD)
void svd(int n, double complex *A, double *lambda, double complex *U, double complex *V);

// overlap matrix between two sets of vectors
void overlap(size_t n, double complex *C1, double complex *C2, double complex *S);

// symmetric (Loewdin) orthogonalization
void loewdin_orth(size_t n, double complex *C, double complex *C_orth, int print);

// inverse matrix
int inv(size_t n, double complex *A, double complex *Ainv);

#endif /* CC_LINALG_H_INCLUDED */
