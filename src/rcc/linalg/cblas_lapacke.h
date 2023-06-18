/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2023 The EXP-T developers.
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

/*
 * this header file provides interface to the linear algebra routines
 * implemented in BLAS/LAPACK and the abstracted by the CBLAS/LAPACKE
 * interface written in C.
 * can be used with both MKL and OpenBLAS.
 */

#ifndef EXPT_BLAS_LAPACK_H
#define EXPT_BLAS_LAPACK_H

typedef enum {CblasRowMajor=101, CblasColMajor=102} CBLAS_LAYOUT;
typedef enum {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113} CBLAS_TRANSPOSE;


/*
 * multiplication of real general matrices
 */
void cblas_dgemm(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

/*
 * multiplication of complex general matrices
 */
void cblas_zgemm(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc);

/*
 * LAPACKE matrix layouts
 */
#define LAPACK_ROW_MAJOR               101
#define LAPACK_COL_MAJOR               102

/*
 * LAPACKE data types
 */
#ifndef lapack_int
#define lapack_int     int
#endif

#ifndef lapack_complex_double
#include <complex.h>
#define lapack_complex_double   double _Complex
#endif

/*
 * eigensolver for real symmetric matrices
 */
lapack_int LAPACKE_dsyev( int matrix_layout, char jobz, char uplo, lapack_int n,
                          double* a, lapack_int lda, double* w );

/*
 * eigensolver for real symmetric tridiagonal matrices
 */
lapack_int LAPACKE_dstev( int matrix_layout, char jobz, lapack_int n, double* d,
                          double* e, double* z, lapack_int ldz );

/*
 * eigensolver for complex general matrices
 */
/*extern void zgeev_(char *jobvl, char *jobvr, int *n, double complex *a,
                   int *lda, double complex *w, double complex *vl, int *ldvl,
                   double complex *vr, int *ldvr, double complex *work, int *lwork,
                   double *rwork, int *info);*/
lapack_int LAPACKE_zgeev( int matrix_layout, char jobvl, char jobvr,
                          lapack_int n, lapack_complex_double* a,
                          lapack_int lda, lapack_complex_double* w,
                          lapack_complex_double* vl, lapack_int ldvl,
                          lapack_complex_double* vr, lapack_int ldvr );

/*
 * solves the system of linear equations A * X = B for general real matrices
 */
lapack_int LAPACKE_dgesv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          double* a, lapack_int lda, lapack_int* ipiv,
                          double* b, lapack_int ldb );

/*
 * LU factorization of a real general M-by-N matrix A
 */
lapack_int LAPACKE_dgetrf( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, lapack_int* ipiv );

/*
 * computes the inverse of a matrix using the LU factorization computed by DGETRF
 */
lapack_int LAPACKE_dgetri( int matrix_layout, lapack_int n, double* a,
                           lapack_int lda, const lapack_int* ipiv );

/*
 * LU factorization of a complex general M-by-N matrix A
 */
lapack_int LAPACKE_zgetrf( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_int* ipiv );

/*
 * computes the inverse of a matrix using the LU factorization computed by ZGETRF
 */
lapack_int LAPACKE_zgetri( int matrix_layout, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           const lapack_int* ipiv );

#endif //EXPT_BLAS_LAPACK_H
