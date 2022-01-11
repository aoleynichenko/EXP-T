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

#include "numerov.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mkl.h"

#include "cubic_spline.h"
#include "input_data.h"
#include "units.h"
#include "utils.h"


#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

double *new_zero_matrix(int n);

double *new_tridiagonal_matrix(int n, double diag_values, double offdiag_values);

void print_matrix(int n, double *A, char *comment);

int inverse_matrix(size_t n, double *A, double *Ainv);

void fd2_matrix_solver(input_data_t *input_data, cubic_spline_t *pot, int J, double emin, double emax, int *nroots, double **eigenvalues, double **wavefunctions);


void numerov_matrix_solver(input_data_t *input_data, cubic_spline_t *pot, int J, double emin, double emax, int *nroots, double **eigenvalues, double **wavefunctions)
{
    int N = input_data->grid_size;
    double rmin = pot->x[0];
    double rmax = pot->x[pot->n - 1];
    double l = rmax - rmin;
    double h = 1.0 / (N + 1);
    double mu = input_data->reduced_mass;
    double gamma2 = 2 * mu * l * l;

    double *A = new_tridiagonal_matrix(N, -2.0 / (h * h), 1.0 / (h * h));
    double *B = new_tridiagonal_matrix(N, 10.0 / 12.0, 1.0 / 12.0);
    double *V = new_zero_matrix(N);

    /*
     * construct potential energy matrix
     */
    for (int i = 0; i < N; i++) {
        double ri = rmin + (i + 1) * h * l;
        double E_rot = 1.0 / (2.0 * mu * ri * ri) * J * (J + 1);
        V[i * N + i] = evaluate_spline(pot, ri) + E_rot;
    }

    /*
     * construct kinetic energy matrix, then hamiltonian
     * V := -1/gamma2 B^-1 A + V
     */
    double *Binv = new_zero_matrix(N);
    inverse_matrix(N, B, Binv);
    cblas_dgemm(CblasRowMajor, CblasNoTrans,  CblasNoTrans, N, N, N, -1.0/gamma2, Binv, N, A, N, 1.0, V, N);

    /*
     * find eigenvalues and eigenvectors
     */
    double *w = (double *) calloc(N, sizeof(double));
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', N, V, N, w);

    /*
     * calculate number of roots in the given energy range [emin,emax]
     * and allocate memory for eigenvalues and eigenvectors
     */

    *nroots = 0;
    for (int i = 0; i < N; i++) {
        if (emin <= w[i] && w[i] <= emax) {
            *nroots = *nroots + 1;
        }
    }
    *eigenvalues = (double *) calloc(*nroots, sizeof(double));
    *wavefunctions = (double *) calloc((*nroots) * N, sizeof(double));

    /*
     * copy eigenvalues & eigenvectors to the output arrays
     */
    int root_count = 0;
    for (int i = 0; i < N; i++) {
        if (emin <= w[i] && w[i] <= emax) {

            (*eigenvalues)[root_count] = w[i];

            for (int j = 0; j < N; j++) {
                (*wavefunctions)[root_count * N + j] = V[j * N + root_count];
            }

            root_count++;
        }
    }

    // clean up
    free(A);
    free(B);
    free(Binv);
    free(V);
    free(w);
}


/*
 * allocates square matrix filled with zeros
 */
double *new_zero_matrix(int n)
{
    return (double *) calloc(n * n, sizeof(double));
}


/*
 * allocates square tridiagonal n x n matrix.
 *
 * elements on the main diagonal are equal to 'diag_values',
 * and elements on the 1st and -1st diagonals are equal to 'offdiag_values'.
 * other elements are zero.
 */
double *new_tridiagonal_matrix(int n, double diag_values, double offdiag_values)
{
    double *A = (double *) calloc(n * n, sizeof(double));

    for (int i = 0; i < n; i++) {
        A[i * n + i] = diag_values;
    }
    for (int i = 1; i < n; i++) {
        A[(i - 1) * n + i] = offdiag_values;
        A[i * n + (i - 1)] = offdiag_values;
    }

    return A;
}


/*
 * prints n x n matrix to stdout
 */
void print_matrix(int n, double *A, char *comment)
{
    printf(" Matrix %dx%d: %s\n", n, n, comment);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.6f", A[i * n + j]);
        }
        printf("\n");
    }
}


/**
 * calculates inverse matrix
 * @param n matrix dimensino
 * @param A square n x n matrix; will not be destroyed
 * @param Ainv inverse matrix Ainv = A^{-1}
 */
int inverse_matrix(size_t n, double *A, double *Ainv)
{
    int *ipiv;
    int ret;

    ipiv = (int *) calloc(n + 1, sizeof(int));

    memmove(Ainv, A, sizeof(double) * n * n);

    ret = LAPACKE_dgetrf(CblasRowMajor, n, n, Ainv, n, ipiv);

    if (ret != 0) {
        free(ipiv);
        return ret;
    }

    ret = LAPACKE_dgetri(CblasRowMajor, n, Ainv, n, ipiv);

    free(ipiv);

    return ret;
}