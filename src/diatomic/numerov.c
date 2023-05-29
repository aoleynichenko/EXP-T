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

#include <stdlib.h>
#include <string.h>

#include "blas_lapack.h"
#include "cubic_spline.h"
#include "input_data.h"
#include "utils.h"


int inverse_matrix(size_t n, double *A, double *Ainv);


/**
 * Matrix Numerov method for integration of the radial Schrodinger equation.
 * Uniform grid is used.
 *
 * For details on the method, see, for example:
 *
 * M. Pillai, J. Goglio, T. G. Walker,
 * Matrix Numerov method for solving Schrodinger’s equation.
 * Am. J. Phys. 80, 1017 (2012)
 * doi: 10.1119/1.4748813
 */
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
     * find eigenvalues and eigenvectors.
     * we declare the column-major order to obtain the eigenvectors stored row-wise.
     */
    double *w = (double *) calloc(N, sizeof(double));
    LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', N, V, N, w);

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

            int offset = root_count * N;
            memcpy(*wavefunctions + offset, V + offset, sizeof(double) * N);

            root_count++;
        }
    }

    /*
     * clean up
     */
    free(A);
    free(B);
    free(Binv);
    free(V);
    free(w);
}


/**
 * calculates inverse matrix
 * @param n matrix dimension
 * @param A square n x n matrix; will not be destroyed
 * @param Ainv inverse matrix Ainv = A^{-1}
 */
int inverse_matrix(size_t n, double *A, double *Ainv)
{
    int *ipiv;
    int ret;

    ipiv = (int *) calloc(n + 1, sizeof(int));

    memmove(Ainv, A, sizeof(double) * n * n);

    ret = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, Ainv, n, ipiv);

    if (ret != 0) {
        free(ipiv);
        return ret;
    }

    ret = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, Ainv, n, ipiv);

    free(ipiv);

    return ret;
}



