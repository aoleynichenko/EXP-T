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

#include "finite_diff.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blas_lapack.h"
#include "cubic_spline.h"
#include "input_data.h"


/**
 * matrix solver of radial Schrodinger equation
 * based on the 2nd order formula for numerical differentiation.
 *
 * See also the paper:
 * V. V. Meshkov, A. V. Stolyarov, R. J. Le Roy
 * Adaptive analytical mapping procedure for efficiently solving the radial Schrodinger equation
 * Phys. Rev. A 78, 052510 (2008)
 * doi: 10.1103/PhysRevA.78.052510
 */
void fd2_matrix_solver(input_data_t *input_data, cubic_spline_t *pot, int J, double emin, double emax, int *nroots, double **eigenvalues, double **wavefunctions)
{
    int N = input_data->grid_size;
    double rmin = pot->x[0];
    double rmax = pot->x[pot->n - 1];
    double l = rmax - rmin;
    double mu = input_data->reduced_mass;

    mapping_t *map = input_data->mapping;
    double ymin = map->r2y(map, rmin);
    double ymax = map->r2y(map, rmax);
    double h = (ymax - ymin) / (N + 1);

    /*
     * precalculation of radial grid, effective potential and the dr/dy function
     */
    double *r = (double *) calloc(N, sizeof(double));
    double *g = (double *) calloc(N, sizeof(double));
    double *U = (double *) calloc(N, sizeof(double));

    for (int i = 0; i < N; i++) {

        double yi = ymin + (i + 1) * h;
        r[i] = map->y2r(map, yi);

        double F = map->F(map, yi);
        g[i] = map->drdy(map, yi);

        double E_rot = 1.0 / (2.0 * mu * r[i] * r[i]) * J * (J + 1);
        U[i] = evaluate_spline(pot, r[i]) + E_rot - 1.0/(2 * mu) * F / (g[i] * g[i]);
    }

    /*
     * construct Hamiltonian
     */
    double *diagonal = (double *) calloc(N, sizeof(double));
    double *off_diag = (double *) calloc(N, sizeof(double));

    for (int i = 0; i < N; i++) {
        diagonal[i] = 2.0 / (2.0 * mu * h * h * g[i] * g[i]) + U[i];
    }
    for (int i = 1; i < N; i++) {
        off_diag[i-1] = -1.0 / (2.0 * mu * h * h * g[i-1] * g[i]);
    }

    /*
     * find eigenvalues and eigenvectors.
     * eigenvalues will be stored in the 'diagonal' array.
     * we declare the column-major order to obtain the eigenvectors stored row-wise.
     *
     */
    double *eigenvec = (double *) calloc(N * N, sizeof(double));
    LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', N, diagonal, off_diag, eigenvec, N);

    /*
     * calculate number of roots in the given energy range [emin,emax]
     * and allocate memory for eigenvalues and eigenvectors
     */

    *nroots = 0;
    for (int i = 0; i < N; i++) {
        if (emin <= diagonal[i] && diagonal[i] <= emax) {
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
        if (emin <= diagonal[i] && diagonal[i] <= emax) {

            (*eigenvalues)[root_count] = diagonal[i];
            int offset = root_count * N;
            memcpy(*wavefunctions + offset, eigenvec + offset, sizeof(double) * N);

            root_count++;
        }
    }

    // clean up
    free(r);
    free(g);
    free(U);
    free(eigenvec);
    free(diagonal);
    free(off_diag);
}
