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

#include "solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "input_data.h"
#include "finite_diff.h"
#include "matrix_element.h"
#include "numerov.h"
#include "units.h"
#include "utils.h"


/**
 * Provides interface to the solvers of the radial Schrodinger equation.
 * Currently implemented solvers are:
 * - Numerov matrix solver
 * - Finite difference 2nd order matrix solver
 *
 * Returns arrays 'nroots', 'energies' and 'wavefunctions'
 * (these arrays are allocated inside the function).
 * Energies and wavefunctions are obtained for J = Jmin ... Jmax.
 * Number of roots obtained for the given J is stored in (*nroots)[J].
 * Energy of the (J,v) state is stored in (*energies)[J][v].
 * Wavefunction values at each point of the radial grid are stored
 * as the one-dimensional array (*wavefunctions)[J][v].
 */
void solve(input_data_t *input_data, cubic_spline_t *pot, int **nroots, double ***energies, double ****wavefunctions)
{
    int v_min = input_data->v_min;
    int v_max = input_data->v_max;
    int J_min = input_data->J_min;
    int J_max = input_data->J_max;
    int n_rot = J_max - J_min + 1;
    int n_vib = v_max - v_min + 1;

    double re, emin;
    find_spline_minimum(pot, &re, &emin, 0);
    double emax = MIN(pot->y[0], pot->y[pot->n-1]);

    int ngrid = input_data->grid_size;

    *nroots = (int *) calloc(n_rot, sizeof(int));
    *energies = new_2d_array(n_rot, n_vib);
    *wavefunctions = new_3d_array(n_rot, n_vib, ngrid);

    for (int J = J_min; J <= J_max; J++) {

        int index_J = J - J_min;

        int nroots_J = 0;
        double *E_J = NULL;
        double *psi_J = NULL;

        if (input_data->solver == SOLVER_NUMEROV) {
            numerov_matrix_solver(input_data, pot, J, emin, emax, &nroots_J, &E_J, &psi_J);
        }
        else {
            // 2nd order central finite difference method
            fd2_matrix_solver(input_data, pot, J, emin, emax, &nroots_J, &E_J, &psi_J);
        }

        for (int v = v_min; v <= MIN(v_max,nroots_J-1); v++) {
            int index_v = v - v_min;

            // save number of vibrational states
            (*nroots)[index_J] = nroots_J;

            // save energy
            (*energies)[index_J][index_v] = E_J[v];

            // save wavefunction
            double *target_psi = (*wavefunctions)[index_J][index_v];
            double *source_psi = psi_J + ngrid * v;
            memcpy(target_psi, source_psi, sizeof(double) * ngrid);
        }

        free(E_J);
        free(psi_J);
    }
}


/**
 * prints beautiful table with energies of the J-v states
 * and expectation values of properties
 */
void print_energy_levels(input_data_t *input_data, cubic_spline_t *pot, int *nroots, double **energies,
                         double ***wavefunctions)
{
    int v_min = input_data->v_min;
    int v_max = input_data->v_max;
    int J_min = input_data->J_min;
    int J_max = input_data->J_max;

    double r_min = pot->x[0];
    double r_max = pot->x[pot->n - 1];
    int n_rot = J_max - J_min + 1;
    int n_vib = v_max - v_min + 1;
    int ngrid = input_data->grid_size;
    double emin = input_data->min_energy;
    double *radial_grid = construct_radial_grid(ngrid, r_min, r_max, input_data->mapping);

    cubic_spline_t *r_spline = construct_cubic_spline(input_data->n_points, input_data->r, input_data->r);

    printf(" > rovibrational energy levels\n\n");
    printf("     J   v       energy, cm^-1        < r >, A        < prop >\n\n");

    for (int J = J_min; J <= J_max; J++) {
        int index_J = J - J_min;

        for (int v = v_min; v <= MIN(v_max, nroots[index_J]-1); v++) {
            int index_v = v - v_min;

            double E_Jv = energies[index_J][index_v];
            double *psi = wavefunctions[index_J][index_v];

            double e_cm = (E_Jv - emin) * ATOMIC_TO_CM;
            double r_v = matrix_element_spline(ngrid, radial_grid, psi, psi, r_spline);
            double prp = 0.0;
            if (input_data->prop != NULL) {
                prp = matrix_element_spline(ngrid, radial_grid, psi, psi, input_data->prop);
            }

            printf("@ %4d%4d%20.4f%16.6f%16.6f\n", J, v, e_cm, r_v * ATOMIC_TO_ANGSTROM, prp);
        }
        printf("\n");
    }

    free(radial_grid);
    delete_spline(r_spline);
}
