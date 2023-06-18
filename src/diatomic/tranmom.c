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


#include "tranmom.h"

#include <stdio.h>
#include <stdlib.h>

#include "units.h"
#include "matrix_element.h"

#define MAX(a, b) (((a)>(b))?(a):(b))
#define MIN(a, b) (((a)<(b))?(a):(b))


/**
 * calculation of transition moments < psi_i | prop | psi_j >.
 * prints table with transition moments.
 */
void calc_transition_moments(input_data_t *input_data,
                             int *nroots1, double **energies1, double ***wavefunctions1,
                             int *nroots2, double **energies2, double ***wavefunctions2)
{
    int v_min = input_data->v_min;
    int v_max = input_data->v_max;
    int J_min = input_data->J_min;
    int J_max = input_data->J_max;
    double min_energy = input_data->min_energy;
    int n_grid = input_data->grid_size;

    cubic_spline_t *pot1 = input_data->pot1;
    cubic_spline_t *pot2 = input_data->pot2;
    cubic_spline_t *prop_fun = input_data->prop;

    double r_min = pot1->x[0];
    double r_max = pot1->x[pot1->n - 1];
    double *radial_grid = construct_radial_grid(n_grid, r_min, r_max, input_data->mapping);
    mapping_t *mapping = input_data->mapping;

    printf(" > rovibronic transition matrix elements\n");

    printf("\n");
    printf("(J',v')    ground level\n");
    printf("(J'',v'')  excited level\n");

    // outer two loops: over J' and J''
    for (int J1 = J_min; J1 <= J_max; J1++) {
        int index_J1 = J1 - J_min;

        for (int J2 = J_min; J2 <= J_max; J2++) {
            int index_J2 = J2 - J_min;

            printf("\n    J' J''  v' v''        e1,cm^-1        e2,cm^-1   delta e,cm^-1            prop        |prop|^2\n");

            // inner two loops: over v' and v''
            for (int v1 = v_min; v1 <= MIN(v_max, nroots1[index_J1]-1); v1++) {
                int index_v1 = v1 - v_min;

                for (int v2 = v_min; v2 <= MIN(v_max, nroots2[index_J2]-1); v2++) {
                    int index_v2 = v2 - v_min;

                    // bra
                    double E1 = energies1[index_J1][index_v1];
                    double E1_cm = (E1 - min_energy) * ATOMIC_TO_CM;
                    double *psi1 = wavefunctions1[index_J1][index_v1];

                    // ket
                    double E2 = energies2[index_J2][index_v2];
                    double E2_cm = (E2 - min_energy) * ATOMIC_TO_CM;
                    double *psi2 = wavefunctions2[index_J2][index_v2];

                    double prp12 = matrix_element_spline(n_grid, radial_grid, psi1, psi2, prop_fun);

                    printf("$ %4d%4d%4d%4d%16.4f%16.4f%16.4f%16.6f%16.6f\n", J1, J2, v1, v2, E1_cm, E2_cm, E2_cm-E1_cm, prp12, prp12*prp12);
                }
            }

        }
    }

    printf("\n");

    free(radial_grid);
}
