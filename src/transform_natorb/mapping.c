/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2024 The EXP-T developers.
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

#include "mapping.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "errquit.h"


/*
 * constructs mapping:
 * active spinor index -> general spinor index
 *
 * spinors are matched using the energy criterion:
 * eps( active[i] ) = eps( general[j] )  ==>  we have mapping: i -> j
 *
 * for some active spinors not molecular spinor with the same energy,
 * but its Kramers partner is to be used.
 *
 * input:
 *   n_act_spinors       number of active spinors
 *   n_all_spinors       total number of molecular spinors (only one Kramers partner)
 *   active_eigenvalues  list of active spinor energies
 *   all_eigenvalues     list of energies of all molecular spinors
 *
 * output:
 *   act_to_all          mapping (will be allocated)
 *   kramers_conj        mapping onto Kramers partner (1) or not (0) (will be allocated)
 */
void construct_spinor_mapping(int n_act_spinors, double *active_eigenvalues,
                              int n_all_spinors, double *all_eigenvalues,
                              int **act_to_all, int **kramers_conj)
{
    const double COMPARE_THRESH = 1e-8;
    int is_kramers_conj = 0;
    int n_matches = 0;

    printf(" > construction of the mapping: active spinors --> all spinors\n");

    *act_to_all = (int *) calloc(n_act_spinors, sizeof(int));
    *kramers_conj = (int *) calloc(n_act_spinors, sizeof(int));

    int ispi = 0;
    for (int iact = 0; iact < n_act_spinors; iact++) {
        double eps_act = active_eigenvalues[iact];

        if (iact == n_act_spinors / 2) {
            ispi = 0;  // now go to the Kramers partners
            is_kramers_conj = 1;
        }

        for (; ispi < n_all_spinors; ispi++) {
            if (fabs(all_eigenvalues[ispi] - eps_act) < COMPARE_THRESH) {
                (*act_to_all)[iact] = ispi;
                (*kramers_conj)[iact] = is_kramers_conj;
                ispi++;
                n_matches++;
                break;
            }
        }
    }

    if (n_matches == 0) {
        errquit("no matches of spinor energies");
    }

    // print mapping
    for (int i = 0; i < n_act_spinors; i++) {
        printf("%6d   --> %4d%2s     eps = %20.12f\n", i + 1, (*act_to_all)[i] + 1, (*kramers_conj)[i] ? "*" : "", active_eigenvalues[i]);
    }
    printf("\n");
}

