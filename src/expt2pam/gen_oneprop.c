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

/**
 * Backend to the OneProp code by Leonid Skripnikov.
 * Generates the basis.L file.
 *
 * 2021 Alexander Oleynichenko
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basis.h"
#include "ecp.h"
#include "molecule.h"


#define MAX_N_PRIMITIVES 200
#define MAX_N_CONTRACTED 200


void print_oneprop_table_of_contractions(FILE *out,
                                         int n_primitives, int n_contracted,
                                         double *exponents, double coeffs[][MAX_N_PRIMITIVES]);


void gen_oneprop(FILE *out, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib)
{
    char elem_symbol[8];
    static double exponents[MAX_N_PRIMITIVES];
    static double coeffs[MAX_N_CONTRACTED][MAX_N_PRIMITIVES];

    fprintf(out, "basis -> basis.L\n");

    FILE *bas_file = fopen("basis.L", "w");

    for (int z = 0; z < N_CHEM_ELEMENTS; z++) {
        basis_t *basis = bas_lib->basis_list[z];
        if (basis == NULL) {
            continue;
        }

        int L_max = basis_get_max_ang_mom(basis);
        get_element_symbol(z, elem_symbol);

        fprintf(bas_file, "%s   0\n", elem_symbol);

        for (int L = 0; L <= L_max; L++) {

            int n_primitives = 0;
            int n_contracted = 0;
            memset(exponents, 0, sizeof(double) * MAX_N_PRIMITIVES);
            memset(coeffs, 0, sizeof(double) * MAX_N_CONTRACTED * MAX_N_PRIMITIVES);

            for (int ifun = 0; ifun < basis->nfun; ifun++) {
                bfn_t *fun = basis->functions[ifun];
                if (fun->l != L) {
                    continue;
                }

                for (int i = 0; i < fun->nprim; i++) {
                    double exponent = fun->e[i];
                    double coef = fun->c[i];

                    // try to find this exponent
                    int found = 0;
                    for (int j = 0; j < n_primitives; j++) {
                        if (fabs(exponents[j] - exponent) < 1e-10) {
                            coeffs[n_contracted][j] = coef;
                            found = 1;
                            break;
                        }
                    }
                    // if exponent not found => add this primitive to the list
                    if (!found) {
                        exponents[n_primitives] = exponent;
                        coeffs[n_contracted][n_primitives] = coef;
                        n_primitives++;
                    }
                }

                n_contracted++;
            }

            // header for the block with L
            for (int i = 0; i < n_contracted; i++) {
                fprintf(bas_file, "%c", angular_momentum_to_char(L));
            }
            fprintf(bas_file, "%4d   1.0\n", n_primitives);

            // print table of exponents and coefficients
            print_oneprop_table_of_contractions(bas_file, n_primitives, n_contracted, exponents, coeffs);
        }

        fprintf(bas_file, "****\n");
    }

    fclose(bas_file);
}


void print_oneprop_table_of_contractions(FILE *out,
                                         int n_primitives, int n_contracted,
                                         double *exponents, double coeffs[][MAX_N_PRIMITIVES])
{
    for (int i = 0; i < n_primitives; i++) {
        fprintf(out, "%25.16f", exponents[i]);
        for (int j = 0; j < n_contracted; j++) {
            fprintf(out, "%25.16f", coeffs[j][i]);
        }
        fprintf(out, "\n");
    }
}
