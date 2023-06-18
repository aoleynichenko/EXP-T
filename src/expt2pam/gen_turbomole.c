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

/**
 * Turbomole backend.
 *
 * 2021 Alexander Oleynichenko
 */

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basis.h"
#include "ecp.h"
#include "elements.h"
#include "molecule.h"

#define MAX_ANG_MOM    10
#define MAX_CONTRACTED 100

void gen_cfour_GENBAS(FILE *out, basis_lib_t *bas_lib);
void gen_cfour_ECPDATA(FILE *out, ecp_lib_t *ecp_lib);
void gen_cfour_ZMAT(FILE *out, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib);
void gen_turbomole_coord(FILE *out, molecule_t *mol);
void gen_turbomole_ecp(FILE *out, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib);
char angular_momentum_to_char(int l);
int basis_get_max_ang_mom(basis_t *bas);
void str_toupper(char *s);
void str_tolower(char *s);


void gen_turbomole(FILE *out, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib)
{
    FILE *f = fopen("coord", "w");
    gen_turbomole_coord(f, mol);
    fprintf(out, "geometry -> coord\n");
    fclose(f);

    gen_turbomole_ecp(out, bas_lib, ecp_lib);
}


void gen_turbomole_coord(FILE *out, molecule_t *mol)
{
    char elem_sym[8];

    fprintf(out, "$coord\n");
    for (int i = 0; i < mol->n_atoms; i++) {
        get_element_symbol(mol->charges[i], elem_sym);
        str_tolower(elem_sym);
        fprintf(out, "%22.14f%22.14f%22.14f\n", mol->x[i], mol->y[i], mol->z[i]);
    }
    fprintf(out, "$end\n");
}


/**
 * different basis+ecp files will be generated for different elements
 */
void gen_turbomole_ecp(FILE *out, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib)
{
    for (int ielement = 0; ielement < N_CHEM_ELEMENTS; ielement++) {

        char elem_sym[8];
        get_element_symbol(ielement, elem_sym);
        str_tolower(elem_sym);

        basis_t *basis = basis_lib_get(bas_lib, elem_sym);
        ecp_t *ecp = ecp_lib_get(ecp_lib, elem_sym);

        if (basis == NULL && ecp == NULL) {
            continue;
        }

        FILE *f = fopen(elem_sym, "w");

        /* 1. basis set */
        if (basis != NULL) {
            printf("basis -> %s\n", elem_sym);

            // header
            fprintf(f, "$basis\n");
            fprintf(f, "*\n");
            fprintf(f, "%s custom\n", elem_sym);
            fprintf(f, "#########################################\n");
            fprintf(f, "# this file was generated automatically #\n");
            fprintf(f, "#########################################\n");
            fprintf(f, "*\n");

            // contracted functions
            for (int ifun = 0; ifun < basis->nfun; ifun++) {
                bfn_t *func = basis->functions[ifun];
                int nprim = func->nprim;
                char am_sym = tolower(angular_momentum_to_char(func->l));
                fprintf(f, "%5d  %c\n", nprim, am_sym);
                for (int i = 0; i < nprim; i++) {
                    fprintf(f, "%18.10f%18.10f\n", func->e[i], func->c[i]);
                }
            }

            // '*' = end of basis
            fprintf(f, "*\n\n");
        }

        /* 2. effective core potential */
        if (ecp != NULL) {
            printf("ecp -> %s\n", elem_sym);

            int len_arep, len_esop;
            ecp_get_len(ecp, &len_arep, &len_esop);
            int lmax = len_arep-1;

            // header
            fprintf(f, "$ecp\n");
            fprintf(f, "*\n");
            fprintf(f, "%s custom\n", elem_sym);
            fprintf(f, "*\n");
            fprintf(f, "          ncore= %3d   lmax=  %2d   lsomax=  %2d\n", ecp->n_core_elec, lmax, len_esop);
            fprintf(f, "#        coefficient   r^n          exponent\n");

            // scalar part
            char lmax_sym = tolower(angular_momentum_to_char(lmax));
            for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
                ecp_expansion_t *arep_f = ecp->arep[i];
                if (arep_f == NULL && i < len_arep) { // no expansion, but it is needed
                    fprintf(f, "%c-%c\n", tolower(angular_momentum_to_char(i - 1)), lmax_sym);
                    fprintf(f, "%25.16E%4d%25.16E\n", 0.0, 0, 1.0);
                    continue;
                }
                else if (arep_f == NULL) {
                    continue;
                }
                if (i == 0) {
                    fprintf(f, "%c\n", lmax_sym);
                }
                else {
                    fprintf(f, "%c-%c\n", tolower(angular_momentum_to_char(i - 1)), lmax_sym);
                }
                for (int j = 0; j < arep_f->nprim; j++) {
                    fprintf(f, "%25.16E%4d%25.16E\n", arep_f->c[j], arep_f->powers[j], arep_f->e[j]);
                }
            }
            fprintf(f, "*\n");

            // spin-orbit part
            for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
                ecp_expansion_t *so_f = ecp->esop[i];
                if (so_f == NULL && i < len_esop && i >= 2) { // no expansion, but it is needed
                    fprintf(f, "%c-SpinOrbit\n", tolower(angular_momentum_to_char(i - 1)));
                    fprintf(f, "%25.16E%4d%25.16E\n", 0, 1.0, 0.0);
                    continue;
                }
                else if (so_f == NULL) {
                    continue;
                }
                else {
                    fprintf(f, "%c-SpinOrbit\n", tolower(angular_momentum_to_char(i - 1)));
                }
                for (int j = 0; j < so_f->nprim; j++) {
                    fprintf(f, "%25.16E%4d%25.16E\n", so_f->c[j], so_f->powers[j], so_f->e[j]);
                }
            }
            fprintf(f, "*\n");
        }

        fprintf(f, "$end\n");

        fclose(f);
    }
}

