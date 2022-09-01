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

/**
 * CFOUR backend.
 * basis set will be uncontracted.
 *
 * 2020 Alexander Oleynichenko
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
char angular_momentum_to_char(int l);
int basis_get_max_ang_mom(basis_t *bas);
void str_toupper(char *s);
void gen_cfour_print_identity_matrix(FILE *out, int n);


void gen_cfour(FILE *out, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib)
{
    FILE *f = fopen("GENBAS", "w");
    gen_cfour_GENBAS(f, bas_lib);
    fprintf(out, "basis set -> GENBAS\n");
    fclose(f);

    f = fopen("ECPDATA", "w");
    gen_cfour_ECPDATA(f, ecp_lib);
    fprintf(out, "ecp -> ECPDATA\n");
    fclose(f);

    f = fopen("ZMAT", "w");
    gen_cfour_ZMAT(f, mol, bas_lib, ecp_lib);
    fprintf(out, "geometry & template for settings -> ZMAT\n");
    fclose(f);
}


void gen_cfour_ZMAT(FILE *out, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib)
{
    char elem_sym[8];

    fprintf(out, "ZMAT file template\n");   // comment line

    for (int i = 0; i < mol->n_atoms; i++) {
        get_element_symbol(mol->charges[i], elem_sym);
        str_toupper(elem_sym);
        fprintf(out, "%2s%8.5f%8.5f%8.5f\n", elem_sym, mol->x[i], mol->y[i], mol->z[i]);
    }
    fprintf(out, "\n");

    // calculation setup
    fprintf(out, "*CFOUR(CALC=CCSD\n");
    fprintf(out, "BASIS=SPECIAL\n");
    fprintf(out, "ECP=ON\n");
    fprintf(out, "SPHERICAL=ON\n");
    fprintf(out, "UNITS=BOHR\n");
    fprintf(out, "SCF_MAXCYC=500\n");
    fprintf(out, "LINDEP_TOL=7\n");
    fprintf(out, "SCF_EXTRAPOLATION=ON\n");
    fprintf(out, "REF=UHF\n");
    fprintf(out, "SCF_DAMPING=100\n");
    fprintf(out, "MULT=<?>\n");
    fprintf(out, "CHARGE=<?>\n");
    fprintf(out, "ABCDTYPE=AOBASIS\n");
    fprintf(out, "CC_PROGRAM=VCC\n");
    fprintf(out, "CC_MAXCYC=1000\n");
    fprintf(out, "LINEQ_MAXCY=1000\n");
    fprintf(out, "SUBGRPAXIS=3\n");
    fprintf(out, "PROPS=1\n");
    fprintf(out, "MEM_UNIT=GB\n");
    fprintf(out, "MEMORY_SIZE=8)\n");
    fprintf(out, "\n");

    // basis sets
    for (int ielem = 0; ielem < N_CHEM_ELEMENTS; ielem++) {
        basis_t *bas = bas_lib->basis_list[ielem];
        if (bas != NULL) {
            get_element_symbol(bas->element, elem_sym);
            str_toupper(elem_sym);
            fprintf(out, "%s:UNCONT\n", elem_sym);
        }
    }
    fprintf(out, "\n");

    // ECPs
    for (int ielem = 0; ielem < N_CHEM_ELEMENTS; ielem++) {
        ecp_t *ecp = ecp_lib->ecp_list[ielem];
        if (ecp != NULL) {
            get_element_symbol(ecp->element, elem_sym);
            str_toupper(elem_sym);
            fprintf(out, "%s:ECP\n", elem_sym);
        }
    }
    fprintf(out, "\n");
}


void gen_cfour_GENBAS(FILE *out, basis_lib_t *bas_lib)
{
    const int MAX_NUM_BLOCKS = 32;
    int n_blocks = 0;
    int block_ang_momentum[MAX_NUM_BLOCKS];
    int block_num_functions[MAX_NUM_BLOCKS];

    for (int ielem = 0; ielem < N_CHEM_ELEMENTS; ielem++) {
        basis_t *bas = bas_lib->basis_list[ielem];
        if (bas == NULL) {
            continue;
        }

        /*
         * label of the chemical element and name of the basis set
         */
        char elem_sym[8];
        get_element_symbol(bas->element, elem_sym);
        str_toupper(elem_sym);
        fprintf(out, "\n%s:UNCONT\n", elem_sym);

        /*
         * determine maximum angular momentum of basis function
         */
        int max_L = basis_get_max_ang_mom(bas);

        /*
         * no more than 10 functions can be placed into one block.
         * if there are > 10 primitive Gaussians, several blocks are created.
         */
        for (int L = 0; L <= max_L; L++) {
            int n_fun_L = 0;
            for (int i = 0; i < bas->nfun; i++) {
                if (bas->functions[i]->l == L) {
                    n_fun_L++;
                }
            }

            while (n_fun_L > 0) {
                int MAX_FUN_PER_BLOCK = 10;

                // maximum number of functions in the block depends on L
                if (L == 6) {
                    MAX_FUN_PER_BLOCK = 4;
                }

                int n_fun_block = 0;
                if (n_fun_L > MAX_FUN_PER_BLOCK) {
                    n_fun_block = MAX_FUN_PER_BLOCK;
                }
                else {
                    n_fun_block = n_fun_L;
                }

                block_num_functions[n_blocks] = n_fun_block;
                block_ang_momentum[n_blocks] = L;
                n_blocks++;

                n_fun_L -= n_fun_block;
            }
        }

        /*
         * header:
         * comment line
         */
        for (int i = 0; i < n_blocks; i++) {
            int L = block_ang_momentum[i];
            fprintf(out, "    %c", tolower(angular_momentum_to_char(L)));
        }
        fprintf(out, "\n\n");

        /*
         * header:
         * number of basis functions in each block
         */
        fprintf(out, "%3d\n", n_blocks);
        for (int i = 0; i < n_blocks; i++) {
            int L = block_ang_momentum[i];
            fprintf(out, "    %1d", L);
        }
        fprintf(out, "\n");

        for (int i = 0; i < n_blocks; i++) {
            fprintf(out, "%5d", block_num_functions[i]);
        }
        fprintf(out, "\n");
        for (int i = 0; i < n_blocks; i++) {
            fprintf(out, "%5d", block_num_functions[i]);
        }
        fprintf(out, "\n");

        /*
         * for each block:
         * print exponents and identity matrix of coeff-s
         */
        for (int L = 0; L <= max_L; L++) {
            int offs = 0;

            for (int iblock = 0; iblock < n_blocks; iblock++) {

                if (L != block_ang_momentum[iblock]) {
                    continue;
                }

                fprintf(out, "\n");

                /*
                 * print exponents
                 */
                int count = 0;
                int i;
                for (i = offs; i < bas->nfun; i++) {
                    if (bas->functions[i]->l == L) {

                        if (count < block_num_functions[iblock]) {
                            double e = bas->functions[i]->e[0];
                            fprintf(out, "%14.8f", e);
                            if ((count+1) % 5 == 0) {
                                fprintf(out, "\n");
                            }
                        }
                        else {
                            break;
                        }

                        count++;
                    }
                }
                offs = i;

                if ((count) % 5 != 0) {
                    fprintf(out, "\n");
                }
                fprintf(out, "\n");

                /*
                 * coefficients
                 */
                gen_cfour_print_identity_matrix(out, block_num_functions[iblock]);
            }
        }
    }
}


void gen_cfour_print_identity_matrix(FILE *out, int n)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                fprintf(out, " 1.");
            }
            else {
                fprintf(out, " 0.");
            }
        }
        fprintf(out, "\n");
    }
}


void gen_cfour_ECPDATA(FILE *out, ecp_lib_t *ecp_lib)
{
    char elem_sym[8];

    for (int ielem = 0; ielem < N_CHEM_ELEMENTS; ielem++) {
        ecp_t *ecp = ecp_lib->ecp_list[ielem];
        if (ecp == NULL) {
            continue;
        }

        get_element_symbol(ecp->element, elem_sym);
        str_toupper(elem_sym);
        fprintf(out, "*\n");
        fprintf(out, "%s:ECP\n", elem_sym);
        fprintf(out, "#\n");
        fprintf(out, "*\n");

        int len_arep, len_esop;
        ecp_get_len(ecp, &len_arep, &len_esop);
        fprintf(out, "  NCORE =   %3d           LMAX =     %1d\n", ecp->n_core_elec, len_arep - 1);
        fprintf(out, "#\n");

        for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
            ecp_expansion_t *f = ecp->arep[i];
            if (f == NULL && i < len_arep) { // no expansion, but it is needed
                fprintf(out, "%24.12e%4d%24.12e\n", 0.0, 0, 1.0);
                continue;
            }
            else if (f == NULL) {
                continue;
            }
            if (i == 0) {
                fprintf(out, "%c\n", tolower(angular_momentum_to_char(len_arep - 1)));
            }
            else {
                fprintf(out, "%c-%c\n", tolower(angular_momentum_to_char(i - 1)),
                        tolower(angular_momentum_to_char(len_arep - 1)));
            }
            for (int j = 0; j < f->nprim; j++) {
                fprintf(out, "%24.12e%4d%24.12e\n", f->c[j], f->powers[j], f->e[j]);
            }
        }
        fprintf(out, "*\n");
    }
}

