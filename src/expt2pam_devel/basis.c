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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basis.h"

#define ZERO_THRESH 1e-12


bfn_t *bfn_new(int ang_mom, int nprim, double *e, double *c)
{
    bfn_t *bf;
    int n_nonzero_coef = 0;

    for (int i = 0; i < nprim; i++) {
        if (fabs(c[i]) > ZERO_THRESH) {
            n_nonzero_coef++;
        }
    }

    if (n_nonzero_coef == 0) {
        return NULL;
    }

    bf = (bfn_t *) malloc(sizeof(bfn_t));
    bf->e = (double *) malloc(sizeof(double) * n_nonzero_coef);
    bf->c = (double *) malloc(sizeof(double) * n_nonzero_coef);

    bf->l = ang_mom;
    bf->nprim = 0;
    for (int i = 0; i < nprim; i++) {
        if (fabs(c[i]) > ZERO_THRESH) {
            bf->e[bf->nprim] = e[i];
            bf->c[bf->nprim] = c[i];
            bf->nprim++;
        }
    }

    return bf;
}


void bfn_delete(bfn_t *bf)
{
    free(bf->e);
    free(bf->c);
    free(bf);
}


int bfn_same_exponents(bfn_t *b1, bfn_t *b2)
{
    if (b1->nprim != b2->nprim) {
        return 0;
    }

    for (int i = 0; i < b1->nprim; i++) {
        if (fabs(b1->e[i] - b2->e[i]) > ZERO_THRESH) {
            return 0;
        }
    }

    return 1;
}


basis_t *basis_new(int element)
{
    basis_t *bas;

    bas = (basis_t *) malloc(sizeof(basis_t));

    bas->element = element;
    bas->nfun = 0;
    bas->functions = NULL;

    return bas;
}


void basis_delete(basis_t *bas)
{
    for (int i = 0; i < bas->nfun; i++) {
        bfn_t *bf = bas->functions[i];
        bfn_delete(bf);
    }
    free(bas->functions);
    free(bas);
}


void basis_add_function(basis_t *bas, int ang_mom, int nprim, double *e, double *c)
{
    bfn_t *bf = bfn_new(ang_mom, nprim, e, c);
    if (bf == NULL) {
        return;
    }

    if (bas->nfun == 0) {
        bas->nfun = 1;
        bas->functions = (bfn_t **) malloc(sizeof(bfn_t *) * 1);
    }
    else {
        bas->nfun++;
        bas->functions = (bfn_t **) realloc(bas->functions, sizeof(bfn_t *) * bas->nfun);
    }

    bas->functions[bas->nfun - 1] = bf;
}


void basis_print(basis_t *bas)
{
    printf("Basis set for E%d\n", bas->element);
    printf("Number of contracted basis functions: %d\n", bas->nfun);
    for (int i = 0; i < bas->nfun; i++) {
        bfn_t *bf = bas->functions[i];
        printf("[%d] L=%d nprim=%d\n", i, bf->l, bf->nprim);
        for (int j = 0; j < bf->nprim; j++) {
            printf("\t%18.8e%18.8e\n", bf->e[j], bf->c[j]);
        }
    }
}


basis_lib_t *basis_lib_new()
{
    basis_lib_t *bas_lib = (basis_lib_t *) malloc(sizeof(basis_lib_t));

    memset(bas_lib->basis_list, 0, sizeof(bas_lib->basis_list));

    return bas_lib;
}


void basis_lib_delete(basis_lib_t *bas_lib)
{
    int sz = sizeof(bas_lib->basis_list) / sizeof(basis_t *);

    for (int i = 0; i < sz; i++) {
        basis_t *bas = bas_lib->basis_list[i];
        if (bas != NULL) {
            basis_delete(bas);
        }
    }

    free(bas_lib);
}


basis_t *basis_lib_set(basis_lib_t *bas_lib, basis_t *bas)
{
    int z = bas->element;
    if (z < 0 || z > N_CHEM_ELEMENTS) {
        return NULL;
    }

    bas_lib->basis_list[z] = bas;
    return bas;
}


basis_t *basis_lib_get(basis_lib_t *bas_lib, char *elem_symbol)
{
    int z = get_element_nuc_charge(elem_symbol);
    if (z == -1) {
        return NULL;
    }

    return bas_lib->basis_list[z];
}


void basis_lib_print(basis_lib_t *bas_lib)
{
    for (int i = 0; i < N_CHEM_ELEMENTS; i++) {
        if (bas_lib->basis_list[i] != NULL) {
            char buf[8];
            get_element_symbol(i, buf);
            printf(" ===================== %s ======================\n", buf);
            basis_print(bas_lib->basis_list[i]);
        }
    }
}
