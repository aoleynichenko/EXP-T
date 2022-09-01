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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ecp.h"


/**
 * constructor
 */
ecp_expansion_t *ecp_expansion_new(int nprim, int *powers, double *e, double *c)
{
    ecp_expansion_t *f = (ecp_expansion_t *) malloc(sizeof(ecp_expansion_t));

    f->powers = (int *) malloc(sizeof(int) * nprim);
    f->e = (double *) malloc(sizeof(double) * nprim);
    f->c = (double *) malloc(sizeof(double) * nprim);

    f->nprim = nprim;
    memcpy(f->powers, powers, sizeof(int) * nprim);
    memcpy(f->e, e, sizeof(double) * nprim);
    memcpy(f->c, c, sizeof(double) * nprim);

    return f;
}


/**
 * destructor
 */
void ecp_expansion_delete(ecp_expansion_t *f)
{
    free(f->powers);
    free(f->e);
    free(f->c);
    free(f);
}


ecp_t *ecp_new(int element)
{
    ecp_t *ecp = (ecp_t *) malloc(sizeof(ecp_t));

    ecp->element = element;
    ecp->n_core_elec = 0;
    memset(ecp->arep, 0, sizeof(ecp->arep));
    memset(ecp->esop, 0, sizeof(ecp->esop));

    return ecp;
}


void ecp_delete(ecp_t *ecp)
{
    int nfun = sizeof(ecp->arep) / sizeof(ecp_expansion_t *);

    for (int i = 0; i < nfun; i++) {
        if (ecp->arep[i]) {
            ecp_expansion_delete(ecp->arep[i]);
        }
        if (ecp->esop[i]) {
            ecp_expansion_delete(ecp->esop[i]);
        }
    }

    free(ecp);
}


void ecp_add_function(ecp_t *ecp, int type, int ang_mom, int nprim, int *powers, double *e, double *c)
{
    assert(type == ECP_AVERAGED || type == ECP_SPIN_ORBIT);

    ecp_expansion_t *f = ecp_expansion_new(nprim, powers, e, c);
    ecp_expansion_t **fun_list = (type == ECP_AVERAGED) ? ecp->arep : ecp->esop;

    assert(ang_mom <= ECP_MAX_ANG_MOM);

    int pos;
    if (ang_mom == ECP_UL) {
        pos = 0;
    }
    else {
        pos = ang_mom + 1;
    }

    if (fun_list[pos] != NULL) {
        ecp_expansion_delete(fun_list[pos]);
    }
    fun_list[pos] = f;
}


void ecp_get_len(ecp_t *ecp, int *len_arep, int *len_esop)
{
    *len_arep = 0;
    *len_esop = 0;

    for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
        ecp_expansion_t *f = ecp->arep[i];
        if (f != NULL) {
            *len_arep = i + 1;
        }
    }
    for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
        ecp_expansion_t *f = ecp->esop[i];
        if (f != NULL) {
            *len_esop = i - 1;  // since begins from L = P
        }
    }
}


void ecp_print(ecp_t *ecp)
{
    printf("\nECP for E%d\n", ecp->element);
    printf(" %d core electrons\n", ecp->n_core_elec);
    printf(" arep:\n");
    for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
        ecp_expansion_t *f = ecp->arep[i];
        if (f == NULL) {
            continue;
        }
        if (i == 0) {
            printf("L = Ul\n");
        }
        else {
            printf("L = %d\n", i - 1);
        }
        for (int j = 0; j < f->nprim; j++) {
            printf("%4d%18.8e%18.8e\n", f->powers[j], f->e[j], f->c[j]);
        }
    }
    printf(" esop:\n");
    for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
        ecp_expansion_t *f = ecp->esop[i];
        if (f == NULL) {
            continue;
        }
        if (i == 0) {
            printf("L = Ul\n");
        }
        else {
            printf("L = %d\n", i - 1);
        }
        for (int j = 0; j < f->nprim; j++) {
            printf("%4d%18.8e%18.8e\n", f->powers[j], f->e[j], f->c[j]);
        }
    }
    printf("\n");
}


ecp_lib_t *ecp_lib_new()
{
    ecp_lib_t *ecp_lib = (ecp_lib_t *) malloc(sizeof(ecp_lib_t));

    memset(ecp_lib->ecp_list, 0, sizeof(ecp_lib->ecp_list));

    return ecp_lib;
}


void ecp_lib_delete(ecp_lib_t *ecp_lib)
{
    int sz = sizeof(ecp_lib->ecp_list) / sizeof(ecp_t *);

    for (int i = 0; i < sz; i++) {
        ecp_t *ecp = ecp_lib->ecp_list[i];
        if (ecp != NULL) {
            ecp_delete(ecp);
        }
    }

    free(ecp_lib);
}


ecp_t *ecp_lib_set(ecp_lib_t *ecp_lib, ecp_t *ecp)
{
    int z = ecp->element;
    if (z < 0 || z > N_CHEM_ELEMENTS) {
        return NULL;
    }

    ecp_lib->ecp_list[z] = ecp;
    return ecp;
}


ecp_t *ecp_lib_get(ecp_lib_t *ecp_lib, char *elem_symbol)
{
    int z = get_element_nuc_charge(elem_symbol);
    if (z == -1) {
        return NULL;
    }

    return ecp_lib->ecp_list[z];
}


void ecp_lib_print(ecp_lib_t *ecp_lib)
{
    for (int i = 0; i < N_CHEM_ELEMENTS; i++) {
        if (ecp_lib->ecp_list[i] != NULL) {
            char buf[8];
            get_element_symbol(i, buf);
            printf(" ===================== %s ======================\n", buf);
            ecp_print(ecp_lib->ecp_list[i]);
        }
    }
}
