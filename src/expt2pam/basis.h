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

#ifndef BASIS_H_INCLUDED
#define BASIS_H_INCLUDED

#include "elements.h"

/**
 * contracted gaussian basis function
 */
typedef struct {
    int l;       /* angular momentum */
    int nprim;   /* length of expansion -- number of primitives */
    double *e;   /* exponents (alpha) */
    double *c;   /* contraction coefficients */
} bfn_t;

bfn_t *bfn_new(int ang_mom, int nprim, double *e, double *c);
void bfn_delete(bfn_t *bf);
int bfn_same_exponents(bfn_t *b1, bfn_t *b2);

typedef struct {
    int element;
    int nfun;
    bfn_t **functions;
} basis_t;

basis_t *basis_new(int element);
void basis_delete(basis_t *bas);
void basis_add_function(basis_t *bas, int ang_mom, int nprim, double *e, double *c);
void basis_print(basis_t *bas);

typedef struct {
    basis_t *basis_list[N_CHEM_ELEMENTS];
} basis_lib_t;

basis_lib_t *basis_lib_new();
void basis_lib_delete(basis_lib_t *);
basis_t *basis_lib_set(basis_lib_t *bas_lib, basis_t *bas);
basis_t *basis_lib_get(basis_lib_t *bas_lib, char *elem_symbol);
int basis_get_max_ang_mom(basis_t *bas);
void basis_lib_print(basis_lib_t *bas_lib);
char angular_momentum_to_char(int l);

#endif /* BASIS_H_INCLUDED */
