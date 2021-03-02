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

/*******************************************************************************
 * symmetry.h
 * ==========
 *
 * Tools for working with direct product tables.
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#ifndef CC_SYMMETRY_H_INCLUDED
#define CC_SYMMETRY_H_INCLUDED

#include "comdef.h"

enum {
    CC_GROUP_REAL = 1,
    CC_GROUP_COMPLEX = 2,
    CC_GROUP_QUATERNION = 4
};

// number of irreps
extern int nsym;

// is this group abelian or not
extern int is_abelian;

// direct product table
extern int ***dir_prod_table;
extern int dir_prod_table_abelian[CC_MAX_NUM_IRREPS][CC_MAX_NUM_IRREPS];

// representation names
extern char **rep_names;

// number of totally symmetry irrep
extern int irrep_a1;

// point group type (Quaternion, Real, Complex)
extern int point_group_nz;

int get_num_irreps();

int get_totally_symmetric_irrep();

int inverse_irrep_abelian(int irep);

void print_symmetry_info();

void symmetry_cleanup();

int get_rep_number(char *name);

char *get_irrep_name(int irep);

int mulrep2_abelian(int irep1, int irep2);

// checks if the direct product of n irreps from the rep_list contains the
// totally symmetry irrep
extern int (*dpd_contains_totsym_rank[CC_MAX_NUM_IRREPS])(int *gamma);

void set_point_group_name(char *name);

#endif /* CC_SYMMETRY_H_INCLUDED */
