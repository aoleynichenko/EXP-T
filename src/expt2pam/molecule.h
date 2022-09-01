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

#ifndef MOLECULE_H_INCLUDED
#define MOLECULE_H_INCLUDED

#include "elements.h"

#define MAX_N_ATOMS 100

enum {
    SYMMETRY_AUTO,
    SYMMETRY_C1,
    SYMMETRY_Ci,
    SYMMETRY_Cs,
    SYMMETRY_C2,
    SYMMETRY_C2v,
    SYMMETRY_C2h,
    SYMMETRY_D2,
    SYMMETRY_D2h,
    SYMMETRY_Cinfv,
    SYMMETRY_Dinfh
};

typedef struct {
    int group;   // group label (one of SYMMETRY_* constants)
    int xyz[3];  // main group element (axis for axial groups and C2h, plane for Cs)
} symgrp_t;

enum {
    CAGE_ICOSAHEDRON,
    CAGE_TETRAHEDRON,
    CAGE_OCTAHEDRON,
    CAGE_DODECAHEDRON
};

typedef struct {
    int n_points;
    double total_charge;
    double *points;
} charged_cage_t;

typedef struct {
    int n_atoms;
    int charges[MAX_N_ATOMS];
    double x[MAX_N_ATOMS];
    double y[MAX_N_ATOMS];
    double z[MAX_N_ATOMS];
    symgrp_t sym_group;
    charged_cage_t cage;
} molecule_t;

molecule_t *molecule_new();
void molecule_delete(molecule_t *mol);
void molecule_add_atom(molecule_t *mol, int nuc_charge, double x, double y, double z);
void molecule_print(molecule_t *mol);
int molecule_n_atom_types(molecule_t *mol);
int molecule_n_atoms_of(molecule_t *mol, int element);

#endif /* MOLECULE_H_INCLUDED */
