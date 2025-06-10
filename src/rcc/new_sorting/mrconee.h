/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2025 The EXP-T developers.
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


#ifndef CC_MRCONEE_H_INCLUDED
#define CC_MRCONEE_H_INCLUDED

#include <stdio.h>

typedef struct {
    int dirac_int_size;       // size of integer numbers used by DIRAC: 4-byte or 8-byte
    int num_spinors;          // total number of spinors
    double nuc_rep_energy;    // core energy (inactive energy + nuclear repulsion)
    int group_arith;          // group type (1 real, 2 complex, 4 quaternion)
    int is_spinfree;          // spinfree formalism
    double scf_energy;        // total SCF energy (from DIRAC)
    int invsym;               // inversion symmetry (yes : 2; no : 1)
    int num_irreps;           // number of fermion irreps in the Abelian subgroup
    char *point_group;        // point group symbol
    int totally_sym_irrep;    // number of a totally symmetric irrep
    char **irrep_names;       // names of these irreps
    int *mult_table;          // multiplication table for direct products in the Abelian subgroup
    int *occ_numbers;         // occupation numbers for each spinor
    int *spinor_irreps;       // irrep in Abelian subgroup
    double *spinor_energies;  // one-electron energies from SCF
    double _Complex *fock;    // Fock matrix
} mrconee_data_t;

mrconee_data_t *read_mrconee(char *path);

void free_mrconee_data(mrconee_data_t *data);

void print_mrconee_data(FILE *out, mrconee_data_t *data);

#endif // CC_MRCONEE_H_INCLUDED
