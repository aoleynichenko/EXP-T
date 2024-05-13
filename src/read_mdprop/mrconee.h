//
// Created by Alexander Oleynichenko on 29.12.2023.
//

#ifndef MRCONEE_H
#define MRCONEE_H

#include <stdio.h>

typedef struct {
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
    int **mult_table;         // multiplication table for direct products in the Abelian subgroup
    int *occ_numbers;         // occupation numbers for each spinor
    int *spinor_irreps;       // irrep in Abelian subgroup
    double *spinor_energies;  // one-electron energies from SCF
    double _Complex *fock;    // Fock matrix
} mrconee_data_t;

mrconee_data_t *read_mrconee(char *path);

void free_mrconee_data(mrconee_data_t *data);

void print_mrconee_data(FILE *out, mrconee_data_t *data);

#endif //MRCONEE_H
