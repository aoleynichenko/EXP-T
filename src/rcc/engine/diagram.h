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

/*
 * Any CC code operates with integrals, amplitudes and their intermediate
 * combinations. These entities can be represented as DIAGRAMS; each DIAGRAM
 * consists of BLOCKS (Direct Product Decomposition + tiling).
 *
 * Data structures implemented in the 'diagram' module:
 *  diagram_t
 *  block_t
 *
 * To learn more about DPD, please see the following papers:
 *                             [basic papers]
 * => P. Čársky, L.J. Schaad, B.A. Hess, M. Urban, J. Noga.
 *    J. Chem. Phys., V. 87, P. 411 (1987)
 * => J.F. Stanton, J. Gauss, J.D. Watts, R.J. Bartlett.
 *    J. Chem. Phys., V. 94, P. 4334 (1991)
 *                           [relativistic DPD]
 * => A. Shee, L. Visscher, T. Saue. J. Chem. Phys., V. 145, P. 184107 (2016)
 * => L. Visscher, T. Lee, K. Dyall. J. Chem. Phys., V. 105, P. 8769 (1996)
 */

#ifndef CC_DIAGRAM_H_INCLUDED
#define CC_DIAGRAM_H_INCLUDED

#include "block.h"

typedef struct diagram {

    // name of the diagram
    char name[CC_DIAGRAM_MAX_NAME];

    // tensor rank: 2, 4, 6 ...
    int rank;

    // symmetry of the operator represented by the diagram
    int symmetry;

    // list of quasiparticles (in natural form)
    // for each spinor
    // of size [rank] (and all the subsequent arrays too)
    int qparts[CC_DIAGRAM_MAX_RANK]; // iii

    // list of valence/inactive (in natural form)
    int valence[CC_DIAGRAM_MAX_RANK]; // iiv
    int t3space[CC_DIAGRAM_MAX_RANK];

    // "reordering" from the "natural" form
    int order[CC_DIAGRAM_MAX_RANK]; // iiu

    // inverted index of blocks
    size_t *inv_index;

    // links to symmetry blocks
    size_t n_blocks;
    block_t **blocks;

    int only_unique;
} diagram_t;


// creates new object of diagram -- "constructor", but don't binds it to the
// singly-linked list "dg_stack"
diagram_t *diagram_new(char *name, char *qparts, char *valence, char *t3space, char *order, int perm_unique, int irrep);

// deallocates all memory associated with this object
void diagram_delete(diagram_t *dg);

// "copy constructor"
diagram_t *diagram_copy(diagram_t *dg);

storage_type_t diagram_get_storage_type(diagram_t *dg);

void diagram_get_valence(diagram_t *diag, char *valence);

void diagram_get_t3space(diagram_t *diag, char *t3space);

void diagram_get_quasiparticles(diagram_t *diag, char *qparts);

void diagram_get_order(diagram_t *diag, char *order);

void diagram_get_memory_used(diagram_t *dg, size_t *ram_used, size_t *disk_used);


// reorder diagrams
diagram_t *diagram_reorder(diagram_t *diag, int *perm);

// prints diagram's metainfo
void diagram_debug_print(diagram_t *dg);

double complex diagram_get(diagram_t *dg, int *idx);

double complex diagram_get_2(diagram_t *dg, int i, int j);

double complex diagram_get_4(diagram_t *dg, int i, int j, int k, int l);

void diagram_set(diagram_t *dg, double complex val, int *idx);

void diagram_set_4(diagram_t *dg, double complex val, int i, int j, int k, int l);

void diagram_summary(diagram_t *diag);

void print_unique_info(diagram_t *diag);

// read & write operations
void diagram_write_binary(diagram_t *dg, char *file_name);

diagram_t *diagram_read_binary(char *file_name);

void diagram_write_formatted(diagram_t *dg, char *file_name);

// simple arithmetics
diagram_t *diagram_clear(diagram_t *dg);

diagram_t *diagram_add(diagram_t *dg_target,
                       double fact1, diagram_t *dg1,
                       double fact2, diagram_t *dg2);

double diagram_max(diagram_t *dg, int *indices);

double diagram_diffmax(diagram_t *dg1, diagram_t *dg2, int *indices);

// divide matrix elements by energy denominators
diagram_t *diagram_diveps(diagram_t *dg);

// diagram contraction
diagram_t *diagram_mult(diagram_t *dg1, diagram_t *dg2, int ncontr, int perm_unique);

// memory management: on disk or in RAM
void diagram_set_storage_type(diagram_t *dg, int storage_type);

block_t *diagram_get_block(diagram_t *dg, int *spinor_blocks_nums);//, size_t *block_index);

void set_order(char *dg_name, char *new_order);

void restore_block(diagram_t *dg, block_t *b);

void restore_diagram(diagram_t *dg);

void destroy_block(block_t *b);

void print_n_restored();

#endif // CC_DIAGRAM_H_INCLUDED
