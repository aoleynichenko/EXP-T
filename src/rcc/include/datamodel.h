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

/*
 * Any CC code operates with integrals, amplitudes and their intermediate
 * combinations. These entities can be represented as DIAGRAMS; each DIAGRAM
 * consists of BLOCKS (Direct Product Decomposition + tiling).
 *
 * Data structures implemented in the 'datamodel' module:
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

#ifndef CC_DATAMODEL_H_INCLUDED
#define CC_DATAMODEL_H_INCLUDED

#include <complex.h>
#include <stdio.h>

#include "comdef.h"

typedef enum {
    CC_DIAGRAM_IN_MEM,
    CC_DIAGRAM_ON_DISK,
    CC_DIAGRAM_DUMMY
} storage_type_t;

struct block {

    // unique (global) ID of this block
    // TODO: check for collisions with diagrams read from disk!
    int64_t id;

    // tensor rank: 2, 4, 6 ...
    int rank;

    // spinor blocks numbers for each dimension
    int spinor_blocks[CC_DIAGRAM_MAX_RANK];

    // sets of spinor indices for each dimension of the tensor
    // first element of each subarray -- total number of indices for this dim
    // example: for dimension 2:
    // indices[2] = {5, 0, 1, 2, 4, 5} // indices[2][0] == len({0,1,2,4,5}) == 5
    // len(indices) == rank
    int **indices;

    // total number of stored integrals
    // = product of all dimensions' sizes
    size_t size;

    // buffer (if in RAM, else zero)
    double complex *buf;

    // file (name) on disk (if needed, else zero)
    // TODO: must use block's unique id
    char *file_name;

    // flag: on disk or in RAM
    int storage_type;

    int is_unique;
    int sign;
    int n_equal_perms;
    int perm_from_unique[CC_DIAGRAM_MAX_RANK];
    int perm_to_unique[CC_DIAGRAM_MAX_RANK];
    int is_compressed;
};

typedef struct block block_t;

// constructor and destructor
block_t *block_new(int rank, int *spinor_blocks_nums, int *qparts, int *valence, int *t3space, int *order, int storage_type,
                   int only_unique);

void block_gen_indices(block_t *block, int *indices);

void block_delete(block_t *block);

void block_copy_data(block_t *dst, block_t *src);

void symblock_set(block_t *block, double complex val, int *idx);

double complex symblock_get(block_t *block, int *idx);

void symblock_print(block_t *);

void as_compound_index(int n /* rank */, int *dims, int lin_idx, int *idx);

int as_linear_index(int n /* rank */, int *dims, int *idx);

void block_clear(block_t *block);

void block_get_dims(block_t *block, int *dims);

void block_load(block_t *block);

void block_unload(block_t *block);

void block_store(block_t *block);

void symblock_write(int fd, block_t *block);

void block_write_formatted(FILE *txt_file, block_t *block);

void compress_triples_rank6(block_t *block);

void decompress_triples_rank6(block_t *block);

block_t *symblock_read(int fd);

struct diagram {

    // unique ID and name of the diagram
    int64_t dg_id;
    char name[CC_DIAGRAM_MAX_NAME];

    // the following metainfo coincides with the same info in each block
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
};

typedef struct diagram diagram_t;

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
diagram_t *diagram_reorder(diagram_t *dg, int *perm, int perm_unique);

// prints diagram's metainfo
void diagram_print(diagram_t *);

double complex diagram_get(diagram_t *dg, int *idx);

double complex diagram_get_2(diagram_t *dg, int i, int j);

double complex diagram_get_4(diagram_t *dg, int i, int j, int k, int l);

void diagram_set(diagram_t *dg, double complex val, int *idx);

void diagram_summary(diagram_t *diag);

void print_unique_info(diagram_t *diag);

// read & write operations
void diagram_write(diagram_t *dg, char *file_name);

diagram_t *diagram_read(char *file_name);

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

// diagrams are organized as a stack:
diagram_t *diagram_stack_push(diagram_t *dg);

diagram_t *diagram_stack_replace(char *name, diagram_t *dg);

diagram_t *diagram_stack_find(char *name);

void diagram_stack_erase(char *name);

int diagram_stack_find_index(char *name);

void diagram_stack_print();

typedef int dg_stack_pos_t;

dg_stack_pos_t get_stack_pos();

void restore_stack_pos(dg_stack_pos_t pos);

void restore_block(diagram_t *dg, block_t *b);

void restore_diagram(diagram_t *dg);

void block_unique(block_t *b, int *qparts, int *valence, int *order);

void destroy_block(block_t *b);

void transform(int n, int *idx, int *out, int *perm, int shift);

void print_n_restored();

#endif /* CC_DATAMODEL_H_INCLUDED */
