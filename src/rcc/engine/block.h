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

#ifndef CC_BLOCK_H_INCLUDED
#define CC_BLOCK_H_INCLUDED

#include <complex.h>
#include <stdio.h>

#include "comdef.h"

typedef enum {
    CC_DIAGRAM_IN_MEM,
    CC_DIAGRAM_ON_DISK,
    CC_DIAGRAM_DUMMY
} storage_type_t;

typedef struct block_t {

    // unique (global) ID of this block
    // TODO: check for collisions with diagrams read from disk!
    int64_t id;

    // tensor rank: 2, 4, 6 ...
    int rank;

    // spinor blocks numbers for each dimension
    int spinor_blocks[CC_DIAGRAM_MAX_RANK];

    // shape of the tensor containing matrix elements / amplitudes
    int *shape;

    // sets of spinor indices for each dimension of the tensor
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
} block_t;

// constructor and destructor
block_t *block_new(int rank, int *spinor_blocks_nums, int *qparts, int *valence, int *t3space, int *order, int storage_type,
                   int only_unique);

void block_gen_indices(block_t *block, int *indices);

void block_delete(block_t *block);

void block_copy_data(block_t *dst, block_t *src);

void block_set_element(block_t *block, double complex val, int *idx);

double complex block_get_element(block_t *block, int *idx);

void block_debug_print(block_t *block);

void tensor_index_to_compound(int rank /* rank */, const int *shape, size_t linear_idx, int *idx);

size_t tensor_index_to_linear(int rank /* rank */, const int *shape, const int *idx);

void block_clear(block_t *block);

void block_load(block_t *block);

void block_unload(block_t *block);

void block_store(block_t *block);

void block_write_binary(int fd, block_t *block);

void block_write_file_formatted(FILE *txt_file, block_t *block);

void compress_triples_rank6(block_t *block);

void decompress_triples_rank6(block_t *block);

block_t *block_read_binary(int fd);

void transform(int n, int *idx, int *out, int *perm, int shift);

#endif // CC_BLOCK_H_INCLUDED
