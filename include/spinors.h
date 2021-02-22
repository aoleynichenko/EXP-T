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
 * spinors.h
 * =========
 *
 * All information about one-particle functions (spin-orbitals or spinors)
 * (and symmetry).
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#ifndef CC_SPINORS_H_INCLUDED
#define CC_SPINORS_H_INCLUDED

#include "comdef.h"

// number of one-particle functions (e.g. spinors)
extern int nspinors;

// attributes of each spinors:
//  - sequential number
//  - irrep no
//  - hole / particle (occupied in reference or not) (1 - hole, 0 - particle)
//  - active or not (1 or 0)
//  - space1: if 1, triples amplitude with this spinor can be non-zero
//    (otherwise it will be set to zero, "restriction of triples")
//  - its energy
typedef struct {
    int seqno;
    int repno;
    int blockno;
    int occ;
    int active;
    int space1;
    double eps;
} spinor_attr_t;

// array of spinor attributes
extern spinor_attr_t *spinor_info;

// list of spinor indices -- of each symmetry
typedef struct {
    int size;
    int repno;
    int *indices;   // indices of spinors
} spinor_block_t;
extern size_t n_spinor_blocks;
extern spinor_block_t *spinor_blocks;

// all arrays n the following struct contain its size as a first element
typedef struct {
    int *holes;      // including active holes
    int *inact_holes;
    int *act_holes;
    int *act_parts;
    int *inact_parts;
    int *parts;      // including active particles
} spinor_types_t;
extern spinor_types_t spinor_types;

// mapping "global spinor index -> local spinor index (in spinor block)"
extern int spinor_index_global2local[CC_MAX_SPINORS];

// helper functions

int is_active(int idx);

int is_hole(int idx);

int is_act_hole(int idx);

int is_inact_hole(int idx);

int is_part(int idx);

int is_act_part(int idx);

int is_inact_part(int idx);

double get_eps(int idx);

int get_num_electrons();

int get_vacuum_irrep();

void get_active_space_size(int *nacth, int *nactp);

void get_active_space(int *nacth, int *nactp, moindex_t *active_holes_indices, moindex_t *active_parts_indices);

void create_spinor_blocks();

void classify_spinors(double actsp_min, double actsp_max,
                      int nacth, int nactp);

void setup_singles_only_space(int *ccs_only_space);

void print_spinor_info();

int get_spinor_block_size(int spinor_block_number);

int get_max_spinor_block_size();

void spinors_cleanup();

extern int (*is_symblock_zerox[CC_DIAGRAM_MAX_RANK])(int *spinor_blocks, int *qparts, int *valence);

int is_symblock_zero(int rank, int *spinor_blocks, int *qparts, int *valence);

#endif /* CC_SPINORS_H_INCLUDED */
