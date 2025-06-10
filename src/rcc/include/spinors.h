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
 * All information about one-particle functions (spin-orbitals or spinors)
 * (and symmetry).
 */

#ifndef CC_SPINORS_H_INCLUDED
#define CC_SPINORS_H_INCLUDED

#include "comdef.h"
#include "options.h"

//#define CC_FLAG_OCCUPIED  0x00000001
#define CC_FLAG_ACTIVE    0x00000002
#define CC_FLAG_MAIN      0x00000004
#define CC_FLAG_T3_SPACE  0x00000008

// attributes of each spinors:
//  - sequential number
//  - irrep no
//  - hole / particle (occupied in reference or not) (1 - hole, 0 - particle)
//  - space_flags => active/inactive, main/intermediate
//  - its energy
typedef struct {
    int repno;
    int blockno;
    double eps;
    int occ;
    uint32_t space_flags;
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

// mapping "global spinor index -> local spinor index (in spinor block)"
extern int spinor_index_global2local[CC_MAX_SPINORS];

void set_occupied(int idx, int occ_number);

int get_num_spinors();

int is_active(int idx);

int is_hole(int idx);

int is_act_hole(int idx);

int is_inact_hole(int idx);

int is_particle(int idx);

int is_act_particle(int idx);

int is_inact_particle(int idx);

int is_t3_space_spinor(int idx);

double get_eps(int idx);

void get_spinor_energies(size_t n, double *eps);

int get_num_electrons();

void get_spinor_indices_occupied(int *occ_idx);

int get_vacuum_irrep();

void get_active_space_size(int *nacth, int *nactp);

void get_active_holes_particles(int *nacth, int *nactp, moindex_t *active_holes_indices, moindex_t *active_parts_indices);

void get_active_space(int sect_h, int sect_p, int *n_active, int *active_spinors);

void create_spinor_info(int n_spinors, int *irrep_numbers, double *spinor_energies, int *occ_numbers);

void create_spinor_blocks(int tilesize);

void setup_occupation_numbers(cc_options_t *options, int num_spinors, spinor_attr_t *spinor_info);

void setup_active_space(cc_options_t *cc_opts);

void setup_fast_access_spinor_lists();

void print_spinor_info_table();

int get_spinor_block_size(int spinor_block_number);

int get_max_spinor_block_size();

void spinors_cleanup();

int is_symblock_zero(int rank, int *spinor_blocks, int *qparts, int *valence, int *t3space);

#endif /* CC_SPINORS_H_INCLUDED */
