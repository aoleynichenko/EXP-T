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
 * Interface to the diagram engine. It includes only some top-level routines,
 * so the user is not required to think about diagrams' internal structure and
 * algorithms used.
 */

#ifndef CC_ENGINE_H_INCLUDED
#define CC_ENGINE_H_INCLUDED

#include <complex.h>

// unified interface to platform-specific functions
#include "platform.h"

// functions for time measurements
#include "timer.h"

// error handling
#include "error.h"

// memory allocator
#include "memory.h"

#include "options.h"

// data structures:
#include "../engine/block.h"    // blocks
#include "../engine/diagram.h"  // diagrams
#include "../engine/dgstack.h"  // stack of diagrams

enum {
    NOT_PERM_UNIQUE = 0,
    IS_PERM_UNIQUE = 1
};

// returns rank of a diagram -- number of dimensions of its tensors
int rank(char *name);

// creates empty diagram and adds it to the diagram stack
void tmplt(char *name, char *qparts, char *valence, char *order, int perm_unique);
void tmplt_sym(char *name, char *qparts, char *valence, char *order, int perm_unique, int irrep);

// sets all matrix elements in the diagram to zero
void clear(char *name);

// find max (abs value) matrix element in the diagram;
void findmax(char *name, double *max_val, int *idx);

// find max (abs value) difference between corresponding matrix elements
void diffmax(char *name1, char *name2, double *max_val, int *idx);

// prints given diagram to stdout
void prt(char *name);

// brief summary about the diagram
void summary(char *name);

// creates a copy of a diagram
void copy(char *src, char *target);

// interchanges vertices in the two-electron diagram
void interchange_electrons(char *src_name, char *dst_name);

// in-place complex conjugation of all matrix elements in the diagram
void conjugate_elements_inplace(char *source_name);

// reorders diagram 'src_diargam_name'; puts result into 'target_diagram_name'
void reorder(char *src_diargam_name, char *target_diagram_name, char *perm_str);

// performs contraction of two diagrams
void mult(char *name1, char *name2, char *target, int ncontr);

void tt_enable();
void tt_disable();

// scalar product of two diagrams (full contraction)
// (mult.c)
double complex scalar_product(char *conj1, char *conj2, char *name1, char *name2);

// performs (elementwise) addition of two diagrams
void add(double fact1, char *name1, double fact2, char *name2, char *target);

// Updates diagram dg1_name: dg1 += factor * dg2
void update(char *dg1_name, double factor, char *dg2_name);

// performs index permutation for the diagram
void perm(char *dg_name, char *perm_str);

// division by energy denominators
void diveps(char *name);

// extracts closed part of the diagram
void closed(char *src_name, char *tgt_name);

// restricts some lines from general to valence
void restrict_valence(char *src_name /*large*/, char *tgt_name /*small*/, char *new_valence, int extract_valence);;

// get symmetry information
void get_rep_name(int irep, char *name);

// get spinor information (by spinor number)
void get_spinor_info(int idx, int *repno, int *occ, int *active, double *eps);

int get_spinor_irrep(int idx);

void expand_diagram(char *name_small, char *name_large);

void restrict_triples(char *dg_name, double erange_1, double erange_2);

void remove_core_correlation(char *diag_name);

void remove_inner_core_correlation(char *diag);

//void restore_diagram(char *dg_name);

void clear_non_unique(char *name);

void check_unique(char *name);

void predict_intruders_for_diagram(char *diagram_name, int nmax);

void print_ampl_vs_denom(char *name, char *out_file);

size_t count_amplitudes(char *diagram_name);

size_t count_amplitudes_in_range(char *diagram_name, double lower_bound, double upper_bound);

void print_amplitude_distribution_analysis(char *diagram_name);

void diagram_conjugate(char *source_name, char *target_name);

#include "../engine/disconnected.h"
#include "../engine/tensor_trains.h"

#endif /* CC_ENGINE_H_INCLUDED */
