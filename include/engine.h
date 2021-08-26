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
 * engine.h
 * ========
 *
 * Interface to the diagram engine. It includes only some top-level routines,
 * so the user is not required to think about diagrams' internal structure and
 * algorithms used.
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

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

enum {
    NOT_PERM_UNIQUE = 0,
    IS_PERM_UNIQUE = 1
};

// returns rank of a diagram -- number of dimensions of its tensors
int rank(char *name);

// creates empty diagram and adds it to the diagram stack
void tmplt(char *name, char *qparts, char *valence, char *order, int perm_unique);

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

// renames diagram in the stack
void rename_diagram(char *old_name, char *new_name);

// creates a copy of a diagram
void copy(char *src, char *target);

// reorders diagram 'src'; puts result into 'target'
void reorder(char *src, char *target, char *perm_str);

// performs contraction of two diagrams
void mult(char *name1, char *name2, char *target, int ncontr);

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

void predict_intruders(char *diagram_name, int nmax);

void print_ampl_vs_denom(char *name, char *out_file);

#endif /* CC_ENGINE_H_INCLUDED */
