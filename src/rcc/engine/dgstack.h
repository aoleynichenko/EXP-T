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
 * diagrams are organized as a stack
 */

#ifndef CC_DGSTACK_H_INCLUDED
#define CC_DGSTACK_H_INCLUDED

#include "diagram.h"

typedef int dg_stack_pos_t;

diagram_t *diagram_stack_push(diagram_t *dg);

diagram_t *diagram_stack_replace(char *name, diagram_t *dg);

diagram_t *diagram_stack_find(char *name);

void diagram_stack_erase(char *name);

int diagram_stack_find_index(char *name);

void diagram_stack_print();

dg_stack_pos_t get_stack_pos();

void restore_stack_pos(dg_stack_pos_t pos);

void assert_diagram_exists(char *name);

void rename_diagram(char *old_name, char *new_name);

#endif // CC_DGSTACK_H_INCLUDED
