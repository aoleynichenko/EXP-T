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
 * basic.c
 * =======
 *
 * Basic operations with diagrams
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "engine.h"
#include "datamodel.h"
#include "error.h"
#include "options.h"
#include "spinors.h"
#include "timer.h"


/*******************************************************************************
 * tmplt
 *
 * Creates empty diagram and adds it to the diagram stack.
 * The diagram will be initialized with zeros.
 *
 * Arguments:
 *   name     symbolic name of a new diagram
 *   qparts   holes/particles list (hh, hp, hhpp, pppp, etc)
 *   valence  which indices are valence (1/0) ('10', '0000', '1100', etc)
 *   order    initial order of tensor dimensions ('1234', '1243', etc)
 ******************************************************************************/
void tmplt(char *name, char *qparts, char *valence, char *order, int perm_unique)
{
    diagram_t *dg;

    timer_new_entry("tmplt", "Diagram template constr (tmplt)");
    timer_start("tmplt");

    // remove old diagram if presented
    // TODO: test it
    dg = diagram_stack_find(name);

    if (dg != NULL) {
        dg = diagram_new(name, qparts, valence, order, perm_unique);
        diagram_stack_replace(name, dg);
    }
    else {
        dg = diagram_new(name, qparts, valence, order, perm_unique);
        diagram_stack_push(dg);
    }

    timer_stop("tmplt");
}


void check_uniquex(char *name)
{
    diagram_t *dg = diagram_stack_find(name);
    assert(dg != NULL);

    if (dg->only_unique == 1) {
        return;
    }
    dg->only_unique = 1;

    for (size_t i = 0; i < dg->n_blocks; i++) {
        block_t *b = dg->blocks[i];
        block_unique(b, dg->qparts, dg->valence, dg->order);
    }
}


void clear_non_uniquex(char *name)
{
    diagram_t *dg = diagram_stack_find(name);
    assert(dg != NULL);

    for (size_t i = 0; i < dg->n_blocks; i++) {
        block_t *b = dg->blocks[i];
        if (b->is_unique) {
            continue;
        }
        memset(b->buf, 0, SIZEOF_WORKING_TYPE * b->size);
    }
}


/*******************************************************************************
 * clear
 *
 * Sets all matrix elements in the diagram to zero.
 *
 * Arguments:
 *   name     symbolic name of a new diagram
 ******************************************************************************/
void clear(char *name)
{
    diagram_t *dg;

    dg = diagram_stack_find(name);
    if (dg == NULL) {
        errquit("clear(): diagram '%s' not found", name);
    }

    diagram_clear(dg);
}


/*******************************************************************************
 * copy
 *
 * Creates a copy of a diagram.
 *
 * Arguments:
 *   src      symbolic name of the diagram to be copied
 *   target   symbolic name of the new diagram -- copy
 ******************************************************************************/
void copy(char *src, char *target)
{
    diagram_t *dg_src, *dg_tgt;

    dg_src = diagram_stack_find(src);
    if (dg_src == NULL) {
        errquit("copy(): diagram '%s' not found", src);
    }

    // copy
    dg_tgt = diagram_copy(dg_src);
    strcpy(dg_tgt->name, target);

    // save new diagram to stack
    // (replace the old diagram named 'target' if needed)
    if (diagram_stack_find(target) != NULL) {
        diagram_stack_replace(target, dg_tgt);
    }
    else {
        diagram_stack_push(dg_tgt);
    }
}


/**
 * Renames diagram in the stack.
 */
void rename_diagram(char *old_name, char *new_name)
{
    diagram_t *dg = diagram_stack_find(old_name);
    if (dg == NULL) {
        errquit("rename_diagram(): diagram '%s' not found", old_name);
    }

    strcpy(dg->name, new_name);
}
