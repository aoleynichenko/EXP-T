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

#include <string.h>

#include "engine.h"
#include "options.h"
#include "symmetry.h"
#include "timer.h"


/**
 * Creates empty diagram and adds it to the diagram stack.
 * The diagram will be initialized with zeros.
 *
 * Arguments:
 *   name     symbolic name of a new diagram
 *   qparts   holes/particles list (hh, hp, hhpp, pppp, etc)
 *   valence  which indices are valence (1/0) ('10', '0000', '1100', etc)
 *   order    initial order of tensor dimensions ('1234', '1243', etc)
 *   irrep    irrep of the operator represented by a diagram
 */
void tmplt_sym(char *name, char *qparts, char *valence, char *order, int perm_unique, int irrep)
{
    timer_new_entry("tmplt", "Diagram template constr (tmplt)");
    timer_start("tmplt");

    // prepare the 't3space' string
    int rank = strlen(qparts);
    char t3space[64];
    for (int i = 0; i < rank; i++) {
        if (rank == 6 && cc_opts->do_restrict_t3) {
            t3space[i] = '1';
        }
        else {
            t3space[i] = '0';
        }
    }
    t3space[rank] = '\0';

    // remove old diagram if presented
    diagram_t *dg = diagram_stack_find(name);

    if (dg != NULL) {
        dg = diagram_new(name, qparts, valence, t3space, order, perm_unique, irrep);
        diagram_stack_replace(name, dg);
    }
    else {
        dg = diagram_new(name, qparts, valence, t3space, order, perm_unique, irrep);
        diagram_stack_push(dg);
    }

    timer_stop("tmplt");
}


/**
 * specific case of a totally symmetric operator
 */
void tmplt(char *name, char *qparts, char *valence, char *order, int perm_unique)
{
    tmplt_sym(name, qparts, valence, order, perm_unique, get_totally_symmetric_irrep());
}
