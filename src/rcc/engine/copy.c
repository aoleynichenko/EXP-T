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


/**
 * Creates a copy of a diagram.
 *
 * Arguments:
 *   src      symbolic name of the diagram to be copied
 *   target   symbolic name of the new diagram -- copy
 */
void copy(char *src, char *target)
{
    assert_diagram_exists(src);
    diagram_t *dg_src = diagram_stack_find(src);

    // copy
    diagram_t *dg_tgt = diagram_copy(dg_src);
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
