/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2020 The EXP-T developers.
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
 * @file info_queries.c
 *
 * Info queries to diagrams: rank, size, valence, quasiparticles etc.
 * It is a part of the user-friendly interface to diagrams which does not
 * require to know anything about the data model employed.
 *
 * @author 2018-2020 Alexander Oleynichenko
 ******************************************************************************/

#include "engine.h"

#include <stdio.h>

#include "error.h"
#include "datamodel.h"


/**
 * Returns rank of a diagram -- number of dimensions of its tensors.
 *
 * @param name symbolic name of a diagram
 * @return rank of the diagram
 */
int rank(char *name)
{
    diagram_t *dg;

    dg = diagram_stack_find(name);
    if (dg == NULL) {
        errquit("rank(): diagram '%s' not found", name);
    }

    return dg->rank;
}


/**
 * prints brief information about the diagram
 *
 * @param name of the diagram
 */
void summary(char *name)
{
    diagram_t *diag;

    diag = diagram_stack_find(name);
    if (diag == NULL) {
        errquit("diagram '%s' not found", name);
    }

    diagram_summary(diag);
}


/**
 * prints given diagram to stdout.
 *
 * Arguments:
 *   name    symbolic name of the diagram to be printed
 */
void prt(char *name)
{
    diagram_t *dg;

    dg = diagram_stack_find(name);
    if (dg == NULL) {
        errquit("prt(): diagram '%s' not found", name);
    }

    diagram_print(dg);
}
