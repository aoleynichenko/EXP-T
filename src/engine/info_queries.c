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
#include "spinors.h"


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


void print_ampl_vs_denom(char *name, char *out_file)
{
    diagram_t *dg;
    FILE *f = fopen(out_file, "w");

    dg = diagram_stack_find(name);
    if (dg == NULL) {
        errquit("print_ampl_vs_denom(): diagram '%s' not found", name);
    }

    int rank = dg->rank;

    for (size_t ib = 0; ib < dg->n_blocks; ib++) {
        if (dg->blocks[ib]->is_unique == 0) {
            continue;
        }

        block_t *block = dg->blocks[ib];

        size_t i, j;
        int *indices = (int *) cc_malloc(block->size * block->rank * sizeof(int));

        symblock_load(block);
        symblock_gen_indices(block, indices);

        for (i = 0; i < block->size; i++) {
            double t;
            int *idx = indices + block->rank * i;
            if (arith == CC_ARITH_COMPLEX) {
                t = cabs(block->buf[i]);
            }
            else {
                t = fabs(((double *)block->buf)[i]);
            }

            if (t >= cc_opts->conv) {
                double denom = 0.0;
                double esum1 = 0.0;
                double esum2 = 0.0;
                for (j = 0; j < block->rank/2; j++) {
                    double eps = spinor_info[idx[j]].eps;
                    esum1 += eps;
                }
                for (j = block->rank/2; j < block->rank; j++) {
                    double eps = spinor_info[idx[j]].eps;
                    esum2 -= eps;
                }
                denom = esum1 + esum2;
                for (int j = 0; j < block->rank; j++) {
                    fprintf(f, "%5d", idx[j]);
                }
                fprintf(f,"%16.8f%16.8f%16.8f%20.12f\n", esum1, esum2, -denom, t);
            }
        }

        symblock_unload(block);

        cc_free(indices);

    }

    fclose(f);

}
