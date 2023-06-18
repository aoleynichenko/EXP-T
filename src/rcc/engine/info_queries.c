/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2023 The EXP-T developers.
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
 * Info queries to diagrams: rank, size, valence, quasiparticles etc.
 * It is a part of the user-friendly interface to diagrams which does not
 * require to know anything about the data model employed.
 */

#include "engine.h"

#include <math.h>
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


/**
 * returns total number of amplitudes stored.
 */
size_t count_amplitudes(char *diagram_name)
{
    size_t num_amplitudes = 0;

    diagram_t *diag = diagram_stack_find(diagram_name);
    if (diag == NULL) {
        errquit("count_amplitudes(): diagram '%s' is not found", diagram_name);
    }

    for (size_t iblock = 0; iblock < diag->n_blocks; iblock++) {
        block_t *block = diag->blocks[iblock];

        // only unique blocks are to be analyzed
        if (block->is_unique == 0) {
            continue;
        }

        num_amplitudes += block->size;
    }

    return num_amplitudes;
}



/**
 * counts amplitudes with absolute values in the given range.
 */
size_t count_amplitudes_in_range(char *diagram_name, double lower_bound, double upper_bound)
{
    size_t num_amplitudes = 0;

    diagram_t *diag = diagram_stack_find(diagram_name);
    if (diag == NULL) {
        errquit("count_amplitudes_in_range(): diagram '%s' is not found", diagram_name);
    }

    for (size_t iblock = 0; iblock < diag->n_blocks; iblock++) {
        block_t *block = diag->blocks[iblock];

        // only unique blocks are to be analyzed
        if (block->is_unique == 0) {
            continue;
        }

        block_load(block);

        for (size_t i = 0; i < block->size; i++) {
            double complex t = (arith == CC_ARITH_COMPLEX) ? block->buf[i] : ((double *) block->buf)[i] + 0.0 * I;
            double abs_t = cabs(t);
            if (lower_bound <= abs_t && abs_t <= upper_bound) {
                num_amplitudes++;
            }
        }

        block_unload(block);
    }

    return num_amplitudes;
}


void print_amplitude_distribution_analysis(char *diagram_name)
{
    double upper_bound = 1.0;
    double lower_bound = 0.1;

    size_t num_total = count_amplitudes(diagram_name);

    printf("\n");
    printf(" diagram: %s\n", diagram_name);
    printf(" total: %ld amplitudes\n", num_total);

    while (upper_bound >= 1e-15) {
        size_t num_in_range = count_amplitudes_in_range(diagram_name, lower_bound, upper_bound);

        printf("  [%6.0e :%6.0e ]%16ld  %4.1f%%\n", lower_bound, upper_bound, num_in_range, 100.0 * ((double) num_in_range) / num_total);

        upper_bound /= 10;
        lower_bound /= 10;
    }

    // all remaining amplitudes
    lower_bound = 0.0;
    size_t num_in_range = count_amplitudes_in_range(diagram_name, lower_bound, upper_bound);
    printf("  [%6.0e :%6.0e ]%16ld  %4.1f%%\n", lower_bound, upper_bound, num_in_range, 100.0 * ((double) num_in_range) / num_total);
    printf("\n");
}
