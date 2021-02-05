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
 * max.c
 * =====
 *
 * max/diffmax operations for diagrams.
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "engine.h"
#include "datamodel.h"
#include "error.h"
#include "linalg.h"
#include "options.h"
#include "spinors.h"
#include "timer.h"


/*******************************************************************************
 * findmax
 *
 * Find max (abs value) matrix element in the diagram.
 *
 * Arguments:
 *   name     [input]  symbolic name of a new diagram
 *   max_val  [output] max abs value
 *   idx      [output] compound index of max value -- MUST be allocated
 *            (array of size 'rank')
 ******************************************************************************/
void findmax(char *name, double *max_val, int *idx)
{
    diagram_t *dg;

    dg = diagram_stack_find(name);
    if (dg == NULL) {
        errquit("findmax(): diagram '%s' not found", name);
    }

    *max_val = diagram_max(dg, idx);
}


/*******************************************************************************
 * diffmax
 *
 * Find max (abs value) difference between corresponding matrix elements
 * of the two diagrams.
 *
 * Arguments:
 *   name1    [input]  symbolic name of the first diagram
 *   name2    [input]  symbolic name of the second diagram
 *   max_val  [output] max abs value
 *   idx      [output] compound index of max value -- MUST be allocated
 *            (array of size 'rank')
 ******************************************************************************/
void diffmax(char *name1, char *name2, double *max_val, int *idx)
{
    diagram_t *dg1, *dg2;

    dg1 = diagram_stack_find(name1);
    if (dg1 == NULL) {
        errquit("findmax(): diagram-1 '%s' not found", name1);
    }
    dg2 = diagram_stack_find(name2);
    if (dg2 == NULL) {
        errquit("findmax(): diagram-2 '%s' not found", name2);
    }

    *max_val = diagram_diffmax(dg1, dg2, idx);
}


/**
 * Returns absolute value of the max matrix element and its indices.
 * NOTE: array 'indices' of size 'dg->rank' must be allocated and initialized!
 *
 * TODO: test more and more! (especially transform 'linear index' -> 'compound')
 */
double diagram_max(diagram_t *dg, int *indices)
{
    int rank;
    double max_val = 0.0;
    size_t max_idx = 0;
    double max_val_local;
    size_t max_idx_local;
    int dims[CC_DIAGRAM_MAX_RANK];
    block_t *maxblock = NULL;

    rank = dg->rank;
    for (int j = 0; j < rank; j++) {
        indices[j] = 0;
    }

    for (size_t isb = 0; isb < dg->n_blocks; isb++) {
        block_t *block = dg->blocks[isb];
        if (block->is_unique == 0) {
            continue;
        }
        symblock_load(block);

        max_idx_local = ixamax(WORKING_TYPE, block->size, block->buf, &max_val_local);
        if (max_val_local > max_val) {
            max_val = max_val_local;
            max_idx = max_idx_local;
            maxblock = block;
        }

        symblock_unload(block);
    }   // end loop over blocks

    // recalculate index of max value from the "linear" form to the "compound"
    if (maxblock != NULL) {
        // get block's dimension
        symblock_get_dims(maxblock, dims);
        // transform 'linear' -> 'compound'
        as_compound_index(rank, dims, max_idx, indices);
        // 'relative' indices into spinor indices
        for (int i = 0; i < rank; i++) {
            indices[i] = maxblock->indices[i][indices[i] + 1];
        }
    }

    return max_val;
}


/**
 * Returns absolute value of the max difference between two diagrams
 * (elementwise) and its indices.
 * NOTE: array 'indices' of size 'dg->rank' must be allocated and initialized!
 *
 * TODO: test more and more! (especially transform 'linear index' -> 'compound')
 */
double diagram_diffmax(diagram_t *dg1, diagram_t *dg2, int *indices)
{
    size_t isb1, isb2;
    block_t *block1, *block2, *maxblock = NULL;
    int rank;
    double max_val = 0.0;
    double max_val_local;
    size_t max_idx = 0, max_idx_local;
    int dims[CC_DIAGRAM_MAX_RANK];

    rank = dg1->rank;
    for (int j = 0; j < rank; j++) {
        indices[j] = 0;
    }

    for (isb1 = 0; isb1 < dg1->n_blocks; isb1++) {
        block1 = dg1->blocks[isb1];
        if (block1->is_unique == 0) {
            continue;
        }

        // get block2 via its spinor_blocks array
        diagram_get_block(dg2, block1->spinor_blocks, &isb2);
        block2 = dg2->blocks[isb2];

        symblock_load(block1);
        symblock_load(block2);

        max_idx_local = ixadiffmax(WORKING_TYPE,
                                   block1->size, block1->buf,
                                   block2->buf, &max_val_local);
        if (max_val_local > max_val) {
            max_val = max_val_local;
            max_idx = max_idx_local;
            maxblock = block1;
        }

        symblock_unload(block1);
        symblock_unload(block2);
    }

    // recalculate index of max value from the "linear" form to the "compound"
    if (maxblock != NULL) {
        // get block's dimension
        symblock_get_dims(maxblock, dims);
        // transform 'linear' -> 'compound'
        as_compound_index(rank, dims, max_idx, indices);
        // 'relative' indices into spinor indices
        for (int i = 0; i < rank; i++) {
            indices[i] = maxblock->indices[i][indices[i] + 1];
        }
    }

    return max_val;
}
