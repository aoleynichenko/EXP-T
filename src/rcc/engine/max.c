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
 * max/diffmax operations for diagrams.
 */

#include "engine.h"
#include "linalg.h"
#include "options.h"


/**
 * Find max (abs value) matrix element in the diagram.
 *
 * Arguments:
 *   name     [input]  symbolic name of a new diagram
 *   max_val  [output] max abs value
 *   idx      [output] compound index of max value -- MUST be allocated
 *            (array of size 'rank')
 */
void findmax(char *name, double *max_val, int *idx)
{
    assert_diagram_exists(name);
    diagram_t *dg = diagram_stack_find(name);
    *max_val = diagram_max(dg, idx);
}


/**
 * Find max (abs value) difference between corresponding matrix elements
 * of the two diagrams.
 *
 * Arguments:
 *   name1    [input]  symbolic name of the first diagram
 *   name2    [input]  symbolic name of the second diagram
 *   max_val  [output] max abs value
 *   idx      [output] compound index of max value -- MUST be allocated
 *            (array of size 'rank')
 */
void diffmax(char *name1, char *name2, double *max_val, int *idx)
{
    assert_diagram_exists(name1);
    assert_diagram_exists(name2);

    diagram_t *dg1 = diagram_stack_find(name1);
    diagram_t *dg2 = diagram_stack_find(name2);

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
    double max_val = 0.0;
    size_t max_idx = 0;
    block_t *max_block = NULL;

    int rank = dg->rank;
    for (int j = 0; j < rank; j++) {
        indices[j] = 0;
    }

    for (size_t isb = 0; isb < dg->n_blocks; isb++) {
        block_t *block = dg->blocks[isb];
        if (block->is_unique == 0) {
            continue;
        }
        block_load(block);

        double max_val_local;
        size_t max_idx_local = max_idx_local = ixamax(WORKING_TYPE, block->size, block->buf, &max_val_local);
        if (max_val_local > max_val) {
            max_val = max_val_local;
            max_idx = max_idx_local;
            max_block = block;
        }

        block_unload(block);
    }   // end loop over blocks

    // recalculate index of max value from the "linear" form to the "compound"
    if (max_block != NULL) {
        // transform 'linear' -> 'compound'
        tensor_index_to_compound(rank, max_block->shape, max_idx, indices);

        // 'relative' indices into spinor indices
        for (int i = 0; i < rank; i++) {
            indices[i] = max_block->indices[i][indices[i]];
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
    double max_val = 0.0;
    size_t max_idx = 0;
    block_t *max_block = NULL;

    int rank = dg1->rank;
    for (int j = 0; j < rank; j++) {
        indices[j] = 0;
    }

    for (size_t isb1 = 0; isb1 < dg1->n_blocks; isb1++) {
        block_t *block1 = dg1->blocks[isb1];
        if (block1->is_unique == 0) {
            continue;
        }

        // get block2 via its spinor_blocks array
        block_t *block2 = diagram_get_block(dg2, block1->spinor_blocks);

        block_load(block1);
        block_load(block2);

        double max_val_local;
        size_t max_idx_local = ixadiffmax(WORKING_TYPE,
                                   block1->size, block1->buf,
                                   block2->buf, &max_val_local);

        if (max_val_local > max_val) {
            max_val = max_val_local;
            max_idx = max_idx_local;
            max_block = block1;
        }

        block_unload(block1);
        block_unload(block2);
    }

    // recalculate index of max value from the "linear" form to the "compound"
    if (max_block != NULL) {
        // transform 'linear' -> 'compound'
        tensor_index_to_compound(rank, max_block->shape, max_idx, indices);

        // 'relative' indices into spinor indices
        for (int i = 0; i < rank; i++) {
            indices[i] = max_block->indices[i][indices[i]];
        }
    }

    return max_val;
}
