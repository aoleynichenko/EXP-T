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
 * relaxation.c
 * ============
 *
 * Sets to zero all Doubles amplitudes containing indices of spinors which
 * were classified as those to be treated only by at the CCS level.
 * (in order to clarify the role of relaxation in properties calculations).
 *
 * 2020-2021 Alexander Oleynichenko
 ******************************************************************************/

#include "datamodel.h"
#include "engine.h"
#include "error.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"


/**
 * sets to zero all amplitudes with TWO HOLE INDICES with eps < inner_core_thresh
 */
void remove_inner_core_correlation(char *diag_name)
{
    diagram_t *diag;
    int idx[CC_DIAGRAM_MAX_RANK];
    int dims[CC_DIAGRAM_MAX_RANK];
    int rank;

    diag = diagram_stack_find(diag_name);
    if (diag == NULL) {
        errquit("remove_core_correlation(): diagram '%s' not found");
        return;
    }

    rank = diag->rank;
    for (size_t isb = 0; isb < diag->n_blocks; isb++) {
        block_t *block = diag->blocks[isb];
        size_t size = block->size;
        if (block->is_unique == 0) {
            continue;
        }
        symblock_get_dims(block, dims);
        symblock_load(block);
        double complex *buf = block->buf;

        // loop over matrix elements
        for (size_t i = 0; i < size; i++) {
            int n_inner_holes = 0;
            as_compound_index(rank, dims, i /* linear index */, idx);
            // relative indices -> spinor indices
            for (size_t j = 0; j < rank; j++) {
                idx[j] = block->indices[j][idx[j] + 1];
                if (spinor_info[idx[j]].eps <= cc_opts->inner_core_thresh) {
                    n_inner_holes++;
                }
            }
            if (n_inner_holes == 2) {
                if (carith) {
                    buf[i] = 0.0 + 0.0 * I;
                }
                else {
                    ((double *) buf)[i] = 0.0;
                }
            }
            /*for (size_t j = 0; j < rank; j++) {
                printf("%2d ", idx[j]);
            }
            printf("  :  %d\n", n_inner_holes);*/
        }

        symblock_store(block);
    }
}


/**
 * Sets to zero all amplitudes corresponding to the "core" and "core valence"
 * ranges of molecular spinors.
 * => formally, only CCS is used for such spinors.
 * To be used in order to estimate contribution of purely relaxation effects
 * to properties.
 */
void remove_core_correlation(char *diag_name)
{
    diagram_t *diag;
    int idx[CC_DIAGRAM_MAX_RANK];
    int dims[CC_DIAGRAM_MAX_RANK];
    int ccs_only_space[CC_MAX_SPINORS];
    int rank;

    diag = diagram_stack_find(diag_name);
    if (diag == NULL) {
        errquit("remove_core_correlation(): diagram '%s' not found");
        return;
    }

    setup_singles_only_space(ccs_only_space);

    rank = diag->rank;
    for (size_t isb = 0; isb < diag->n_blocks; isb++) {
        block_t *block = diag->blocks[isb];
        size_t size = block->size;
        symblock_get_dims(block, dims);
        symblock_load(block);
        double complex *buf = block->buf;

        // loop over matrix elements
        for (size_t i = 0; i < size; i++) {
            as_compound_index(rank, dims, i /* linear index */, idx);
            // relative indices -> spinor indices
            for (size_t j = 0; j < rank; j++) {
                idx[j] = block->indices[j][idx[j] + 1];
                if (ccs_only_space[idx[j]]) {
                    if (carith) {
                        buf[i] = 0.0 + 0.0 * I;
                    }
                    else {
                        ((double *) buf)[i] = 0.0;
                    }
                    break;
                }
            }
        }

        symblock_store(block);
    }
}
