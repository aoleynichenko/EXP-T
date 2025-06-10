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
 * Extracts closed part of the diagram.
 */

#include <assert.h>
#include <complex.h>
#include <string.h>

#include "engine.h"
#include "memory.h"
#include "spinors.h"
#include "timer.h"
#include "utils.h"

void clear_valence(diagram_t *dg);


/**
 * Turns the diagram's "general" (valence + inactive) lines (index ranges) into
 * valence ones.
 *
 * 'extract_valence': if 1, all amplitudes "extracted" from the 'src_name'
 * to the 'tgt_name' diagram will be set to zero in the source diagram.
 * This make the code more universal.
 * extract_valence = 0: integrals are copied from src to tgt
 * extract_valence = 1: the same + remove all integrals which were copied
 *
 * Example: new_valence = 1100
 *           "pphp" -> "vvgp"
 * [valence]  0000      1110
 */
void restrict_valence(char *src_name /*large*/, char *tgt_name /*small*/, char *new_valence, int extract_valence)
{
    char qparts[CC_DIAGRAM_MAX_RANK];
    char order[CC_DIAGRAM_MAX_RANK];

    assert_diagram_exists(src_name);
    diagram_t *src = diagram_stack_find(src_name);
    int rank = src->rank;

    if (strlen(new_valence) != rank) {
        errquit("restrict_valence(): length of the 'valence' string (%d) must "
                "be equal to the diagrams's rank (%d)", strlen(new_valence), rank);
    }

    diagram_get_quasiparticles(src, qparts);
    diagram_get_order(src, order);

    tmplt_sym(tgt_name, qparts, new_valence, order, NOT_PERM_UNIQUE, src->symmetry);
    diagram_t *tgt = diagram_stack_find(tgt_name);

    // allocate temporary working arrays
    size_t maxsize = 0;
    for (int isb1 = 0; isb1 < tgt->n_blocks; isb1++) {
        if (tgt->blocks[isb1]->size > maxsize) {
            maxsize = tgt->blocks[isb1]->size;
        }
    }

    int *tgt_indices = (int *) cc_malloc(sizeof(int) * maxsize * rank);

    // for each block in the target (smaller) diagram:
    // 1. find corresponding block in the source (larger) diagram
    // 2. generate indices [i,j,...] for this 'target' (small) block
    // 3. recalculate linear index IDX for [i,j,...] in 'source' (large) diagram
    // 4. target[i,j,...] = source[IDX];

    // loop over blocks in the target diagram
    for (size_t isb2 = 0; isb2 < tgt->n_blocks; isb2++) {
        block_t *target_block = tgt->blocks[isb2];

        // obtain corresponding block from the source diagram
        block_t *source_block = diagram_get_block(src, target_block->spinor_blocks);

        if (source_block->is_unique == 0) {
            restore_block(src, source_block);
        }
        block_load(source_block);
        block_load(target_block);

        block_gen_indices(target_block, tgt_indices);

        for (size_t i = 0; i < target_block->size; i++) { // loop over target diagram
            int *idx = tgt_indices + rank * i;

            if (arith == CC_ARITH_COMPLEX) {
                target_block->buf[i] = block_get_element(source_block, idx);
            }
            else {
                double v = creal(block_get_element(source_block, idx));
                ((double *) target_block->buf)[i] = v;
            }
        }

        block_store(target_block);
        if (extract_valence) {
            block_store(source_block);
        }
        else {
            block_unload(source_block);
        }
        if (source_block->is_unique == 0) {
            destroy_block(source_block);
        }
    }

    if (extract_valence) {
        clear_valence(src);
    }

    // final cleanup
    cc_free(tgt_indices);
}


/**
 * Sets matrix elements with all indices are active (valence) to zero.
 */
void clear_valence(diagram_t *dg)
{
    assert(dg != NULL);

    int idx[CC_DIAGRAM_MAX_RANK];
    int rank = dg->rank;

    for (size_t iblock = 0; iblock < dg->n_blocks; iblock++) {
        block_t *block = dg->blocks[iblock];
        size_t size = block->size;

        if (block->is_unique == 0) {
            continue;
        }

        block_load(block);

        // loop over matrix elements
        for (size_t i = 0; i < size; i++) {
            tensor_index_to_compound(rank, block->shape, i /* linear index */, idx);

            // count number of valence indices
            int nvalence = 0;
            for (int j = 0; j < rank; j++) {
                int spinor_index = block->indices[j][idx[j]]; // relative index -> spinor index
                if (is_active(spinor_index)) {
                    nvalence++;
                }
            }

            if (nvalence == rank) {
                if (arith == CC_ARITH_COMPLEX) {
                    block->buf[i] = 0.0 + 0.0 * I;
                }
                else {
                    ((double *) block->buf)[i] = 0.0;
                }
            }
        }
        block_store(block);
    }
}


/**
 * Extracts closed part of the diagram; puts result into 'tgt_name'
 *
 * Arguments:
 *   src_name     symbolic name of the diagram from which closed part must be
 *                extracted
 *   tgt_name     symbolic name of the new diagram -- result
 */
void closed(char *src_name, char *tgt_name)
{
    char valence[CC_DIAGRAM_MAX_RANK];

    timer_new_entry("closed", "Extraction of a closed part");
    timer_start("closed");

    diagram_t *dg_src = diagram_stack_find(src_name);
    if (dg_src == NULL) {
        errquit("closed(): diagram '%s' not found", src_name);
    }

    // create a new diagram with the same list of quasiparticles and order,
    // but with indices which all are valence
    int rank = dg_src->rank;
    for (int i = 0; i < rank; i++) {
        valence[i] = '1';
    }
    valence[rank] = '\0';

    restrict_valence(src_name, tgt_name, valence, 1);

    timer_stop("closed");
}


/**
 * Expands one or more valence lines of the diagram (valence -> general).
 */
void expand_diagram(char *name_small, char *name_large)
{
    assert_diagram_exists(name_small);
    assert_diagram_exists(name_large);

    diagram_t *dg_small = diagram_stack_find(name_small);
    diagram_t *dg_large = diagram_stack_find(name_large);

    if (dg_small->rank != dg_large->rank) {
        errquit("expand_diagram(): diagrams must have equal ranks (%d and %d)", dg_small->rank, dg_large->rank);
    }
    if (intcmp(dg_small->rank, dg_small->order, dg_large->order) != 0) {
        errquit("expand_diagram(): diagrams must have equal order");
    }
    if (intcmp(dg_small->rank, dg_small->qparts, dg_large->qparts) != 0) {
        errquit("expand_diagram(): diagrams must have equal holes/particles composition");
    }

    for (size_t ib1 = 0; ib1 < dg_small->n_blocks; ib1++) {
        block_t *source_block = dg_small->blocks[ib1];
        int *indices = (int *) cc_malloc(source_block->size * source_block->rank * sizeof(int));
        block_gen_indices(source_block, indices);

        // obtain corresponding block from the large diagram
        block_t *target_block = diagram_get_block(dg_large, source_block->spinor_blocks);

        if (target_block->is_unique == 0) {
            restore_block(dg_large, target_block);
        }

        block_load(source_block);
        block_load(target_block);

        if (arith == CC_ARITH_COMPLEX) {
            for (size_t i = 0; i < source_block->size; i++) {
                int *idx = indices + source_block->rank * i;
                double complex val = source_block->buf[i];
                block_set_element(target_block, val, idx);
            }
        }
        else {
            for (size_t i = 0; i < source_block->size; i++) {
                int *idx = indices + source_block->rank * i;
                double val = ((double *) source_block->buf)[i];
                block_set_element(target_block, val + 0.0 * I, idx);
            }
        }

        block_unload(source_block);
        block_store(target_block);

        if (target_block->is_unique == 0) {
            destroy_block(target_block);
        }

        cc_free(indices);
    }
}

