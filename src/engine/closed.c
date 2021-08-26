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
 * closed.c
 *
 * Extracts closed part of the diagram.
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <assert.h>
#include <complex.h>
#include <string.h>

#include "engine.h"
#include "datamodel.h"
#include "memory.h"
#include "spinors.h"
#include "timer.h"
#include "utils.h"

void clear_valence(diagram_t *dg);


/*******************************************************************************
 * restrict_valence
 *
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
 ******************************************************************************/
void restrict_valence(char *src_name /*large*/, char *tgt_name /*small*/, char *new_valence, int extract_valence)
{
    diagram_t *src, *tgt;
    char qparts[CC_DIAGRAM_MAX_RANK];
    char order[CC_DIAGRAM_MAX_RANK];

    src = diagram_stack_find(src_name);
    if (src == NULL) {
        errquit("restrict_valence(): diagram '%s' not found", src_name);
    }

    if (strlen(new_valence) != src->rank) {
        errquit("restrict_valence(): length of the 'valence' string (%d) must "
                "be equal to the diagrams's rank (%d)", strlen(new_valence), src->rank);
    }

    diagram_get_quasiparticles(src, qparts);
    diagram_get_order(src, order);

    tmplt(tgt_name, qparts, new_valence, order, NOT_PERM_UNIQUE);
    tgt = diagram_stack_find(tgt_name);

    // code is adapted from closed.c
    // allocate temporary working arrays
    size_t maxsize = 0;
    for (int isb1 = 0; isb1 < tgt->n_blocks; isb1++) {
        if (tgt->blocks[isb1]->size > maxsize) {
            maxsize = tgt->blocks[isb1]->size;
        }
    }

    int *tgt_indices = (int *) cc_malloc(sizeof(int) * maxsize * src->rank);

    // for each block in the target (smaller) diagram:
    // 1. find corresponding block in the source (larger) diagram
    // 2. generate indices [i,j,...] for this 'target' (small) block
    // 3. recalculate linear index IDX for [i,j,...] in 'source' (large) diagram
    // 4. target[i,j,...] = source[IDX];

    // loop over blocks in the target diagram
    for (size_t isb2 = 0; isb2 < tgt->n_blocks; isb2++) {
        block_t *sb2 = tgt->blocks[isb2];

        // obtain corresponding block from the target diagram
        block_t *sb1 = diagram_get_block(src, sb2->spinor_blocks);

        if (sb1->is_unique == 0) {
            restore_block(src, sb1);
        }
        symblock_load(sb1);
        symblock_load(sb2);

        symblock_gen_indices(sb2, tgt_indices);

        for (size_t i = 0; i < sb2->size; i++) { // loop over target diagram
            if (carith) {
                sb2->buf[i] = symblock_get(sb1, &tgt_indices[src->rank * i]);
            }
            else {
                double v = creal(symblock_get(sb1, &tgt_indices[src->rank * i]));
                ((double *) sb2->buf)[i] = v;
            }
        }

        symblock_store(sb2);
        if (extract_valence) {
            symblock_store(sb1);
        }
        else {
            symblock_unload(sb1);
        }
        if (sb1->is_unique == 0) {
            destroy_block(sb1);
        }
    }
    if (extract_valence) {
        clear_valence(src);
    }

    // final cleanup
    cc_free(tgt_indices);
}


void clear_valence(diagram_t *dg)
{
    assert(dg != NULL);

    int idx[CC_DIAGRAM_MAX_RANK];
    int dims[CC_DIAGRAM_MAX_RANK];
    int rank = dg->rank;

    for (size_t isb = 0; isb < dg->n_blocks; isb++) {
        block_t *block = dg->blocks[isb];
        size_t size = block->size;

        if (block->is_unique == 0) {
            continue;
        }

        symblock_get_dims(block, dims);
        symblock_load(block);

        // loop over matrix elements
        for (size_t i = 0; i < size; i++) {
            as_compound_index(rank, dims, i /* linear index */, idx);
            // relative indices -> spinor indices
            for (int j = 0; j < rank; j++) {
                idx[j] = block->indices[j][idx[j] + 1];
            }

            // count number of valence indices
            int nvalence = 0;
            for (int j = 0; j < rank; j++) {
                if (is_active(idx[j])) {
                    nvalence++;
                }
            }
            if (nvalence == rank) {
                // set value
                if (carith) {
                    block->buf[i] = 0.0 + 0.0 * I;
                }
                else {
                    ((double *) block->buf)[i] = 0.0;
                }
            }
        }
        symblock_store(block);
    }
}


/*******************************************************************************
 * closed
 *
 * Extracts closed part of the diagram; puts result into 'tgt_name'
 *
 * Arguments:
 *   src_name     symbolic name of the diagram from which closed part must be
 *                extracted
 *   tgt_name     symbolic name of the new diagram -- result
 ******************************************************************************/
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


/*******************************************************************************
 * expand_diagram
 *
 * Expands one or more valence lines of the diagram (valence -> general).
 *
 ******************************************************************************/
void expand_diagram(char *name_small, char *name_large)
{
    diagram_t *dgs = diagram_stack_find(name_small);
    diagram_t *dgl = diagram_stack_find(name_large);

    if (dgs == NULL) {
        errquit("expand_diagram(): diagram '%s' was not found", name_small);
    }
    if (dgl == NULL) {
        errquit("expand_diagram(): diagram '%s' was not found", name_large);
    }

    if (dgs->rank != dgl->rank) {
        errquit("expand_diagram(): diagrams must have equal ranks (%d and %d)", dgs->rank, dgl->rank);
    }
    if (intcmp(dgs->rank, dgs->order, dgl->order) != 0) {
        errquit("expand_diagram(): diagrams must have equal order");
    }
    if (intcmp(dgs->rank, dgs->qparts, dgl->qparts) != 0) {
        errquit("expand_diagram(): diagrams must have equal holes/particles composition");
    }

    for (size_t isb1 = 0; isb1 < dgs->n_blocks; isb1++) {
        block_t *sb1 = dgs->blocks[isb1];
        int *indices = (int *) cc_malloc(sb1->size * sb1->rank * sizeof(int));
        symblock_gen_indices(sb1, indices);

        // obtain corresponding block from the large diagram
        block_t *sb2 = diagram_get_block(dgl, sb1->spinor_blocks);

        if (sb2->is_unique == 0) {
            restore_block(dgl, sb2);
        }

        symblock_load(sb1);
        symblock_load(sb2);

        if (carith) {
            for (size_t i = 0; i < sb1->size; i++) {
                int *idx = &indices[sb1->rank * i];
                double complex val = sb1->buf[i];
                symblock_set(sb2, val, idx);
            }
        }
        else {
            for (size_t i = 0; i < sb1->size; i++) {
                int *idx = &indices[sb1->rank * i];
                double val = ((double *) sb1->buf)[i];
                symblock_set(sb2, val + 0.0 * I, idx);
            }
        }

        symblock_unload(sb1);
        symblock_store(sb2);

        if (sb2->is_unique == 0) {
            destroy_block(sb2);
        }

        cc_free(indices);
    }
}

