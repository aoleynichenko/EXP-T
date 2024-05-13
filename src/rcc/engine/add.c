/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2024 The EXP-T developers.
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
 * Subroutines for addition of diagrams
 */

#include <stdio.h>
#include <string.h>

#include "engine.h"
#include "linalg.h"
#include "options.h"
#include "utils.h"

#ifndef COMPILER_CLANG
#include "omp.h"
#else

void omp_set_num_threads(int num_threads);

void openblas_set_num_threads(int num_threads);

#endif


/**
 * Performs (elementwise) addition of two diagrams:
 * target = fact1 * name1 + fact2 * name2
 * 'target', 'name1', 'name2' are allowed to coincide.
 *
 * Arguments:
 *   fact1    numerical factor for the first diagram
 *   name1    symbolic name of the diagram -- first term
 *   fact2    numerical factor for the second diagram
 *   name2    symbolic name of the diagram -- second term
 *   target   symbolic name of the resulting diagram
 */
void add(double fact1, char *name1, double fact2, char *name2, char *target)
{
    // A = a*A + b*B
    if (strcmp(name1, target) == 0) {
        dg_stack_pos_t pos = get_stack_pos();
        copy(name1, "buf");
        clear(name1);
        update(name1, fact1, "buf");
        update(name1, fact2, name2);
        restore_stack_pos(pos);
        return;
    }

    // B = a*A + b*B
    if (strcmp(name2, target) == 0) {
        dg_stack_pos_t pos = get_stack_pos();
        copy(name2, "buf");
        clear(name2);
        update(name2, fact1, name1);
        update(name2, fact2, "buf");
        restore_stack_pos(pos);
        return;
    }

    // C = a*A + b*B
    copy(name1, target);
    clear(target);
    update(target, fact1, name1);
    update(target, fact2, name2);
}


/**
 * Updates diagram dg1_name:
 * dg1 += factor * dg2
 */
void update(char *dg1_name, double factor, char *dg2_name)
{
    timer_new_entry("update", "Diagram addition (update)");
    timer_start("update");

    assert_diagram_exists(dg1_name);
    assert_diagram_exists(dg2_name);

    diagram_t *dg1 = diagram_stack_find(dg1_name);
    diagram_t *dg2 = diagram_stack_find(dg2_name);

    // diagrams cannot be added if they represent operators belonging to different irreps
    if (dg1->symmetry != dg2->symmetry) {
        errquit("update(): operators to be added have different symmetries");
    }

    // ranks must coincide
    if (dg1->rank != dg2->rank) {
        errquit("update(): ranks must coincide (%s:%d != %s:%d)", dg1_name, dg1->rank, dg2_name, dg2->rank);
    }

    // 'valence' and 'hole-particle' strings must coincide
    if (intcmp(dg1->rank, dg1->valence, dg2->valence) != 0) {
        summary(dg1_name);
        summary(dg2_name);
        errquit("update(): 'valence' strings must coincide");
    }
    if (intcmp(dg1->rank, dg1->qparts, dg2->qparts) != 0) {
        summary(dg1_name);
        summary(dg2_name);
        errquit("update(): 'qparts' strings must coincide");
    }
    if (intcmp(dg1->rank, dg1->t3space, dg2->t3space) != 0) {
        summary(dg1_name);
        summary(dg2_name);
        errquit("update(): 't3space' strings must coincide");
    }

    if (cc_opts->nthreads > 1) {
        omp_set_num_threads(cc_opts->nthreads);
    }

    for (size_t isb1 = 0; isb1 < dg1->n_blocks; isb1++) {
        block_t *block1 = dg1->blocks[isb1];
        if (block1->is_unique == 0) {
            continue;
        }

        block_t *block2 = diagram_get_block(dg2, block1->spinor_blocks);
        if (block1->is_unique && (block2->is_unique == 0)) {
            restore_block(dg2, block2);
        }

        block_load(block1);
        block_load(block2);

        if (block1->size != block2->size) {
            printf("block1->size = %ld\n", block1->size);
            printf("block2->size = %ld\n", block2->size);
            for (int i = 0; i < block1->rank; i++) {
                printf("%d ", block1->spinor_blocks[i]);
            }
            printf("\n");
            for (int i = 0; i < block2->rank; i++) {
                printf("%d ", block2->spinor_blocks[i]);
            }
            printf("\n");
            summary(dg1_name);
            summary(dg2_name);
            errquit("update(): size mismatch");
        }

        // internal threading is redundant here and typically slows down calculations
        xaxpy(WORKING_TYPE, block1->size, factor, block2->buf, block1->buf);

        block_unload(block2);
        block_store(block1);

        if (block1->is_unique && (block2->is_unique == 0)) {
            destroy_block(block2);
        }
    }

    timer_stop("update");
}

