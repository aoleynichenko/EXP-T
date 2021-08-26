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

/**
 * @file add.c
 * Subroutines for addition of Goldstone diagrams
 *
 * @author 2018-2021 Alexander Oleynichenko
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "engine.h"
#include "datamodel.h"
#include "linalg.h"
#include "options.h"
#include "utils.h"

#include "omp.h"


/*******************************************************************************
 * add
 *
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
 ******************************************************************************/
void add(double fact1, char *name1, double fact2, char *name2, char *target)
{
    timer_new_entry("add", "Diagram addition");
    timer_start("add");

    // A = a*A + b*B
    if (strcmp(name1, target) == 0) {
        dg_stack_pos_t pos = get_stack_pos();
        copy(name1, "buf");
        clear(name1);
        update(name1, fact1, "buf");
        update(name1, fact2, name2);
        restore_stack_pos(pos);
    }
        // B = a*A + b*B
    else if (strcmp(name2, target) == 0) {
        dg_stack_pos_t pos = get_stack_pos();
        copy(name2, "buf");
        clear(name2);
        update(name2, fact1, name1);
        update(name2, fact2, "buf");
        restore_stack_pos(pos);
    }
        // C = a*A + b*B
        // TODO: test this code
    else {
        copy(name1, target);
        clear(target);
        update(target, fact1, name1);
        update(target, fact2, name2);
    }

    timer_stop("add");
}


/*******************************************************************************
 * update
 *
 * Updates diagram dg1_name:
 * dg1 += factor * dg2
 ******************************************************************************/
void update(char *dg1_name, double factor, char *dg2_name)
{
    diagram_t *dg1, *dg2;
    double complex *buf1, *buf2;
    size_t size;
    block_t *block1, *block2;

    timer_new_entry("add-2", "Diagram addition (update)");
    timer_start("add-2");

    dg1 = diagram_stack_find(dg1_name);
    if (dg1 == NULL) {
        errquit("update(): diagram '%s' not found", dg1_name);
    }
    dg2 = diagram_stack_find(dg2_name);
    if (dg2 == NULL) {
        errquit("update(): diagram '%s' not found", dg2_name);
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
        block1 = dg1->blocks[isb1];

        // diagram_get_block возвращает почему-то неправильный указатель, в котором
        // первые два байта забиты единицами
        block2 = diagram_get_block(dg2, block1->spinor_blocks);

        if (block1->is_unique == 0) {
            continue;
        }

        // дабы не заполнять уникальный блок ерундой
        if (block1->is_unique && (block2->is_unique == 0)) {
            restore_block(dg2, block2);
        }

        symblock_load(block1);
        symblock_load(block2);

        buf1 = block1->buf;
        buf2 = block2->buf;
        size = block1->size;
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

            // enable internal threading
#if defined BLAS_MKL
        mkl_set_num_threads_local(cc_opts->nthreads);
#elif defined BLAS_OPENBLAS
        openblas_set_num_threads(cc_opts->nthreads);
#endif

        xaxpy(WORKING_TYPE, size, factor, buf2, buf1);

#if defined BLAS_MKL
        mkl_set_num_threads_local(1);
#elif defined BLAS_OPENBLAS
        openblas_set_num_threads(1);
#endif

        symblock_unload(block2);
        symblock_store(block1);

        if (block1->is_unique && (block2->is_unique == 0)) {
            destroy_block(block2);
        }
    }

    timer_stop("add-2");
}
