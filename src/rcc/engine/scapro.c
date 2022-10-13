/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2022 The EXP-T developers.
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
 *
 *  2018-2021 Alexander Oleynichenko
 */

#include "engine.h"

#include <complex.h>

#include "datamodel.h"
#include "error.h"
#include "linalg.h"


/**
 * Performs (elementwise) multiplication of two diagrams with subsequent
 * summation. Operands can be complex-conjugated.
 * Is used, for example, to calculate MP2 or CCSD correlation energy.
 * (for evaluation of closed diagrams)
 *
 * Arguments:
 *   conj1    specifies the form of the 'name1' diagram to be used in the scalar
 *            product:
 *            "N" -- no complex conjugation
 *            "C" -- perform complex conjugation
 *   conj2    the same for the second operand
 *   name1    symbolic name of the diagram -- first term
 *   name2    symbolic name of the diagram -- second term
 */
double complex scalar_product(char *conj1, char *conj2, char *name1, char *name2)
{
    diagram_t *dg1 = diagram_stack_find(name1);
    if (dg1 == NULL) {
        errquit("scapro(): diagram-1 '%s' not found", name1);
    }
    diagram_t *dg2 = diagram_stack_find(name2);
    if (dg2 == NULL) {
        errquit("scapro(): diagram-2 '%s' not found", name2);
    }

    if (dg1->symmetry != dg2->symmetry) {
        return 0.0 + 0.0 * I;
    }

    //printf("\nSCAPRO\n\n");

    double complex scal_prod = 0.0 + 0.0 * I;
    for (size_t isb1 = 0; isb1 < dg1->n_blocks; isb1++) {

        //printf(" isb = %d\n", isb1);

        block_t *block1 = dg1->blocks[isb1];
        block_t *block2 = diagram_get_block(dg2, block1->spinor_blocks);

        if (block2 == NULL) {
            // cannot calculate scalar product, diagrams are completely different
            errquit("scalar_product(): cannot calculate scalar product of diagrams '%s' and '%s'", name1, name2);
        }

        /*printf(" spinor blocks = ");
        for (int i = 0; i < block1->rank; i++) {
            printf("%d ", block1->spinor_blocks[i]);
        }
        printf("\n");

        printf(" unique-1 = %d     unique-2 = %d\n", block1->is_unique, block2->is_unique);
        printf(" sign = %d\n", block1->sign);*/

        /*
        if (block1->is_unique == 0 && block2->is_unique == 1) {
            restore_block(dg1, block1);
        }
        else if (block1->is_unique == 1 && block2->is_unique == 0) {
            restore_block(dg2, block2);
        }
        else if (block1->is_unique == 0 && block2->is_unique == 0) {
            continue;
        }*/

        int is_unique_1 = block1->is_unique;
        int is_unique_2 = block2->is_unique;
        int the_same_block = (block1 == block2);

        if (is_unique_1 == 0) {
            restore_block(dg1, block1);
        }
        if (is_unique_2 == 0 && !the_same_block) {
            restore_block(dg2, block2);
        }

        block_load(block1);
        block_load(block2);

        double complex *buf1 = block1->buf;
        double complex *buf2 = block2->buf;
        size_t size = block1->size;


        //printf(" uniq-1=%d    uniq-2=%d    nperm-1=%d    nperm-2=%d\n", block1->is_unique, block2->is_unique, block1->n_equal_perms, block2->n_equal_perms);
        //printf(" n_eq_perms = %d\n", block1->n_equal_perms);

        scal_prod += /*block1->n_equal_perms * */ xdot(WORKING_TYPE, conj1, conj2, size, buf1, buf2);

        if (is_unique_1 == 0) {
            destroy_block(block1);
        }
        if (is_unique_2 == 0 && !the_same_block) {
            destroy_block(block2);
        }

        /*
        if (block1->is_unique == 0 && block2->is_unique == 1) {
            destroy_block(block1);
        }
        else if (block1->is_unique == 1 && block2->is_unique == 0) {
            destroy_block(block2);
        }
        */

        block_unload(block1);
        block_unload(block2);

        //printf("\n");
    }

    return scal_prod;
}
