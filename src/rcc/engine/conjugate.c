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
 */

#include "engine.h"

#include "datamodel.h"
#include "error.h"


/**
 * Constructs diagram which is Hermitian conjugate to the given one.
 */
void diagram_conjugate(char *source_name, char *target_name)
{
    diagram_t *src_diagram = diagram_stack_find(source_name);
    if (src_diagram == NULL) {
        errquit("diagram_conjugate(): diagram '%s' not found", source_name);
    }

    /*
     * interchange bra and ket parts. for example:
     * <ij|ab> -> <ab|ij>
     */
    int rank = src_diagram->rank;
    char order_ket_bra[CC_DIAGRAM_MAX_RANK];
    for (int i = 0; i < rank / 2; i++) {
        order_ket_bra[i] = '1' + i + rank / 2;
    }
    for (int i = rank / 2; i < rank; i++) {
        order_ket_bra[i] = '1' + i - rank / 2;
    }
    order_ket_bra[rank] = '\0';

    reorder(source_name, target_name, order_ket_bra);

    /*
     * new diagram should have 'normal order': 12, 1234, 123456, etc
     */
    char new_order[CC_DIAGRAM_MAX_RANK];
    for (int i = 0; i < rank; i++) {
        new_order[i] = '1' + i;
    }
    new_order[rank] = '\0';
    set_order(target_name, new_order);

    /*
     * complex conjugation of all matrix elements.
     * (obviously is not required for the case of real arithmetics)
     */
    if (arith == CC_ARITH_COMPLEX) {
        diagram_t *tgt_diagram = diagram_stack_find(target_name);
        for (size_t iblock = 0; iblock < tgt_diagram->n_blocks; iblock++) {

            block_t *block = tgt_diagram->blocks[iblock];

            if (block->storage_type == CC_DIAGRAM_DUMMY) {
                continue;
            }

            block_load(block);
            conj_vector(block->size, block->buf);
            block_store(block);
        }
    }
}
