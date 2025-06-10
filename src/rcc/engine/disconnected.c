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

#include "disconnected.h"

#include "engine.h"
#include "symmetry.h"

/**
 * Constructor of disconnected terms in lambda 0h0p equations.
 * <ab|ij> = <a|i> <b|j>
 */
void construct_disconnected_rank2_rank2(char *src1_name, char *src2_name, char *dst_name)
{
    /*
     * get source diagrams
     */
    diagram_t *src_diagram_1 = diagram_stack_find(src1_name);
    diagram_t *src_diagram_2 = diagram_stack_find(src2_name);

    if (src_diagram_1 == NULL) {
        errquit("construct_disconnected_rank2_rank2(): source diagram '%s' not found", src1_name);
    }
    if (src_diagram_2 == NULL) {
        errquit("construct_disconnected_rank2_rank2(): source diagram '%s' not found", src2_name);
    }

    /*
     * construct template of the destination diagram
     */
    char dst_qparts[] = "xxxx\0";
    dst_qparts[0] = src_diagram_1->qparts[0];
    dst_qparts[1] = src_diagram_2->qparts[0];
    dst_qparts[2] = src_diagram_1->qparts[1];
    dst_qparts[3] = src_diagram_2->qparts[1];

    char dst_valence[] = "xxxx\0";
    dst_valence[0] = src_diagram_1->valence[0] + '0';
    dst_valence[1] = src_diagram_2->valence[0] + '0';
    dst_valence[2] = src_diagram_1->valence[1] + '0';
    dst_valence[3] = src_diagram_2->valence[1] + '0';

    int symmetry_1 = src_diagram_1->symmetry;
    int symmetry_2 = src_diagram_2->symmetry;
    int prod_symmetry = mulrep2_abelian(symmetry_1, symmetry_2);

    tmplt_sym(dst_name, dst_qparts, dst_valence, "1234", NOT_PERM_UNIQUE, prod_symmetry);
    diagram_t *dst_diagram = diagram_stack_find(dst_name);

    /*
     * for each block:
     * for each matrix element:
     * <ab|ij> = <a|i> * <b|j>
     */
    int rank = 4;

    for (size_t iblock = 0; iblock < dst_diagram->n_blocks; iblock++) {
        block_t *block = dst_diagram->blocks[iblock];
        block_load(block);
        if (block->is_unique == 0) {
            continue;
        }

        int *indices = (int *) cc_malloc(block->size * rank * sizeof(int));
        block_gen_indices(block, indices);

        for (size_t i = 0; i < block->size; i++) {
            int idx_a = indices[4 * i + 0];
            int idx_b = indices[4 * i + 1];
            int idx_i = indices[4 * i + 2];
            int idx_j = indices[4 * i + 3];

            double complex t_ai = diagram_get_2(src_diagram_1, idx_a, idx_i);
            double complex t_bj = diagram_get_2(src_diagram_2, idx_b, idx_j);
            double complex t_abij = t_ai * t_bj;

            if (arith == CC_ARITH_COMPLEX) {
                block->buf[i] = t_abij;
            }
            else {
                ((double *) block->buf)[i] = creal(t_abij);
            }
        }

        cc_free(indices);

        block_store(block);
    }
}