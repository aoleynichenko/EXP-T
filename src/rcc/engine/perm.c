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

/*
 * Permutation operators acting on diagrams.
 *
 * 2018-2021 Alexander Oleynichenko
 */

#include <ctype.h>
#include <string.h>

#include "datamodel.h"
#include "engine.h"
#include "error.h"
#include "templates.h"

void reorder_block_double(block_t *, block_t *, int *);
void reorder_block_double_complex_t(block_t *, block_t *, int *);

#define MAX_PERM_TASKS 64

void reverse_perm(int n, int *a, int *ainv);


/**
 * simple hand-coded strncpy implementation.
 * to prevent usage of compiler intrinsics and buffer overflow bugs.
 *
 * @param dst destination buffer
 * @param src source buffer
 * @param n number of chars to be copied
 */
void safe_strncpy(char *dst, char *src, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        dst[i] = src[i];
    }
    dst[n] = '\0';
}


/*
 * works only with elementary permutations of types
 * (xx), (x/xx), (xxx)
 * 0123  012345  01234
 */
void elementary_perm(char *src_name, char *perm_str)
{
    diagram_t *d_src, *d_tgt;
    int rk;

    typedef struct {
        char perm_str[CC_DIAGRAM_MAX_RANK];
        int sign;
    } perm_task_t;
    perm_task_t perm_tasks[MAX_PERM_TASKS];
    int n_perm_tasks;

    // rank-4 && '(xx)'
    if (rank(src_name) == 4 &&
        perm_str[0] == '(' && isdigit(perm_str[1]) && isdigit(perm_str[2]) && perm_str[3] == ')') {
        char p[] = "1234\0";
        int i1 = perm_str[1] - '0' - 1;
        int i2 = perm_str[2] - '0' - 1;
        char t = p[i1];
        p[i1] = p[i2];
        p[i2] = t;
        n_perm_tasks = 1;
        perm_tasks[0].sign = -1;
        safe_strncpy(perm_tasks[0].perm_str, p, 4);
    }
        // rank-6 && '(xx)'
    else if (rank(src_name) == 6 &&
             perm_str[0] == '(' && isdigit(perm_str[1]) && isdigit(perm_str[2]) && perm_str[3] == ')') {
        char p[] = "123456\0";
        int i1 = perm_str[1] - '0' - 1;
        int i2 = perm_str[2] - '0' - 1;
        char t = p[i1];
        p[i1] = p[i2];
        p[i2] = t;
        n_perm_tasks = 1;
        perm_tasks[0].sign = -1;
        safe_strncpy(perm_tasks[0].perm_str, p, 6);
    }
        // '(x/xx)'
    else if (rank(src_name) == 6 &&
             perm_str[0] == '(' && isdigit(perm_str[1]) && perm_str[2] == '/' &&
             isdigit(perm_str[3]) && isdigit(perm_str[4]) && perm_str[5] == ')') {
        char p1[] = "123456\0";
        char p2[] = "123456\0";
        int i1 = perm_str[1] - '0' - 1;
        int i2 = perm_str[3] - '0' - 1;
        int i3 = perm_str[4] - '0' - 1;
        char t = p1[i1];
        p1[i1] = p1[i2];
        p1[i2] = t;
        t = p2[i1];
        p2[i1] = p2[i3];
        p2[i3] = t;
        n_perm_tasks = 2;
        perm_tasks[0].sign = -1;
        safe_strncpy(perm_tasks[0].perm_str, p1, 6);
        perm_tasks[1].sign = -1;
        safe_strncpy(perm_tasks[1].perm_str, p2, 6);
    }
        // '(xxx)'
    else if (rank(src_name) == 6 &&
             perm_str[0] == '(' && isdigit(perm_str[1]) && isdigit(perm_str[2]) &&
             isdigit(perm_str[3]) && perm_str[4] == ')') {
        if (strcmp(perm_str, "(456)") == 0) {
            n_perm_tasks = 5;
            perm_tasks[0].sign = -1;
            safe_strncpy(perm_tasks[0].perm_str, "123465", 6);
            perm_tasks[1].sign = -1;
            safe_strncpy(perm_tasks[1].perm_str, "123546", 6);
            perm_tasks[2].sign = 1;
            safe_strncpy(perm_tasks[2].perm_str, "123564", 6);
            perm_tasks[3].sign = 1;
            safe_strncpy(perm_tasks[3].perm_str, "123645", 6);
            perm_tasks[4].sign = -1;
            safe_strncpy(perm_tasks[4].perm_str, "123654", 6);
        }
        else if (strcmp(perm_str, "(123)") == 0) {
            n_perm_tasks = 5;
            perm_tasks[0].sign = -1;
            safe_strncpy(perm_tasks[0].perm_str, "213456", 6);
            perm_tasks[1].sign = 1;
            safe_strncpy(perm_tasks[1].perm_str, "231456", 6);
            perm_tasks[2].sign = -1;
            safe_strncpy(perm_tasks[2].perm_str, "132456", 6);
            perm_tasks[3].sign = 1;
            safe_strncpy(perm_tasks[3].perm_str, "312456", 6);
            perm_tasks[4].sign = -1;
            safe_strncpy(perm_tasks[4].perm_str, "321456", 6);
        }
        else {
            errquit("permute(): wrong permutation: %s", perm_str);
        }
    }
    else {
        errquit("permute(): wrong permutation: %s", perm_str);
    }

    timer_new_entry("permute", "Permutation operators");
    timer_start("permute");

    dg_stack_pos_t pos = get_stack_pos();
    copy(src_name, "_buf_dg");

    d_tgt = diagram_stack_find(src_name);
    if (d_tgt == NULL) {
        errquit("permute(): diagram '%s' not found", src_name);
    }
    d_src = diagram_stack_find("_buf_dg");
    rk = d_tgt->rank;

    for (size_t ib = 0; ib < d_tgt->n_blocks; ib++) {
        block_t *b_tgt = d_tgt->blocks[ib];
        block_load(b_tgt);

        for (int itask = 0; itask < n_perm_tasks; itask++) {

            int perm[CC_DIAGRAM_MAX_RANK];
            int inv_perm[CC_DIAGRAM_MAX_RANK];
            int sign = perm_tasks[itask].sign;

            for (int i = 0; i < rk; i++) {
                perm[i] = perm_tasks[itask].perm_str[i] - '0' - 1;
            }
            reverse_perm(rk, perm, inv_perm);
            for (int i = 0; i < rk; i++) {
                inv_perm[i] -= 1;
            }

            int perm_spinor_blocks[CC_DIAGRAM_MAX_RANK];

            block_t *b_direct = diagram_get_block(d_src, b_tgt->spinor_blocks);

            for (int i = 0; i < rk; i++) {
                perm_spinor_blocks[i] = b_tgt->spinor_blocks[perm[i]];
            }
            block_t *b_perm = diagram_get_block(d_src, perm_spinor_blocks);

            block_t *b_reordered = block_new(rk, b_direct->spinor_blocks, d_src->qparts, d_src->valence, d_src->t3space,
                                             d_src->order, CC_DIAGRAM_IN_MEM, 0);

            if (arith == CC_ARITH_COMPLEX) {
                TEMPLATE(reorder_block, double_complex_t)(b_perm, b_reordered, inv_perm);

                block_load(b_reordered);
                for (size_t i = 0; i < b_direct->size; i++) {
                    b_tgt->buf[i] += sign * b_reordered->buf[i];
                }
                block_store(b_reordered);
            }
            else {
                TEMPLATE(reorder_block, double)(b_perm, b_reordered, inv_perm);

                block_load(b_reordered);
                double *b_tgt_buf = (double *) b_tgt->buf;
                double *b_reordered_buf = (double *) b_reordered->buf;
                for (size_t i = 0; i < b_direct->size; i++) {
                    b_tgt_buf[i] += sign * b_reordered_buf[i];
                }
                block_store(b_reordered);
            }

            block_delete(b_reordered);
        }
        block_store(b_tgt);
    }
    restore_stack_pos(pos);

    timer_stop("permute");
}


/**
 * Performs index permutation for the diagram.
 * P(ij) = 1 - Pij
 * P(ab) = 1 - Pab
 * P(ij|ab) = P(ij)P(ab) = 1 - Pij - Pab + Pij*Pab
 * P(i/jk) = 1 - Pij - Pik
 * P(i/jk|a/bc) = P(i/jk)P(a/bc)
 *
 * TODO: permutations for 3-particle diagrams
 * TODO: refactoring, patternmatching in a regexp-like manner
 */
void perm(char *dg_name, char *perm_str)
{
    int rk = rank(dg_name);

    dg_stack_pos_t pos = get_stack_pos();
    if (rk == 4 && strcmp(perm_str, "(12)") == 0) {
        elementary_perm(dg_name, "(12)");
    }
    else if (rk == 4 && strcmp(perm_str, "(34)") == 0) {
        elementary_perm(dg_name, "(34)");
    }
    else if (rk == 4 && strcmp(perm_str, "(12|34)") == 0) {
        elementary_perm(dg_name, "(12)");
        elementary_perm(dg_name, "(34)");
    }
        // case of (xx) permutations
        // note: for 6-rank diagram only!
    else if (rank(dg_name) == 6 &&
             perm_str[0] == '(' && perm_str[3] == ')' && isdigit(perm_str[1]) && isdigit(perm_str[2])) {
        elementary_perm(dg_name, perm_str);
    }
        // case of (x/xx) permutations
        // note: for 6-rank diagram only!
    else if (strlen(perm_str) == 6 &&
             perm_str[0] == '(' && perm_str[5] == ')' && perm_str[2] == '/' &&
             isdigit(perm_str[1]) && isdigit(perm_str[3]) && isdigit(perm_str[4])) {
        elementary_perm(dg_name, perm_str);
    }
        // case of (xxx) permutations
    else if (strcmp(perm_str, "(456)") == 0) {
        elementary_perm(dg_name, perm_str);
    }
    else if (strcmp(perm_str, "(123)") == 0) {
        elementary_perm(dg_name, perm_str);
    }
        // case of (xx|yy) permutations
    else if (strlen(perm_str) == 7 && perm_str[0] == '(' && perm_str[6] == ')' && perm_str[3] == '|') {
        char perm1[] = "(xx)";
        char perm2[] = "(yy)";
        perm1[1] = perm_str[1];
        perm1[2] = perm_str[2];
        perm2[1] = perm_str[4];
        perm2[2] = perm_str[5];
        elementary_perm(dg_name, perm1);
        elementary_perm(dg_name, perm2);
    }
        // case of (x/xx|y/yy) permutations
    else if (strlen(perm_str) == 11 && perm_str[0] == '(' && perm_str[10] == ')' && perm_str[5] == '|') {
        char perm1[] = "(x/xx)";
        char perm2[] = "(y/yy)";
        perm1[1] = perm_str[1];
        perm1[3] = perm_str[3];
        perm1[4] = perm_str[4];
        perm2[1] = perm_str[6];
        perm2[3] = perm_str[8];
        perm2[4] = perm_str[9];
        elementary_perm(dg_name, perm1);
        elementary_perm(dg_name, perm2);
    }
        // case of (xx|y/yy) permutations
    else if (strlen(perm_str) == 9 && perm_str[0] == '(' && perm_str[8] == ')' && perm_str[3] == '|') {
        char perm1[] = "(xx)";
        char perm2[] = "(y/yy)";
        perm1[1] = perm_str[1];
        perm1[2] = perm_str[2];
        perm2[1] = perm_str[4];
        perm2[3] = perm_str[6];
        perm2[4] = perm_str[7];
        elementary_perm(dg_name, perm1);
        elementary_perm(dg_name, perm2);
    }
        //         012345678
        // case of (x/xx|yy) permutations
    else if (strlen(perm_str) == 9 && perm_str[0] == '(' && perm_str[8] == ')' && perm_str[5] == '|') {
        char perm1[] = "(x/xx)";
        char perm2[] = "(yy)";
        perm1[1] = perm_str[1];
        perm1[3] = perm_str[3];
        perm1[4] = perm_str[4];
        perm2[1] = perm_str[6];
        perm2[2] = perm_str[7];
        elementary_perm(dg_name, perm1);
        elementary_perm(dg_name, perm2);
    }
        //         0123456789
        // case of (x/xx|yyy) permutations
    else if (strlen(perm_str) == 10 && perm_str[0] == '(' && perm_str[9] == ')' && perm_str[5] == '|' &&
             perm_str[2] == '/') {
        //              012345
        char perm1[] = "(x/xx)";
        //              01234
        char perm2[] = "(yyy)";
        perm1[1] = perm_str[1];
        perm1[3] = perm_str[3];
        perm1[4] = perm_str[4];
        perm2[1] = perm_str[6];
        perm2[2] = perm_str[7];
        perm2[3] = perm_str[8];
        elementary_perm(dg_name, perm1);
        elementary_perm(dg_name, perm2);
    }
        //         0123456789
        // case of (xxx|y/yy) permutations
    else if (strlen(perm_str) == 10 && perm_str[0] == '(' && perm_str[9] == ')' && perm_str[4] == '|' &&
             perm_str[6] == '/') {
        //              01234
        char perm1[] = "(xxx)";
        //              012345
        char perm2[] = "(y/yy)";
        perm1[1] = perm_str[1];
        perm1[2] = perm_str[2];
        perm1[3] = perm_str[3];
        perm2[1] = perm_str[5];
        perm2[3] = perm_str[7];
        perm2[4] = perm_str[8];
        elementary_perm(dg_name, perm1);
        elementary_perm(dg_name, perm2);
    }
        //         01234567
        // case of (xx|yyy) permutations
    else if (strlen(perm_str) == 8 && perm_str[0] == '(' && perm_str[7] == ')' && perm_str[3] == '|') {
        //              0123
        char perm1[] = "(xx)";
        //              01234
        char perm2[] = "(yyy)";
        perm1[1] = perm_str[1];
        perm1[2] = perm_str[2];
        perm2[1] = perm_str[4];
        perm2[2] = perm_str[5];
        perm2[3] = perm_str[6];
        elementary_perm(dg_name, perm1);
        elementary_perm(dg_name, perm2);
    }
        //         01234567
        // case of (xxx|yy) permutations
    else if (strlen(perm_str) == 8 && perm_str[0] == '(' && perm_str[7] == ')' && perm_str[4] == '|') {
        //              0123
        char perm1[] = "(xxx)";
        //              01234
        char perm2[] = "(yy)";
        perm1[1] = perm_str[1];
        perm1[2] = perm_str[2];
        perm1[3] = perm_str[3];
        perm2[1] = perm_str[5];
        perm2[2] = perm_str[6];
        elementary_perm(dg_name, perm1);
        elementary_perm(dg_name, perm2);
    }
    else {
        errquit("perm(): wrong permutation string %s", perm_str);
    }
    restore_stack_pos(pos);
}
