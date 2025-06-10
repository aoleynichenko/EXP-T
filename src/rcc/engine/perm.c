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
 * Permutation operators acting on diagrams.
 *
 * Operators implemented:
 * P(ij), P(ab), P(ij|ab)
 * P(i/kj), P(ijk), P(a/bc), P(abc) + permutations + their products P(i/jk|abc) etc
 */

#include <ctype.h>
#include <string.h>

#include "engine.h"
#include "error.h"

int match_permutation_string(char *str, char *pattern);
void reorder_block(block_t *source_block, block_t *target_block, int *perm);
void elementary_perm(char *src_name, char *perm_str);
void safe_strncpy(char *dst, char *src, size_t n);

#define MAX_PERM_TASKS 64

void reverse_perm(int n, const int *direct_perm, int *inv_perm);


/**
 * Performs index permutation for the diagram.
 * P(ij) = 1 - Pij
 * P(ab) = 1 - Pab
 * P(ij|ab) = P(ij)P(ab) = 1 - Pij - Pab + Pij*Pab
 * P(i/jk) = 1 - Pij - Pik
 * P(i/jk|a/bc) = P(i/jk)P(a/bc)
 *
 * TODO: permutations for 3-particle diagrams
 * TODO: refactoring, pattern matching in a regexp-like manner
 */
void perm(char *dg_name, char *perm_str)
{
    int rk = rank(dg_name);

    dg_stack_pos_t pos = get_stack_pos();

    if (rk == 4) {
        if (strcmp(perm_str, "(12)") == 0) {
            elementary_perm(dg_name, "(12)");
        }
        else if (strcmp(perm_str, "(34)") == 0) {
            elementary_perm(dg_name, "(34)");
        }
        else if (strcmp(perm_str, "(12|34)") == 0) {
            elementary_perm(dg_name, "(12)");
            elementary_perm(dg_name, "(34)");
        }
    }
    else if (rk == 6) {
        if (match_permutation_string(perm_str, "(xx)")) {
            elementary_perm(dg_name, perm_str);
        }
        else if (match_permutation_string(perm_str, "(x/xx)")) {
            elementary_perm(dg_name, perm_str);
        }
        else if (strcmp(perm_str, "(456)") == 0) {
            elementary_perm(dg_name, perm_str);
        }
        else if (strcmp(perm_str, "(123)") == 0) {
            elementary_perm(dg_name, perm_str);
        }
        else if (match_permutation_string(perm_str, "(xx|yy)")) {
            char perm1[] = "(xx)";
            char perm2[] = "(yy)";
            perm1[1] = perm_str[1];
            perm1[2] = perm_str[2];
            perm2[1] = perm_str[4];
            perm2[2] = perm_str[5];
            elementary_perm(dg_name, perm1);
            elementary_perm(dg_name, perm2);
        }
        else if (match_permutation_string(perm_str, "(x/xx|y/yy)")) {
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
        else if (match_permutation_string(perm_str, "(xx|y/yy)")) {
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
        else if (match_permutation_string(perm_str, "(x/xx|yy)")) {
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
        else if (match_permutation_string(perm_str, "(x/xx|yyy)")) {
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
        else if (match_permutation_string(perm_str, "(xxx|y/yy)")) {
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
        else if (match_permutation_string(perm_str, "(xx|yyy)")) {
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
        else if (match_permutation_string(perm_str, "(xxx|yy)")) {
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
    }
    else {
        errquit("perm(): permutation operators for rank = %d are not yet implemented", rk);
    }

    restore_stack_pos(pos);
}


/*
 * works only with elementary permutations of types
 * (xx), (x/xx), (xxx)
 * 0123  012345  01234
 */
void elementary_perm(char *src_name, char *perm_str)
{
    typedef struct {
        char perm_str[CC_DIAGRAM_MAX_RANK];
        int sign;
    } perm_task_t;
    perm_task_t perm_tasks[MAX_PERM_TASKS];
    int n_perm_tasks;

    assert_diagram_exists(src_name);
    int rk = rank(src_name);

    if (rk == 4 && match_permutation_string(perm_str, "(xx)")) {
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
    else if (rk == 6 && match_permutation_string(perm_str, "(xx)")) {
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
    else if (rk == 6 && match_permutation_string(perm_str, "(x/xx)")) {
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
    else if (rk == 6 && strcmp(perm_str, "(456)") == 0) {
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
    else if (rk == 6 && strcmp(perm_str, "(123)") == 0) {
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

    timer_new_entry("permute", "Permutation operators");
    timer_start("permute");

    dg_stack_pos_t pos = get_stack_pos();
    copy(src_name, "_buf_dg");

    diagram_t *d_tgt = diagram_stack_find(src_name);
    diagram_t *d_src = diagram_stack_find("_buf_dg");

    int nthreads = 1;
    if (cc_opts->openmp_algorithm == CC_OPENMP_ALGORITHM_EXTERNAL) {
        nthreads = cc_opts->nthreads;
    }

    //#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
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

            // why -1 ??
            for (int i = 0; i < rk; i++) {
                inv_perm[i] -= 1;
            }

            int perm_spinor_blocks[CC_DIAGRAM_MAX_RANK];

            for (int i = 0; i < rk; i++) {
                perm_spinor_blocks[i] = b_tgt->spinor_blocks[perm[i]];
            }
            block_t *b_perm = diagram_get_block(d_src, perm_spinor_blocks);

            block_t *b_reordered = block_new(
                rk, b_tgt->spinor_blocks, d_src->qparts,
                d_src->valence, d_src->t3space,
                d_src->order, CC_DIAGRAM_IN_MEM, 0
            );

            reorder_block(b_perm, b_reordered, inv_perm);
            block_load(b_reordered);

            xaxpy(arith == CC_ARITH_COMPLEX ? CC_DOUBLE_COMPLEX : CC_DOUBLE,
                  b_tgt->size, sign, b_reordered->buf, b_tgt->buf);

            block_store(b_reordered);
            block_delete(b_reordered);
        }
        block_store(b_tgt);
    }
    restore_stack_pos(pos);

    timer_stop("permute");
}


/**
 * 'x' stands for a digit
 */
int match_permutation_string(char *str, char *pattern)
{
    if (strlen(str) == 0 || strlen(pattern) == 0) {
        return 0;
    }

    if (strlen(str) != strlen(pattern)) {
        return 0;
    }

    for (int i = 0; i < strlen(str); i++) {
        if (str[i] == '(' && pattern[i] == '(') {
            continue;
        }
        else if (str[i] == '/' && pattern[i] == '/') {
            continue;
        }
        else if (str[i] == '|' && pattern[i] == '|') {
            continue;
        }
        else if (str[i] == ')' && pattern[i] == ')') {
            continue;
        }
        else if (isdigit(str[i]) && (pattern[i] == 'x' || pattern[i] == 'y')) {
            continue;
        }
        else {
            return 0;
        }
    }

    return 1;
}


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



