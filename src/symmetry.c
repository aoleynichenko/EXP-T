/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2020 The EXP-T developers.
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
 * symmetry.c
 * ==========
 *
 * Tools for working with direct product tables (and other symmetry-related
 * things).
 *
 * 2018 Alexander Oleynichenko
 ******************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "engine.h"
#include "memory.h"
#include "symmetry.h"

// number of irreps
int nsym;

// is this group abelian or not
int is_abelian = 0;

// direct product table
int ***dir_prod_table;

int dir_prod_table_abelian[CC_MAX_NREP][CC_MAX_NREP];

// representation names
char **rep_names;

// number of fully symmetric irrep
int irrep_a1;

// point group type (Quaternion, Real, Complex)
int point_group_nz;

char point_group_name[64] = {'\0'};

void detect_point_group();

/*******************************************************************************
 * get_irrep_name
 *
 * returns pointer to the name of representation number 'irep'
 ******************************************************************************/
char *get_irrep_name(int irep)
{
    return rep_names[irep];
}


/*******************************************************************************
 * get_rep_number
 *
 * finds irrep named 'name'
 * returns -1 if this irrep not found
 ******************************************************************************/
int get_rep_number(char *name)
{
    for (size_t i = 0; i < nsym; i++) {
        if (strcmp(rep_names[i], name) == 0) {
            return i;
        }
    }
    return -1;
}


/*******************************************************************************
 * mulrep
 *
 * returns irreps which is a direct product of irreps 1 and 2.
 * this works only for abelian groups (all irreps are 1-dimensional)
 ******************************************************************************/
int mulrep2_abelian(int irep1, int irep2)
{
    return dir_prod_table[irep1][irep2][1];
}


int inverse_irrep_abelian(int irep)
{
    for (int irep2 = 0; irep2 < nsym; irep2++) {
        if (mulrep2_abelian(irep, irep2) == irrep_a1) {
            return irep2;
        }
    }
}


inline int get_totally_symmetric_irrep()
{
    return irrep_a1;
}


/*******************************************************************************
 * dpd_prod_to_sum
 *
 * direct product decomposition into direct sum.
 * arguments:
 *   gamma_inp[n_inp]   irreps to be multiplied
 *   gamma_out[n_out]   output -- direct sum
 * multipliers in the direct sum are omitted!
 * we use smth like "inverse polish calculator" to perform multiplications.
 * elements of the current sum are stored in the "stack" array.
 ******************************************************************************/
void dpd_prod_to_sum(int n_inp, int *gamma_inp, int *n_out, int *gamma_out)
{
    int op1, op2;
    int i, j, sum_size;
    int ptr_op;
    struct {
        int top;
        int elem[CC_MAX_DPD_STACK_SIZE];
    } stack, result;

    // is valid for dirac only:
    if (n_inp == 2) {
        op1 = gamma_inp[0];
        op2 = gamma_inp[1];
        gamma_out[0] = dir_prod_table[op1][op2][1];
        *n_out = 1;
        return;
    }

    ptr_op = n_inp - 1;    // last irrep in the input
    stack.top = 0;

    while (ptr_op >= 0) {
        op1 = gamma_inp[ptr_op--];   // get next irrep
        result.top = 0;

        if (stack.top == 0) {
            stack.elem[stack.top++] = op1;
        }
        else {
            for (i = 0; i < stack.top; i++) {
                op2 = stack.elem[i];
                sum_size = dir_prod_table[op1][op2][0];
                for (j = 1; j <= sum_size; j++) {
                    result.elem[result.top++] = dir_prod_table[op1][op2][j];
                }
            }
            stack.top = result.top;
            memcpy(stack.elem, result.elem, sizeof(int) * result.top);
        }
    }

    // return answer
    *n_out = stack.top;
    for (i = 0; i < stack.top; i++) {
        gamma_out[i] = stack.elem[i];
    }
}


/**
 * checks if the direct product of irreps [(<bra|)* x |ket>] contains the
 * totally symmetry irrep.
 * gamma = [bra, ket]   // single array
 *
 * @note implemented only for the case of abelian symmetry and diagrams
 * of rank <= 6 (triples).
 *
 * TODO: refactor name
 */

static int dpd_rank2(int *gamma)
{
    return gamma[0] == gamma[1];
}


static int dpd_rank4(int *gamma)
{
    int op1, op2;
    op1 = gamma[0];
    op2 = gamma[1];
    int left_rep = dir_prod_table_abelian[op1][op2];
    op1 = gamma[2];
    op2 = gamma[3];
    int right_rep = dir_prod_table_abelian[op1][op2];
    return left_rep == right_rep;
}


static int dpd_rank6(int *gamma)
{
    int op1, op2;

    op1 = gamma[0];
    op2 = gamma[1];
    op1 = dir_prod_table_abelian[op1][op2];
    op2 = gamma[2];
    int left_rep = dir_prod_table_abelian[op1][op2];

    op1 = gamma[3];
    op2 = gamma[4];
    op1 = dir_prod_table_abelian[op1][op2];
    op2 = gamma[5];
    int right_rep = dir_prod_table_abelian[op1][op2];
    return left_rep == right_rep;
}


// fast DPD is implemented only for models up to CCSDT (max rank == 6)
static int dpd_rank_other(int *gamma)
{
    int rank = 8;
    assert(rank <= 6);
}


int (*dpd_contains_totsym_rank[CC_MAX_NREP])(int *gamma) = {
        dpd_rank2, dpd_rank4, dpd_rank6, dpd_rank_other
};


/*******************************************************************************
 * dpd_contains_tot_symm_2
 *
 * checks if the direct product of irreps [(<bra|)* x |ket>] contains the
 * totally symmetry irrep.
 * we use smth like "inverse polish calculator" to perform multiplications.
 * elements of the current sum are stored in the "stack" array.
 *
 * @note it is a general implementation of the DPD scheme for the case of
 * abelian symmetry and diagrams of rank <= 6 (triples).
 *
 * TODO: refactor name
 ******************************************************************************/
int dpd_contains_tot_symm_2(int n_bra, int *gamma_bra, int n_ket, int *gamma_ket)
{
    int bra_sum[CC_MAX_DPD_STACK_SIZE];
    int ket_sum[CC_MAX_DPD_STACK_SIZE];
    int n_bra_sum, n_ket_sum;
    int i, j;

    dpd_prod_to_sum(n_bra, gamma_bra, &n_bra_sum, bra_sum);
    dpd_prod_to_sum(n_ket, gamma_ket, &n_ket_sum, ket_sum);

    for (i = 0; i < n_bra_sum; i++) {
        for (j = 0; j < n_ket_sum; j++) {
            if (bra_sum[i] == ket_sum[j]) {
                return 1; // contains totally symmetric rep
            }
        }
    }

    return 0; // no totally symmetric rep
}


void print_symmetry_info()
{
    int i, j, k;
    int prt_lvl;
    int repname_len = 0;
    // decorative strings
    char minus_field[64];   // '----'
    char equal_field[64];  // '===='
    char space_field[64];  // '    '

    prt_lvl = cc_opts->print_level;

    if (prt_lvl < CC_PRINT_HIGH) {
        return;
    }

    printf("\n");
    printf("\t\tSymmetry & irreducible representations\n");
    printf("\t\t--------------------------------------\n");
    printf("\n");
    printf(" Point group: %s\n", point_group_name);
    if (point_group_nz == 1) {
        printf(" Group type: real\n");
    }
    else if (point_group_nz == 2) {
        printf(" Group type: complex\n");
    }
    else if (point_group_nz == 4) {
        printf(" Group type: quaternion\n");
    }
    printf(" Arithmetic: %s\n", carith == 1 ? "complex" : "real");
    printf(" Is abelian: %d\n", is_abelian);
    printf(" Number of irreps: %d\n", nsym);
    printf(" Totally symmetric irrep: %s\n", get_irrep_name(get_totally_symmetric_irrep()));
    printf(" Representation names:\n");
    for (i = 0; i < nsym; i++) {
        if (strlen(rep_names[i]) > repname_len) {
            repname_len = strlen(rep_names[i]);
        }
    }
    for (i = 0; i < nsym; i++) {
        printf(" %2d %-*s%s", i, repname_len+1, rep_names[i], (i+1)%8 == 0 ? "\n" : "");
    }

    for (i = 0; i < repname_len; i++) {
        minus_field[i] = '-';
        equal_field[i] = '=';
        space_field[i] = ' ';
    }
    minus_field[repname_len] = '\0';
    equal_field[repname_len] = '\0';
    space_field[repname_len] = '\0';

    if (prt_lvl > CC_PRINT_HIGH) {
        printf("direct product table:\n");
        for (i = 0; i < nsym; i++) {
            for (j = 0; j < nsym; j++) {
                printf("%-*s  (x)  %-*s  =  ", repname_len, rep_names[i], repname_len, rep_names[j]);
                printf("%-*s", repname_len, rep_names[dir_prod_table[i][j][1]]);
                k = 2;
                while (dir_prod_table[i][j][k] != -1) {
                    printf("  (+)  %-*s", repname_len, rep_names[dir_prod_table[i][j][k]]);
                    k++;
                }
                printf("\n");
            }
        }
    }
    // finally print nice table
    if (is_abelian && nsym <= 32) {
        printf("\nmultiplication table (for abelian only):\n\n");
        // header
        printf("%s||", space_field);
        for (i = 0; i < nsym; i++) {
            printf("%-*s|", repname_len, rep_names[i]);
        }
        printf("\n");
        // separator
        printf("%s++", equal_field);
        for (i = 0; i < nsym; i++) {
            printf("%s+", equal_field);
        }
        printf("\n");
        for (i = 0; i < nsym; i++) {
            printf("%-*s||", repname_len, rep_names[i]);
            for (j = 0; j < nsym; j++) {
                printf("%-*s|", repname_len, rep_names[dir_prod_table[i][j][1]]);
            }
            printf("\n");
            // separator
            printf("%s++", minus_field);
            for (j = 0; j < nsym; j++) {
                printf("%s+", minus_field);
            }
            printf("\n");
        }
    }
    printf("\n");
}


void set_point_group_name(char *name)
{
    strcpy(point_group_name, name);
}


void symmetry_cleanup()
{
    size_t i, j;

    for (i = 0; i < nsym; i++) {
        cc_free(rep_names[i]);
    }
    cc_free(rep_names);
    rep_names = NULL;

    for (i = 0; i < nsym; i++) {
        for (j = 0; j < nsym; j++) {
            cc_free(dir_prod_table[i][j]);
        }
        cc_free(dir_prod_table[i]);
    }
    cc_free(dir_prod_table);
    dir_prod_table = NULL;
}
