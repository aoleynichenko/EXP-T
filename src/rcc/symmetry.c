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
 * Tools for working with direct product tables (and other symmetry-related
 * things).
 */

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
int **dir_prod_table_abelian;

// representation names
char **rep_names;

// number of fully symmetric irrep
int irrep_a1;

// point group type (Quaternion, Real, Complex)
int point_group_nz;

char point_group_name[64] = {'\0'};


void setup_symmetry(int group_type, char *group_name, int num_irreps, char **irrep_names, int fully_symmetric_rep, int *mult_table)
{
    point_group_nz = group_type;
    strcpy(point_group_name, group_name);
    nsym = num_irreps;
    rep_names = irrep_names;
    irrep_a1 = fully_symmetric_rep;

    /* general case */
    dir_prod_table = (int ***) cc_malloc(num_irreps * sizeof(int **));
    for (int i = 0; i < num_irreps; i++) {
        dir_prod_table[i] = (int **) cc_malloc(num_irreps * sizeof(int *));
        for (int j = 0; j < num_irreps; j++) {
            dir_prod_table[i][j] = (int *) cc_malloc((num_irreps + 1) * sizeof(int));
            for (int k = 0; k < num_irreps + 1; k++) {
                dir_prod_table[i][j][k] = -1;
            }
            // get values from the dirac_data structure (from MRCONEE files)
            dir_prod_table[i][j][0] = 1;
            dir_prod_table[i][j][1] = mult_table[i * num_irreps + j];
        }
    }

    /* special case: abelian groups */
    is_abelian = 1;
    dir_prod_table_abelian = (int **) cc_malloc(sizeof(int *) * num_irreps);
    for (int i = 0; i < num_irreps; i++) {
        dir_prod_table_abelian[i] = (int *) cc_malloc(sizeof(int) * num_irreps);
    }
    for (int i = 0; i < num_irreps; i++) {
        for (int j = 0; j < num_irreps; j++) {
            dir_prod_table_abelian[i][j] = mult_table[i * num_irreps + j];
        }
    }
}


/**
 * returns pointer to the name of representation number 'irep'
 */
char *get_irrep_name(int irep)
{
    return rep_names[irep];
}


/**
 * finds irrep named 'name'
 * returns -1 if this irrep not found
 */
int get_rep_number(char *name)
{
    for (size_t i = 0; i < nsym; i++) {
        if (strcmp(rep_names[i], name) == 0) {
            return i;
        }
    }
    return -1;
}


/**
 * returns irreps which is a direct product of irreps 1 and 2.
 * this works only for abelian groups (all irreps are 1-dimensional)
 */
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


int get_num_irreps()
{
    return nsym;
}


inline int get_totally_symmetric_irrep()
{
    return irrep_a1;
}


/**
 * direct product decomposition into direct sum.
 *
 * arguments:
 *   gamma_inp[n_inp]   irreps to be multiplied
 *   gamma_out[n_out]   output -- direct sum
 * multipliers in the direct sum are omitted!
 * we use smth like "inverse polish calculator" to perform multiplications.
 * elements of the current sum are stored in the "stack" array.
 */
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


int (*dpd_contains_totsym_rank[CC_MAX_NUM_IRREPS])(int *gamma) = {
        dpd_rank2, dpd_rank4, dpd_rank6, dpd_rank_other
};


static int dpd_op_sym_rank2(int *gamma, int op_sym)
{
    int g1 = gamma[0];
    int g2 = gamma[1];

    int op_g2 = dir_prod_table_abelian[op_sym][g2];

    return g1 == op_g2;
}


static int dpd_op_sym_rank4(int *gamma, int op_sym)
{
    int g1 = gamma[0];
    int g2 = gamma[1];
    int rep_bra = dir_prod_table_abelian[g1][g2];

    int g3 = gamma[2];
    int g4 = gamma[3];
    int rep_ket = dir_prod_table_abelian[g3][g4];

    int op_ket = dir_prod_table_abelian[op_sym][rep_ket];

    return rep_bra == op_ket;
}


static int dpd_op_sym_rank6(int *gamma, int op_sym)
{
    int g1 = gamma[0];
    int g2 = gamma[1];
    int g3 = gamma[2];
    int g4 = gamma[3];
    int g5 = gamma[4];
    int g6 = gamma[5];

    // irrep bra
    int g1_g2 = dir_prod_table_abelian[g1][g2];
    int rep_bra = dir_prod_table_abelian[g1_g2][g3];

    // irrep ket
    int g4_g5 = dir_prod_table_abelian[g4][g5];
    int rep_ket = dir_prod_table_abelian[g4_g5][g6];

    int op_ket = dir_prod_table_abelian[op_sym][rep_ket];

    return rep_bra == op_ket;
}


// fast DPD is implemented only for models up to CCSDT (max rank == 6)
static int dpd_op_sym_rank_other(int *gamma, int op_sym)
{
    int rank = 8;
    assert(rank <= 6);
}


int (*dpd_contains_totsym_rank_op_sym[CC_MAX_NUM_IRREPS])(int *gamma, int op_sym) = {
        dpd_op_sym_rank2, dpd_op_sym_rank4, dpd_op_sym_rank6, dpd_op_sym_rank_other
};


void print_symmetry_info()
{
    int i, j, k;
    int prt_lvl;
    int repname_len = 0;
    // decorative strings
    char minus_field[CC_MAX_NUM_IRREPS];   // '----'
    char equal_field[CC_MAX_NUM_IRREPS];  // '===='
    char space_field[CC_MAX_NUM_IRREPS];  // '    '

    prt_lvl = cc_opts->print_level;

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
    printf(" Arithmetic: %s\n", (arith == CC_ARITH_COMPLEX) ? "complex" : "real");
    printf(" Is abelian: %s\n", is_abelian ? "yes" : "no");
    printf(" Number of irreps: %d\n", nsym);
    printf(" Totally symmetric irrep: %s\n", get_irrep_name(get_totally_symmetric_irrep()));
    printf(" Representation names:\n");
    for (i = 0; i < nsym; i++) {
        if (strlen(rep_names[i]) > repname_len) {
            repname_len = strlen(rep_names[i]);
        }
    }
    for (i = 0; i < nsym; i++) {
        printf(" %3d %-*s%s", i, repname_len + 1, rep_names[i], (i + 1) % 8 == 0 ? "\n" : "");
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
                printf("[%d]\n", dir_prod_table[i][j][1]);
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
    if (is_abelian && nsym < 32) {
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
                int prod_rep = dir_prod_table_abelian[i][j];
                printf("%-*s|", repname_len, (prod_rep != -1) ? rep_names[prod_rep] : "--");
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


char *get_point_group_name()
{
    return point_group_name;
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


/**
 * returns irrep to which the vacuum (0h0p) electronic state belongs
 */
int get_vacuum_irrep()
{
    return get_totally_symmetric_irrep();
}


/**
 * Generates list of irreducible representations of the Cinfv group.
 * For example, list of irreps generated for Cinfv in DIRAC:
 *
  0 1/2+    1 1/2-    2 3/2+    3 3/2-    4 5/2+    5 5/2-    6 7/2+    7 7/2-
  8 9/2+    9 9/2-   10 11/2+  11 11/2-  12 13/2+  13 13/2-  14 15/2+  15 15/2-
 16 17/2+  17 17/2-  18 19/2+  19 19/2-  20 21/2+  21 21/2-  22 23/2+  23 23/2-
 24 25/2+  25 25/2-  26 27/2+  27 27/2-  28 29/2+  29 29/2-  30 31/2+  31 31/2-
 32 0      33 1+     34 1-     35 2+     36 2-     37 3+     38 3-     39 4+
 40 4-     41 5+     42 5-     43 6+     44 6-     45 7+     46 7-     47 8+
 48 8-     49 9+     50 9-     51 10+    52 10-    53 11+    54 11-    55 12+
 56 12-    57 13+    58 13-    59 14+    60 14-    61 15+    62 15-    63 16+
 */
char **generate_irreps_Cinfv(int max_omega_x2, int *n_irreps)
{
    *n_irreps = max_omega_x2 * 2 + 1;

    char **irrep_names = (char **) cc_malloc((*n_irreps) * sizeof(char *));
    for (int i = 0; i < (*n_irreps); i++) {
        irrep_names[i] = (char *) cc_malloc(MAX_IRREP_NAME * sizeof(char));
    }

    int j = 0;
    for (int i = 1; i < max_omega_x2; i += 2) {
        sprintf(irrep_names[j++], "%d/2+", i);
        sprintf(irrep_names[j++], "%d/2-", i);
    }
    sprintf(irrep_names[j++], "0");
    for (int i = 1; i <= max_omega_x2/2; i++) {
        sprintf(irrep_names[j++], "%d+", i);
        sprintf(irrep_names[j++], "%d-", i);
    }

    return irrep_names;
}


/*
 * Generates list of irreducible representations of the Dinfh group.
 * For example, list of irreps generated for Dinfh in DIRAC:
 *
  0 1/2g+    1 1/2g-    2 3/2g+    3 3/2g-    4 5/2g+    5 5/2g-    6 7/2g+    7 7/2g-
  8 9/2g+    9 9/2g-   10 11/2g+  11 11/2g-  12 13/2g+  13 13/2g-  14 15/2g+  15 15/2g-
 16 1/2u+   17 1/2u-   18 3/2u+   19 3/2u-   20 5/2u+   21 5/2u-   22 7/2u+   23 7/2u-
 24 9/2u+   25 9/2u-   26 11/2u+  27 11/2u-  28 13/2u+  29 13/2u-  30 15/2u+  31 15/2u-
 32 0g      33 1g+     34 1g-     35 2g+     36 2g-     37 3g+     38 3g-     39 4g+
 40 4g-     41 5g+     42 5g-     43 6g+     44 6g-     45 7g+     46 7g-     47 8g+
 48 0u      49 1u+     50 1u-     51 2u+     52 2u-     53 3u+     54 3u-     55 4u+
 56 4u-     57 5u+     58 5u-     59 6u+     60 6u-     61 7u+     62 7u-     63 8u+
 */
char **generate_irreps_Dinfh(int max_omega_x2, int *n_irreps)
{
    *n_irreps = max_omega_x2 * 4 + 2;

    char **irrep_names = (char **) cc_malloc((*n_irreps) * sizeof(char *));
    for (int i = 0; i < (*n_irreps); i++) {
        irrep_names[i] = (char *) cc_malloc(MAX_IRREP_NAME * sizeof(char));
    }

    int j = 0;

    // half-integer irreps, gerade
    for (int i = 1; i < max_omega_x2; i += 2) {
        sprintf(irrep_names[j++], "%d/2g+", i);
        sprintf(irrep_names[j++], "%d/2g-", i);
    }
    // half-integer irreps, ungerade
    for (int i = 1; i < max_omega_x2; i += 2) {
        sprintf(irrep_names[j++], "%d/2u+", i);
        sprintf(irrep_names[j++], "%d/2u-", i);
    }

    // integer irreps, gerade
    sprintf(irrep_names[j++], "0g");
    for (int i = 1; i <= max_omega_x2/2; i++) {
        sprintf(irrep_names[j++], "%dg+", i);
        sprintf(irrep_names[j++], "%dg-", i);
    }

    // integer irreps, ungerade
    sprintf(irrep_names[j++], "0u");
    for (int i = 1; i <= max_omega_x2/2; i++) {
        sprintf(irrep_names[j++], "%du+", i);
        sprintf(irrep_names[j++], "%du-", i);
    }

    return irrep_names;
}


void parse_infty_irrep_name(char *irrep_name, int *omega_x2, int *sign, int *gerade)
{
    char buf[MAX_IRREP_NAME];

    *sign = strchr(irrep_name, '-') ? -1 : 1;

    if (strchr(irrep_name, 'g')) {
        *gerade = +1;
    }
    else if (strchr(irrep_name, 'u')) {
        *gerade = -1;
    }
    else {
        *gerade = 0;
    }

    strcpy(buf, irrep_name);
    char *p_omega_x2 = strtok(buf, "/gu+-");
    *omega_x2 = atoi(p_omega_x2);
    int is_half_integer = strchr(irrep_name, '/') != NULL;
    if (!is_half_integer) {
        *omega_x2 = (*omega_x2) * 2;
    }
}


void multiply_irreps_Cinfv(char *irrep_1, char *irrep_2, char *prod_irrep)
{
    int omega_x2_1, sign_1, gerade_1;
    int omega_x2_2, sign_2, gerade_2;

    parse_infty_irrep_name(irrep_1, &omega_x2_1, &sign_1, &gerade_1);
    parse_infty_irrep_name(irrep_2, &omega_x2_2, &sign_2, &gerade_2);

    int prod_omega_x2 = sign_1 * omega_x2_1 + sign_2 * omega_x2_2;
    int prod_sign = (prod_omega_x2 < 0) ? -1 : +1;
    prod_omega_x2 = abs(prod_omega_x2);

    if (prod_omega_x2 == 0) { // zero irrep
        sprintf(prod_irrep, "0");
    }
    else if (prod_omega_x2 % 2 != 0) { // half-integer irrep
        sprintf(prod_irrep, "%d/2%s", prod_omega_x2, prod_sign == 1 ? "+" : "-");
    }
    else { // integer irrep
        sprintf(prod_irrep, "%d%s", prod_omega_x2 / 2, prod_sign == 1 ? "+" : "-");
    }
}


void multiply_irreps_Dinfh(char *irrep_1, char *irrep_2, char *prod_irrep)
{
    int omega_x2_1, sign_1, gerade_1;
    int omega_x2_2, sign_2, gerade_2;

    parse_infty_irrep_name(irrep_1, &omega_x2_1, &sign_1, &gerade_1);
    parse_infty_irrep_name(irrep_2, &omega_x2_2, &sign_2, &gerade_2);

    int prod_omega_x2 = sign_1 * omega_x2_1 + sign_2 * omega_x2_2;
    int prod_sign = (prod_omega_x2 < 0) ? -1 : +1;
    prod_omega_x2 = abs(prod_omega_x2);
    int prod_gerade = gerade_1 * gerade_2;

    if (prod_omega_x2 == 0) { // zero irrep
        sprintf(prod_irrep, "0%s", (prod_gerade > 0) ? "g" : "u");
    }
    else if (prod_omega_x2 % 2 != 0) { // half-integer irrep
        sprintf(prod_irrep, "%d/2%s%s", prod_omega_x2, (prod_gerade > 0) ? "g" : "u", prod_sign == 1 ? "+" : "-");
    }
    else { // integer irrep
        sprintf(prod_irrep, "%d%s%s", prod_omega_x2 / 2, (prod_gerade > 0) ? "g" : "u", prod_sign == 1 ? "+" : "-");
    }
}


int detect_operator_symmetry(char *bra_irrep_name, char *ket_irrep_name)
{
    int bra_irrep_num = get_rep_number(bra_irrep_name);
    int ket_irrep_num = get_rep_number(ket_irrep_name);

    // Cinfv/Dinfh or other groups?
    int rel_inf_group = 0;
    if (strcmp(point_group_name, "Cinfv") == 0 || strcmp(point_group_name, "Dinfh") == 0) {
        rel_inf_group = 1;
    }

    int min_omega_x2 = 1000;
    int irrep_with_min_omega_x2 = -1;

    for (int irep = 0; irep < get_num_irreps(); irep++) {


        if (mulrep2_abelian(bra_irrep_num, irep) == ket_irrep_num) {
        //if (mulrep2_abelian(irep, ket_irrep_num) == bra_irrep_num) {

            if (rel_inf_group) {
                int omega_x2, sign, gerade;
                parse_infty_irrep_name(get_irrep_name(irep), &omega_x2, &sign, &gerade);
                if (omega_x2 < min_omega_x2) {
                    min_omega_x2 = omega_x2;
                    irrep_with_min_omega_x2 = irep;
                }
            }
            else {
                irrep_with_min_omega_x2 = irep;
                break;
            }
        }
    }

    return irrep_with_min_omega_x2;
}


int search_string(char *str, char **str_list, int list_len)
{
    for (int i = 0; i < list_len; i++) {
        if (strcmp(str_list[i], str) == 0) {
            return i;
        }
    }

    return 0;
}


/*
 * Constructs irrep multiplication table.
 * The table will be stored in a linear array 'mul_table'.
 *
 * mult_table[i*n_irreps+j] = prod_function(irrep_i, irrep_j)
 *
 * returns pointer to the multiplication table
 */
int *construct_direct_product_table(int n_irreps, char **irrep_names, void (*prod_function)(char *a, char *b, char *prod))
{
    char prod_irrep[MAX_IRREP_NAME];

    int *mult_table = cc_malloc(sizeof(int) * n_irreps * n_irreps);

    for (int i = 0; i < n_irreps; i++) {
        for (int j = 0; j < n_irreps; j++) {
            char *irrep_1 = irrep_names[i];
            char *irrep_2 = irrep_names[j];

            prod_function(irrep_1, irrep_2, prod_irrep);

            mult_table[i * n_irreps + j] = search_string(prod_irrep, irrep_names, n_irreps);
        }
    }

    return mult_table;
}
