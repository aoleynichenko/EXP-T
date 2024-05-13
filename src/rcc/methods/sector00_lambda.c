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
 * Lambda equations in the 0h0p Fock space sector.
 *
 * Diagrams for lambda equations and density matrix can be found in:
 * I. Shavitt, R. J. Bartlett, Many-Body Methods in Chemistry and Physics,
 * Cambridge University Press (2009)
 *
 * For details on lambda equations in the relativistic framework, see:
 * A. Shee, L. Visscher, T. Saue.
 * Analytic one-electron properties at the 4-component relativistic
 * coupled cluster level with inclusion of spin-orbit coupling.
 * J. Chem. Phys. 145, 184107 (2016)
 * doi: 10.1063/1.4966643
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "ccutils.h"
#include "engine.h"
#include "diis.h"
#include "options.h"
#include "symmetry.h"
#include "utils.h"

void construct_lambda_singles_0h0p();

void construct_lambda_doubles_0h0p();

void construct_lambda_triples_0h0p();

void sector_0h0p_lambda_equations_guess();

void construct_disconnected_rank2_rank2(char *src1_name, char *src2_name, char *dst);

void sector_0h0p_analytic_density_matrix_lambda();

void sector_0h0p_calculate_natural_spinors();

void calc_L00_perturbative_triples_correction();

void diagram_clear_off_diagonal(char *diagram_name);

int diveps_lambda = 0;


/**
 * Iterative solution of lambda equations.
 * Finally density matrix is constructed and written to disk.
 */
int sector_0h0p_lambda_equations()
{
    int triples = triples_enabled();

    printf("\n");
    printf(" Lambda equations in the 0h0p sector\n");
    printf("\n");

    /*
     * initialize L1 and L2 amplitudes (diagrams L00_1c and L00_2c)
     */
    sector_0h0p_lambda_equations_guess();

    /*
     * for the CCSD(T) model additional one- and two-particle contributions
     * to lambda equations must be calculated and then added at each iteration
     */
    if (cc_opts->cc_model == CC_MODEL_CCSD_T3) {
        calc_L00_perturbative_triples_correction();
    }

    /*
     * main loop
     */
    diveps_lambda = 1;
    int success = solve_amplitude_equations(
            0, 0,
            "L00_1c", "L00_1nw",
            "L00_2c", "L00_2nw",
            NULL, NULL, NULL,
            construct_lambda_singles_0h0p,
            construct_lambda_doubles_0h0p,
            construct_lambda_triples_0h0p,
            NULL
    );
    diveps_lambda = 0;
    if (success == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    /*
     * analyze lambda amplitudes and save them to disk
     */
    print_cluster_operator_analysis(0, 0, "L00_1c", "L00_2c", NULL);
    save_cluster_amplitudes(0, 0, "L00_1c", "L00_2c", NULL, NULL);
    if (cc_opts->do_flush_amplitudes_txt) {
        save_cluster_amplitudes_formatted(0, 0, "L00_1c", "L00_2c", NULL);
    }

    return EXIT_SUCCESS;
}


/**
 * Constructs initial guess for Lambda amplitudes
 * (or read them from disk)
 */
void sector_0h0p_lambda_equations_guess()
{
    int calc_L1 = 1;
    int calc_L2 = 1;
    int calc_L3 = 1;

    if (cc_opts->reuse_amplitudes[0][0]) {
        printf(" Trying to read amplitudes from disk ...\n");
        if (diagram_read_binary("L00_1c.dg") != NULL) {
            printf(" L1 amplitudes successfully read from disk\n");
            calc_L1 = 0;
        }
        if (diagram_read_binary("L00_2c.dg") != NULL) {
            printf(" L2 amplitudes successfully read from disk\n");
            calc_L2 = 0;
        }
        printf("\n");
        /*if (triples) {
            if (diagram_read_binary("t3c.dg") != NULL) {
                printf(" T3 amplitudes successfully read from disk\n");
                calc_t3 = 0;
            }
            else {
                printf(" T3 amplitudes will be calculated\n");
            }
        }*/
    }

    if (calc_L1) {
        copy("ph", "L00_1c");
        diveps("L00_1c");
    }

    if (calc_L2) {
        copy("pphh", "L00_2c");
        diveps("L00_2c");
    }
}


/**
 * Lambda equations for Lambda-singles: lam_{a}^{i}
 */
void construct_lambda_singles_0h0p()
{
    tmplt("L00_1nw", "ph", "00", "12", NOT_PERM_UNIQUE);

    /*
     * in case of CCSD(T): additional contributions to lambda equations
     */
    if (cc_opts->cc_model == CC_MODEL_CCSD_T3) {
        update("L00_1nw", 1.0, "L00_1_t3corr");
    }

    // L00_S1
    update("L00_1nw", 1.0, "ph");

    dg_stack_pos_t pos = get_stack_pos();

    // L00_S2
    reorder("pphh", "r1", "1342");
    mult("r1", "t1c", "r2", 2);
    update("L00_1nw", 1.0, "r2");
    restore_stack_pos(pos);

    // L00_S3
    reorder("L00_1c", "r1", "21");
    mult("pp", "r1", "r2", 1);
    update("L00_1nw", 1.0, "r2");
    restore_stack_pos(pos);

    // L00_S4
    reorder("hh", "r1", "21");
    mult("L00_1c", "r1", "r2", 1);
    update("L00_1nw", -1.0, "r2");
    restore_stack_pos(pos);

    // L00_S5
    reorder("hpph", "r1", "2431");
    mult("r1", "L00_1c", "r2", 2);
    update("L00_1nw", 1.0, "r2");
    restore_stack_pos(pos);

    // L00_S6
    reorder("ph", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("L00_1c", "r2", "r3", 1);
    update("L00_1nw", -1.0, "r3");
    restore_stack_pos(pos);

    // L00_S7
    reorder("L00_1c", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("ph", "r2", "r3", 1);
    update("L00_1nw", -1.0, "r3");
    restore_stack_pos(pos);

    // L00_S8
    reorder("pphp", "r1", "2431");
    reorder("L00_1c", "r2", "21");
    mult("r1", "t1c", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update("L00_1nw", 1.0, "r4");
    restore_stack_pos(pos);

    // L00_S9
    reorder("hphh", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("L00_1c", "r2", "r3", 1);
    update("L00_1nw", -1.0, "r3");
    restore_stack_pos(pos);

    // L00_S10
    reorder("t1c", "r1", "21");
    reorder("pphp", "r2", "1324");
    mult("r1", "L00_1c", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update("L00_1nw", 1.0, "r4");
    restore_stack_pos(pos);

    // L00_S11
    reorder("L00_1c", "r1", "21");
    reorder("hphh", "r2", "2413");
    mult("r1", "t1c", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update("L00_1nw", -1.0, "r4");
    restore_stack_pos(pos);

    // L00_S12
    reorder("pphh", "r1", "2413");
    reorder("t2c", "r2", "4231");
    mult("r2", "L00_1c", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update("L00_1nw", 1.0, "r4");
    restore_stack_pos(pos);

    // L00_S13
    reorder("pphh", "r1", "3412");
    mult("r1", "t2c", "r2", 3);
    mult("L00_1c", "r2", "r3", 1);
    update("L00_1nw", -0.5, "r3");
    restore_stack_pos(pos);

    // L00_S14
    reorder("L00_1c", "r1", "21");
    reorder("t2c", "r2", "3412");
    mult("pphh", "r2", "r3", 3);
    mult("r3", "r1", "r4", 1);
    update("L00_1nw", -0.5, "r4");
    restore_stack_pos(pos);

    // L00_S15
    reorder("t1c", "r1", "21");
    reorder("pphh", "r2", "2431");
    mult("r1", "L00_1c", "r3", 1);
    mult("t1c", "r3", "r4", 1);
    mult("r2", "r4", "r5", 2);
    update("L00_1nw", -1.0, "r5");
    restore_stack_pos(pos);

    // L00_S16
    reorder("pphh", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("r2", "t1c", "r3", 1);
    mult("L00_1c", "r3", "r4", 1);
    update("L00_1nw", -1.0, "r4");
    restore_stack_pos(pos);

    // L00_S17
    reorder("pphh", "r1", "1342");
    reorder("t1c", "r2", "21");
    reorder("L00_1c", "r3", "21");
    mult("r1", "t1c", "r4", 2);
    mult("r4", "r2", "r5", 1);
    mult("r5", "r3", "r6", 1);
    update("L00_1nw", -1.0, "r6");
    restore_stack_pos(pos);

    // L00_S18
    reorder("L00_2c", "r1", "3412");
    mult("phpp", "r1", "r2", 3);
    update("L00_1nw", 0.5, "r2");
    restore_stack_pos(pos);

    // L00_S19
    reorder("hhhp", "r1", "3412");
    mult("L00_2c", "r1", "r2", 3);
    update("L00_1nw", -0.5, "r2");
    restore_stack_pos(pos);

    // L00_S20
    reorder("L00_2c", "r1", "3412");
    mult("pppp", "r1", "r2", 2);
    reorder("r2", "r3", "1342");
    mult("r3", "t1c", "r4", 2);
    update("L00_1nw", 0.5, "r4");
    restore_stack_pos(pos);

    // L00_S21
    reorder("t1c", "r1", "21");
    mult("hhhh", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("L00_2c", "r3", "r4", 3);
    update("L00_1nw", 0.5, "r4");
    restore_stack_pos(pos);

    // L00_S22
    reorder("hpph", "r1", "4132");
    reorder("L00_2c", "r2", "2314");
    mult("r1", "t1c", "r3", 1);
    mult("r2", "r3", "r4", 3);
    update("L00_1nw", -1.0, "r4");
    restore_stack_pos(pos);

    // L00_S23
    reorder("L00_2c", "r1", "4312");
    reorder("hpph", "r2", "2134");
    mult("r1", "t1c", "r3", 1);
    mult("r2", "r3", "r4", 3);
    update("L00_1nw", -1.0, "r4");
    restore_stack_pos(pos);

    // L00_S24
    reorder("ph", "r1", "21");
    reorder("L00_2c", "r2", "2341");
    mult("r1", "t2c", "r3", 1);
    mult("r2", "r3", "r4", 3);
    update("L00_1nw", -0.5, "r4");
    restore_stack_pos(pos);

    // L00_S25
    reorder("L00_2c", "r1", "4312");
    reorder("t2c", "r2", "2134");
    mult("r1", "r2", "r3", 3);
    mult("ph", "r3", "r4", 1);
    update("L00_1nw", -0.5, "r4");
    restore_stack_pos(pos);

    // L00_S26
    reorder("t2c", "r1", "3412");
    reorder("pphp", "r2", "1324");
    mult("r1", "L00_2c", "r3", 3);
    mult("r2", "r3", "r4", 2);
    update("L00_1nw", 0.5, "r4");
    restore_stack_pos(pos);

    // L00_S27
    reorder("L00_2c", "r1", "3412");
    reorder("phhh", "r3", "1324");
    mult("r1", "t2c", "r2", 3);
    mult("r3", "r2", "r4", 2);
    update("L00_1nw", -0.5, "r4");
    restore_stack_pos(pos);

    // L00_S28
    reorder("pphp", "r1", "2431");
    reorder("t2c", "r2", "1324");
    reorder("L00_2c", "r3", "4231");
    mult("r1", "r2", "r4", 2);
    mult("r4", "r3", "r5", 3);
    update("L00_1nw", 1.0, "r5");
    restore_stack_pos(pos);

    // L00_S29
    reorder("L00_2c", "r0", "2431");
    reorder("phhh", "r1", "4231");
    reorder("t2c", "r2", "1324");
    mult("r1", "r2", "r3", 2);
    mult("r0", "r3", "r4", 3);
    update("L00_1nw", -1.0, "r4");
    restore_stack_pos(pos);

    // L00_S30
    reorder("pphp", "r1", "3412");
    mult("r1", "t2c", "r2", 2);
    mult("L00_2c", "r2", "r3", 3);
    update("L00_1nw", -0.25, "r3");
    restore_stack_pos(pos);

    // L00_S31
    reorder("L00_2c", "r1", "3412");
    mult("r1", "t2c", "r2", 2);
    mult("phhh", "r2", "r3", 3);
    update("L00_1nw", 0.25, "r3");
    restore_stack_pos(pos);

    // L00_S32
    reorder("pphp", "r1", "2431");
    reorder("t1c", "r2", "21");
    reorder("L00_2c", "r3", "4312");
    mult("t1c", "r1", "r4", 1);
    mult("r4", "r2", "r5", 1);
    reorder("r5", "r6", "2143");
    mult("r6", "r3", "r7", 3);
    update("L00_1nw", -1.0, "r7");
    restore_stack_pos(pos);

    // L00_S33
    reorder("hphh", "r1", "3142");
    reorder("t1c", "r2", "21");
    mult("t1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    reorder("r4", "r5", "2431");
    mult("L00_2c", "r5", "r6", 3);
    update("L00_1nw", 1.0, "r6");
    restore_stack_pos(pos);

    // L00_S34
    reorder("pphp", "r1", "3412");
    mult("t1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    mult("L00_2c", "r4", "r5", 3);
    update("L00_1nw", -0.5, "r5");
    restore_stack_pos(pos);

    // L00_S35
    reorder("t1c", "r1", "21");
    reorder("t1c", "r2", "21");
    reorder("L00_2c", "r3", "3412");
    mult("r1", "phhh", "r4", 1);
    mult("r2", "r4", "r5", 1);
    reorder("r5", "r6", "3412");
    mult("r6", "r3", "r7", 3);
    update("L00_1nw", 0.5, "r7");
    restore_stack_pos(pos);

    // L00_S36
    reorder("pphh", "r1", "3142");
    reorder("L00_2c", "r2", "2341");
    reorder("t2c", "r3", "4123");
    mult("r1", "t1c", "r4", 2);
    mult("r2", "r3", "r5", 3);
    mult("r5", "r4", "r6", 1);
    update("L00_1nw", -0.5, "r6");
    restore_stack_pos(pos);

    // L00_S37
    reorder("t2c", "r1", "2134");
    reorder("L00_2c", "r2", "4312");
    reorder("pphh", "r3", "1342");
    mult("r2", "r1", "r4", 3);
    mult("r3", "t1c", "r5", 2);
    mult("r5", "r4", "r6", 1);
    update("L00_1nw", -0.5, "r6");
    restore_stack_pos(pos);

    // L00_S38
    reorder("pphh", "r1", "3421");
    reorder("t2c", "r2", "2413");
    reorder("L00_2c", "r3", "1342");
    mult("t1c", "r1", "r4", 1);
    reorder("r4", "r5", "2134");
    mult("r5", "r2", "r6", 2);
    mult("r3", "r6", "r7", 3);
    update("L00_1nw", -1.0, "r7");
    restore_stack_pos(pos);

    // L00_S39
    reorder("t1c", "r0", "21");
    reorder("pphh", "r1", "1342");
    reorder("t2c", "r2", "2413");
    reorder("L00_2c", "r3", "3142");
    mult("r1", "r2", "r4", 2);
    reorder("r4", "r5", "3412");
    mult("r5", "r0", "r6", 1);
    reorder("r6", "r7", "3412");
    mult("r7", "r3", "r8", 3);
    update("L00_1nw", -1.0, "r8");
    restore_stack_pos(pos);

    // L00_S40
    reorder("pphh", "r1", "3412");
    reorder("t1c", "r2", "21");
    mult("t2c", "r1", "r3", 2);
    mult("r3", "r2", "r4", 1);
    reorder("r4", "r5", "3412");
    mult("L00_2c", "r5", "r6", 3);
    update("L00_1nw", 0.25, "r6");
    restore_stack_pos(pos);

    // L00_S41
    reorder("L00_2c", "r0", "3412");
    reorder("pphh", "r1", "1342");
    reorder("t2c", "r2", "3412");
    mult("t1c", "r1", "r3", 1);
    mult("r2", "r3", "r4", 2);
    reorder("r4", "r5", "4312");
    mult("r5", "r0", "r6", 3);
    update("L00_1nw", 0.25, "r6");
    restore_stack_pos(pos);

    // L00_S42
    reorder("t2c", "r1", "3412");
    reorder("pphh", "r2", "1324");
    mult("r1", "L00_2c", "r3", 3);
    mult("r3", "t1c", "r4", 1);
    mult("r2", "r4", "r5", 2);
    update("L00_1nw", -0.5, "r5");
    restore_stack_pos(pos);

    // L00_S43
    reorder("L00_2c", "r1", "3412");
    reorder("t1c", "r2", "21");
    reorder("pphh", "r3", "1324");
    mult("t2c", "r1", "r4", 3);
    mult("r2", "r4", "r5", 1);
    mult("r3", "r5", "r6", 2);
    update("L00_1nw", -0.5, "r6");
    restore_stack_pos(pos);

    // L00_S44
    reorder("pphh", "r1", "3412");
    reorder("t1c", "r2", "21");
    mult("t1c", "r1", "r3", 1);
    mult("t1c", "r3", "r4", 1);
    mult("r2", "r4", "r5", 1);
    reorder("r5", "r6", "4123");
    mult("L00_2c", "r6", "r7", 3);
    update("L00_1nw", 0.5, "r7");
    restore_stack_pos(pos);

    // L00_S45
    reorder("pphh", "r1", "3412");
    reorder("L00_2c", "r2", "3412");
    mult("r1", "t1c", "r3", 1);
    reorder("r3", "r4", "3412");
    mult("t1c", "r2", "r5", 1);
    mult("t1c", "r5", "r6", 1);
    reorder("r6", "r7", "3412");
    mult("r4", "r7", "r8", 3);
    update("L00_1nw", 0.5, "r8");
    restore_stack_pos(pos);
}


/**
 * Lambda equations for Lambda-doubles: lam_{ab}^{ij}
 */
void construct_lambda_doubles_0h0p()
{
    tmplt("L00_2nw", "pphh", "0000", "1234", NOT_PERM_UNIQUE);

    /*
     * in case of CCSD(T): additional contributions to lambda equations
     */
    if (cc_opts->cc_model == CC_MODEL_CCSD_T3) {
        update("L00_2nw", 1.0, "L00_2_t3corr");
    }

    // L00_D1
    update("L00_2nw", 1.0, "pphh");

    dg_stack_pos_t pos = get_stack_pos();

    // L00_D2
    construct_disconnected_rank2_rank2("ph", "L00_1c", "r1");
    perm("r1", "(12|34)");
    update("L00_2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // L00_D3
    reorder("pphh", "r1", "1342");
    mult("r1", "t1c", "r2", 2);
    construct_disconnected_rank2_rank2("r2", "L00_1c", "r3");
    perm("r3", "(12|34)");
    update("L00_2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // L00_D4
    reorder("L00_1c", "r1", "21");
    mult("pphp", "r1", "r2", 1);
    perm("r2", "(34)");
    update("L00_2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // L00_D5
    reorder("hphh", "r1", "2341");
    mult("L00_1c", "r1", "r2", 1);
    perm("r2", "(12)");
    update("L00_2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // L00_D6
    reorder("L00_1c", "r1", "21");
    reorder("pphh", "r2", "1243");
    mult("r1", "t1c", "r3", 1);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "1243");
    perm("r5", "(34)");
    update("L00_2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // L00_D7
    reorder("t1c", "r1", "21");
    reorder("pphh", "r2", "2341");
    mult("L00_1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    perm("r4", "(12)");
    update("L00_2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // L00_D8
    reorder("L00_2c", "r1", "3412");
    mult("r1", "pp", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update("L00_2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // L00_D9
    reorder("hh", "r1", "21");
    mult("L00_2c", "r1", "r2", 1);
    perm("r2", "(34)");
    update("L00_2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // L00_D10
    reorder("L00_2c", "r1", "1342");
    reorder("hpph", "r2", "2413");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(12|34)");
    update("L00_2nw", 1.0, "r4");
    restore_stack_pos(pos);

    // L00_D11
    reorder("L00_2c", "r1", "3412");
    mult("pppp", "r1", "r2", 2);
    update("L00_2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // L00_D12
    reorder("hhhh", "r1", "3412");
    mult("L00_2c", "r1", "r2", 2);
    update("L00_2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // L00_D13
    reorder("ph", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("L00_2c", "r2", "r3", 1);
    perm("r3", "(34)");
    update("L00_2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // L00_D14
    reorder("t1c", "r1", "21");
    reorder("L00_2c", "r2", "3412");
    mult("ph", "r1", "r3", 1);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update("L00_2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // L00_D15
    reorder("pphp", "r1", "2431");
    reorder("L00_2c", "r2", "2341");
    mult("r1", "t1c", "r3", 2);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "4123");
    perm("r5", "(12)");
    update("L00_2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // L00_D16
    reorder("hphh", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("L00_2c", "r2", "r3", 1);
    perm("r3", "(34)");
    update("L00_2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // L00_D17
    reorder("pphp", "r1", "1342");
    reorder("L00_2c", "r2", "2413");
    mult("r1", "t1c", "r3", 1);
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12|34)");
    update("L00_2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // L00_D18
    reorder("L00_2c", "r1", "1342");
    reorder("hphh", "r2", "2413");
    reorder("t1c", "r3", "21");
    mult("r2", "r3", "r4", 1);
    mult("r1", "r4", "r5", 2);
    reorder("r5", "r6", "1324");
    perm("r6", "(12|34)");
    update("L00_2nw", -1.0, "r6");
    restore_stack_pos(pos);

    // L00_D19
    reorder("L00_2c", "r1", "3421");
    mult("r1", "t1c", "r2", 1);
    reorder("r2", "r3", "1243");
    mult("pphp", "r3", "r4", 2);
    update("L00_2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // L00_D20
    reorder("phhh", "r1", "2341");
    mult("t1c", "r1", "r3", 1);
    reorder("r3", "r4", "3412");
    mult("L00_2c", "r4", "r5", 2);
    update("L00_2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // L00_D21
    reorder("pphh", "r1", "3412");
    mult("r1", "t2c", "r2", 3);
    mult("L00_2c", "r2", "r3", 1);
    perm("r3", "(34)");
    update("L00_2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // L00_D22
    reorder("t2c", "r1", "3412");
    reorder("L00_2c", "r2", "3412");
    mult("pphh", "r1", "r3", 3);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update("L00_2nw", -0.5, "r5");
    restore_stack_pos(pos);

    // L00_D23
    reorder("t2c", "r1", "3412");
    reorder("pphh", "r2", "3412");
    mult("L00_2c", "r1", "r3", 3);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update("L00_2nw", -0.5, "r5");
    restore_stack_pos(pos);

    // L00_D24
    reorder("L00_2c", "r1", "3412");
    mult("r1", "t2c", "r2", 3);
    mult("pphh", "r2", "r3", 1);
    perm("r3", "(34)");
    update("L00_2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // L00_D25
    reorder("L00_2c", "r1", "1342");
    reorder("t2c", "r2", "1324");
    reorder("pphh", "r3", "2431");
    mult("r3", "r2", "r4", 2);
    mult("r1", "r4", "r5", 2);
    reorder("r5", "r6", "1324");
    perm("r6", "(12|34)");
    update("L00_2nw", 1.0, "r6");
    restore_stack_pos(pos);

    // L00_D26
    reorder("pphh", "r1", "3412");
    mult("r1", "t2c", "r2", 2);
    mult("L00_2c", "r2", "r3", 2);
    update("L00_2nw", 0.25, "r3");
    restore_stack_pos(pos);

    // L00_D27
    reorder("L00_2c", "r1", "3412");
    mult("r1", "t2c", "r2", 2);
    mult("pphh", "r2", "r3", 2);
    update("L00_2nw", 0.25, "r3");
    restore_stack_pos(pos);

    // L00_D28
    reorder("pphh", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("r2", "t1c", "r3", 1);
    mult("L00_2c", "r3", "r4", 1);
    perm("r4", "(34)");
    update("L00_2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // L00_D29
    reorder("L00_2c", "r1", "3412");
    reorder("t1c", "r2", "21");
    reorder("pphh", "r3", "1342");
    mult("r3", "t1c", "r4", 2);
    mult("r4", "r2", "r5", 1);
    mult("r1", "r5", "r6", 1);
    reorder("r6", "r7", "3412");
    perm("r7", "(12)");
    update("L00_2nw", -1.0, "r7");
    restore_stack_pos(pos);

    // L00_D30
    reorder("L00_2c", "r1", "1342");
    reorder("pphh", "r2", "2431");
    reorder("t1c", "r3", "21");
    mult("t1c", "r2", "r4", 1);
    mult("r4", "r3", "r5", 1);
    reorder("r5", "r6", "2314");
    mult("r1", "r6", "r7", 2);
    reorder("r7", "r8", "1324");
    perm("r8", "(12|34)");
    update("L00_2nw", -1.0, "r8");
    restore_stack_pos(pos);

    // L00_D31
    reorder("pphh", "r1", "3412");
    mult("t1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    mult("L00_2c", "r4", "r5", 2);
    update("L00_2nw", 0.5, "r5");
    restore_stack_pos(pos);

    // L00_D32
    reorder("L00_2c", "r1", "3412");
    mult("t1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    mult("pphh", "r4", "r5", 2);
    update("L00_2nw", 0.5, "r5");
    restore_stack_pos(pos);

    /*
     * contributions from the 0h1p sector if needed
     */
    if (cc_opts->sector_h == 0 && cc_opts->sector_p == 1) {

        // intermediate diagram - "cross"
        copy("dm1", "cross");

        /*reorder("L01_2c", "r1", "3412");
        mult("s2c", "r1", "r2", 3);
        update("cross", -0.5, "r2");
        restore_stack_pos(pos);*/

        // L00_DE1
        reorder("vphh", "r1", "2341");
        mult("cross", "r1", "r2", 1);
        tmplt("r3", "pphh", "0000", "1234", NOT_PERM_UNIQUE);
        expand_diagram("r2", "r3");
        perm("r3", "(12)");
        update("L00_2nw", -1.0, "r3");
        restore_stack_pos(pos);

        // L00_01_D9
        reorder("L01_2c", "r1", "1243");
        reorder("vh", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        perm("r4", "(34)");
        update("L00_2nw", -1.0, "r4");
        restore_stack_pos(pos);

        // L00_01_D10
        reorder("L01_2c", "r0", "2143"); // interchange electrons
        reorder("r0", "r1", "1342");
        reorder("vpph", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        perm("r4", "(12|34)");
        update("L00_2nw", 1.0, "r4");
        restore_stack_pos(pos);

        // L00_01_D12
        reorder("vhhh", "r1", "3412");
        mult("L01_2c", "r1", "r2", 2);
        update("L00_2nw", 0.5, "r2");
        restore_stack_pos(pos);

        // L00_D21
        reorder("pphh", "r1", "3412");
        reorder("L01_2c", "L01_2c_21", "2143"); // interchange electrons
        mult("r1", "s2c", "r2", 3);
        mult("L01_2c_21", "r2", "r3", 1);
        perm("r3", "(34)");
        update("L00_2nw", -0.5, "r3");
        restore_stack_pos(pos);

        // L00_D23
        reorder("s2c", "r1", "3412");
        reorder("pphh", "r2", "3412");
        mult("L01_2c", "r1", "r3", 3);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "3412");
        perm("r5", "(12)");
        update("L00_2nw", -0.5, "r5");
        restore_stack_pos(pos);

        // L00_D24
        reorder("L01_2c", "L01_2c_21", "2143"); // interchange electrons
        reorder("s2c", "s2c_21", "2143"); // interchange electrons
        reorder("L01_2c_21", "r1", "3412");
        mult("r1", "s2c_21", "r2", 3);
        mult("pphh", "r2", "r3", 1);
        perm("r3", "(34)");
        update("L00_2nw", -0.5, "r3");
        restore_stack_pos(pos);

        // L00_D25
        reorder("L01_2c", "L01_2c_21", "2143");
        reorder("L01_2c_21", "r1", "1342");
        reorder("s2c", "r2", "1324");
        reorder("pphh", "r3", "2431");
        mult("r3", "r2", "r4", 2);
        mult("r1", "r4", "r5", 2);
        reorder("r5", "r6", "1324");
        perm("r6", "(12|34)");
        update("L00_2nw", 1.0, "r6");
        restore_stack_pos(pos);

        // L00_D26
        reorder("pphh", "r1", "3412");
        mult("r1", "s2c", "r2", 2);
        mult("L01_2c", "r2", "r3", 2);
        update("L00_2nw", 0.25, "r3");
        restore_stack_pos(pos);
    }
}


/**
 * Lambda equations for Lambda-triples: lam_{abc}^{ijk}
 */
void construct_lambda_triples_0h0p()
{

}


/**
 * CCSD(T) perturbative triples contributions to lambda equations
 */
void calc_L00_perturbative_triples_correction()
{
    double time_start = abs_time();
    printf("\n calculate triples contribution to lambda equations\n");

    tmplt("L00_1_t3corr", "ph", "00", "12", NOT_PERM_UNIQUE);
    tmplt("L00_2_t3corr", "pphh", "0000", "1234", NOT_PERM_UNIQUE);

    diagram_conjugate("t3c", "t3c+");
    diagram_conjugate("t3(d)", "t3(d)+");

    dg_stack_pos_t pos = get_stack_pos();

    // contribution to singles

    reorder("t3c+", "r1", "145623");
    mult("r1", "hhpp", "r2", 4);
    update("L00_1_t3corr", 0.25, "r2");
    restore_stack_pos(pos);

    // contribution to doubles

    reorder("t3c+", "r1", "145623");
    mult("r1", "phpp", "r2", 3);
    reorder("r2", "r3", "1423");
    perm("r3", "(12)");
    update("L00_2_t3corr", 1.0, "r3");
    restore_stack_pos(pos);

    reorder("t3c+", "r1", "124356");
    reorder("hhhp", "r2", "3412");
    mult("r1", "r2", "r3", 3);
    perm("r3", "(34)");
    update("L00_2_t3corr", -1.0, "r3");
    restore_stack_pos(pos);

    reorder("t3(d)+", "r1", "145623");
    mult("r1", "phpp", "r2", 3);
    reorder("r2", "r3", "1423");
    perm("r3", "(12)");
    update("L00_2_t3corr", 0.5, "r3");
    restore_stack_pos(pos);

    reorder("t3(d)+", "r1", "124356");
    reorder("hhhp", "r2", "3412");
    mult("r1", "r2", "r3", 3);
    perm("r3", "(34)");
    update("L00_2_t3corr", -0.5, "r3");
    restore_stack_pos(pos);

    reorder("t3(d)+", "r1", "124563");
    mult("r1", "hp", "r2", 2);
    update("L00_2_t3corr", 1.0, "r2");
    restore_stack_pos(pos);

    printf(" done in %.1f seconds\n\n", abs_time() - time_start);
}


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

