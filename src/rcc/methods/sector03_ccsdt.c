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
 * Multireference Fock-space coupled cluster method for model spaces
 * with three particles (Fock space sector 0h3p)
 *
 * Symbolic names used for diagrams:
 * z3c   current approximation to T^{0h3p}_3 amplitudes
 * z3nw  next (new) approximation to T^{0h3p}_3 amplitudes
 *
 * The CCSD model does not require solution of any amplitude equations.
 * Heff is constructed in one step from molecular integrals and Singles and
 * Doubles amplitudes from lower sectors:
 *   Heff(0h3p,CCSD) = f(T1,T2)
 *
 * This code was obtained by modification of the CC(0,0) program:
 * all we need is to flip hole creation lines to turn them to valence particle
 * annihilation lines. See Kaldor, J. Comp. Chem. V. 8, P.448 (1987) for details.
 * For the detailed description of the CCSD(0h,3p) model see:
 *   Hughes, Kaldor, Chem Phys Lett, V194, P99 (1992)
 *   Hughes, Kaldor, Chem Phys Lett, V204, P339 (1993)
 *
 * For details, see also:
 * L. V. Skripnikov, A. V. Oleynichenko, A. V. Zaitsevskii, D. E. Maison, A. E. Barzakh.
 * Relativistic Fock space coupled-cluster study of bismuth electronic structure
 * to extract the Bi nuclear quadrupole moment.
 * Phys. Rev. C, 104(3), 034316 (2021)
 * DOI: 10.1103/physrevc.104.034316
 *
 * 2019-2022 Alexander Oleynichenko
 */

#include "methods.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ccutils.h"
#include "engine.h"
#include "datamodel.h"
#include "diis.h"
#include "heff.h"
#include "options.h"
#include "sort.h"

void sort_integrals_0h3p();

void init_amplitudes_0h3p();

void const_terms_0h3p();

void construct_triples_0h3p();

void construct_folded_0h3p();

void heff03_ccsd();

void heff03_triples();

void restrict_valence(char *src_name /*large*/, char *tgt_name /*small*/, char *new_valence, int extract_valence);

void diag_T3_0h3p_c_ab();

void diag_T3_0h3p_a_bc();

void diag_T3_0h3p_abc();

void diag_T3_0h3p_k_ij_c_ab();

void diag_T3_0h3p_i_jk();

void diag_T3_0h3p_ijk();

void diag_T3_0h3p_k_ij();

void diag_T3_0h3p_k_ij_a_bc();

void diag_T3_0h3p_k_ij_abc();

void diag_T3_0h3p_i_jk_c_ab();

void diag_T3_0h3p_ijk_c_ab();


/**
 * Solves amplitude equations and find correlation energy for ground state
 * (Fock space sector 0h3p).
 */
int sector_0h3p(cc_options_t *opts)
{
    print_sector_banner(0, 3);
    opts->curr_sector_h = 0;
    opts->curr_sector_p = 3;

    sort_integrals_0h3p();

    /* CCSD, CCSDT-1 */
    if (cc_opts->cc_model < CC_MODEL_CCSDT_2) {
        heff03_ccsd();
        if (cc_opts->cc_model >= CC_MODEL_CCSD_T3) {
            heff03_triples();
        }
    }
        /* CCSDT-2, CCSDT-3, CCSDT */
    else {

        // setup shifts for the IH-like technique by A. V. Zaitsevskii
        if (intham_imms_in_sector(0, 3)) {
            intham_imms_setup(0, 3);
        }

        const_terms_0h3p();
        init_amplitudes_0h3p();
        intruder_state_analysis(NULL, NULL, "z3c");

        /*
         * main loop
         */
        int success = solve_amplitude_equations(
                0, 3,
                NULL, NULL,
                NULL, NULL,
                "z3c", "z3nw",
                "veff03",
                NULL,
                NULL,
                construct_triples_0h3p,
                construct_folded_0h3p
        );
        if (success == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }

        /*
         * analyze converged amplitudes and save them to disk
         */
        print_cluster_operator_analysis(0, 3, NULL, NULL, "z3c");
        save_cluster_amplitudes(0, 3, NULL, NULL, "z3c", "veff03");
        if (opts->do_flush_amplitudes_txt) {
            save_cluster_amplitudes_formatted(0, 3, NULL, NULL, "z3c");
        }
    }

    // следующий код загружает всю диаграмму в память для построения эфф гам-на:
    diagram_t *dg = diagram_stack_find("veff03");
    for (size_t isb = 0; isb < dg->n_blocks; isb++) {
        block_t *sb = dg->blocks[isb];
        block_load(sb);
        sb->storage_type = CC_DIAGRAM_IN_MEM;
    }

    /*
     * construct and diagonalize effective Hamiltonian matrix,
     * analyze its eigenvectors & eigenvalues
     */
    heff_analysis(0, 3, "veff01", "veff02", "veff03");

    /*
     * Ionization potentials
     */
    double ion_pot_1 = cc_opts->ground_energy_0h2p - cc_opts->ground_energy_0h3p;
    double ion_pot_2 = cc_opts->ground_energy_0h1p - cc_opts->ground_energy_0h2p;
    double ion_pot_3 = -cc_opts->ground_energy_0h1p;

    print_ionization_potential("Ionization potential 0h3p -> 0h2p", ion_pot_1);
    print_ionization_potential("Ionization potential 0h2p -> 0h1p", ion_pot_2);
    print_ionization_potential("Ionization potential 0h1p -> 0h0p", ion_pot_3);


    /*
     * calculate quasi-natural orbitals and model-space estimates of properties
     */
    calculate_model_space_properties(0, 3);
    calculate_model_space_natural_spinors(0, 3);

    return EXIT_SUCCESS;
}


/**
 * Sorting of integrals used in the CC(0h3p) amplitude equations.
 */
void sort_integrals_0h3p()
{
    request_sorting("vphp", "pphp", "1000", "1234");
    request_sorting("ppph", "ppph", "0000", "1234");
    request_sorting("vphh", "pphh", "1000", "1234");
    request_sorting("pvvv", "pppp", "0111", "1234");
    request_sorting("vvhv", "pphp", "1101", "1234");
    request_sorting("pphv", "pphp", "0001", "1234");
    request_sorting("ppvvr", "pppp", "0011", "3412");
    request_sorting("vphv", "pphp", "1001", "1234");
    request_sorting("ppvh", "ppph", "0010", "1234");
    request_sorting("vvvv", "pppp", "1111", "1234");
    request_sorting("pvph", "ppph", "0100", "1234");
    request_sorting("pvhv", "pphp", "0101", "1234");
    perform_sorting();

    if (cc_opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }
}


void const_terms_0h3p()
{
    double t1, t2;

    printf(" Construction of S^(0,3)-independent contributions to the FSCC-equations ...\n");
    timer_new_entry("const_s03", "Constant part of 0h3p amplitudes");
    timer_start("const_s03");

    tmplt("z3_0", "pppppp", "111000", "123456", IS_PERM_UNIQUE);

    t1 = abs_time();
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_2) {
        diag_T3_0h3p_k_ij_a_bc();
        diag_T3_0h3p_k_ij_abc();
        diag_T3_0h3p_i_jk_c_ab();
        diag_T3_0h3p_ijk_c_ab();
    }
    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        diag_T3_0h3p_k_ij_c_ab();
        diag_T3_0h3p_i_jk();
        diag_T3_0h3p_ijk();
        diag_T3_0h3p_k_ij();
    }
    t2 = abs_time();
    printf(" done in %.2f sec\n", t2 - t1);

    timer_stop("const_s03");
}


void init_amplitudes_0h3p()
{
    int calc_t3 = 1, calc_veff = 1;

    printf("\n Initialization of T{0h3p}_3 amplitudes ...\n");

    if (cc_opts->reuse_amplitudes[0][3]) {
        printf(" Trying to read amplitudes (sector 0h3p) from disk ...\n");
        if (diagram_read("z3c.dg") != NULL) {
            printf(" T{0h3p}_3 amplitudes successfully read from disk\n");
            calc_t3 = 0;
        }
        else {
            printf(" T{0h3p}_3 amplitudes will be calculated\n");
        }
    }

    if (calc_t3) {
        copy("z3_0", "z3c");
        closed("z3c", "veff03");
        diveps("z3c");
    }

    printf(" done\n\n");
}


// non-iterative construction of Heff at the CCSD level of theory
// diagrams are based on CCSDT-3 diagrams (lines were simply "reversed")
void heff03_ccsd()
{
    dg_stack_pos_t pos;

    timer_new_entry("03-CCSD-Heff", "0h3p -- CCSD part of eff Hamiltonian");
    timer_start("03-CCSD-Heff");

    printf("\n restrict valence ... ");
    restrict_valence("t1c", "t1-v", "01", 0);
    restrict_valence("t2c", "t2-v1", "0010", 0);
    restrict_valence("t2c", "t2-v2", "0001", 0);
    restrict_valence("t2c", "t2-v12", "0011", 0);
    restrict_valence("s2c", "s2-v1", "1010", 0);
    restrict_valence("s2c", "s2-v2", "1001", 0);
    restrict_valence("s2c", "s2-v12", "1011", 0);
    restrict_valence("x2c", "x2-v1", "1110", 0);
    restrict_valence("x2c", "x2-v2", "1101", 0);
    printf("done\n");

    printf(" T1 and T2 contributions to Heff(0h3p)\n");

    tmplt("veff03", "pppppp", "111111", "123456", NOT_PERM_UNIQUE);
    pos = get_stack_pos();

    // T1a
    printf(" [ 1/22] V * T{02}2\n");
    reorder("pvvv", "r1", "2341");
    mult("x2-v1", "r1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|4/56)");
    update("veff03", 1.0, "r3");
    restore_stack_pos(pos);

    // T1b
    printf(" [ 2/22] V * T{01}2\n");
    reorder("s2-v12", "r1", "3412");
    reorder("vvhv", "r2", "1243");
    mult("r1", "r2", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(1/23|6/45)");
    update("veff03", -1.0, "r4");
    restore_stack_pos(pos);

    // T3a
    printf(" [ 3/22] f * T{01}2 * T{02}2\n");
    reorder("ph", "phr_", "21");
    reorder("s2-v12", "s2_21", "2143");
    set_order("s2_21", "1234");
    reorder("s2_21", "r2", "2341");
    mult("x2-v1", "phr_", "r1", 1);
    mult("r1", "r2", "r3", 1);
    reorder("r3", "r4", "124356");
    perm("r4", "(3/12|4/56)");
    update("veff03", -1.0, "r4");
    restore_stack_pos(pos);

    // T3b
    printf(" [ 4/22] V * T{01}2 * T{02}2\n");
    reorder("pphv", "r1", "4231");
    reorder("s2-v1", "r2", "1324");
    reorder("x2-v2", "r3", "1243");
    mult("r2", "r1", "r4", 2);
    mult("r4", "r3", "r5", 1);
    reorder("r5", "r6", "145236");
    perm("r6", "(1/23|456)");
    update("veff03", 1.0, "r6");
    restore_stack_pos(pos);

    // T3c
    printf(" [ 5/22] V * T{01}2^2\n");
    reorder("s2-v1", "r1", "1324");
    reorder("pvhh", "r2", "2431");
    reorder("s2-v12", "s2_21", "2143");
    set_order("s2_21", "1234");
    reorder("s2_21", "r3", "2341");
    mult("r1", "r2", "r4", 2);
    mult("r4", "r3", "r5", 1);
    reorder("r5", "r6", "134256");
    perm("r6", "(123|4/56)");
    update("veff03", -1.0, "r6");
    restore_stack_pos(pos);

    // T3d
    // NOTE: wrong order of r5! 312645
    printf(" [ 6/22] V * T{01}2 * T{02}2\n");
    reorder("s2-v12", "r1", "3412");
    reorder("pphv", "r2", "4312");
    mult("x2c", "r2", "r3", 2);
    mult("r3", "r1", "r4", 1);
    reorder("r4", "r5", "612453");
    perm("r5", "(1/23|6/45)");
    update("veff03", -0.5, "r5");
    restore_stack_pos(pos);

    // T3e
    printf(" [ 7/22] V * T2 * T{02}2\n");
    reorder("t2-v12", "r1", "3412");
    reorder("pvhh", "r2", "2134");
    mult("r1", "r2", "r3", 2);
    mult("x2-v1", "r3", "r4", 1);
    reorder("r4", "r5", "126345");
    perm("r5", "(3/12|4/56)");
    update("veff03", 0.5, "r5");
    restore_stack_pos(pos);

    // T4a
    printf(" [ 8/22] V * T{01}1 * T{02}2\n");
    mult("s1c", "ppvvr", "r1", 1);
    mult("x2-v1", "r1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|4/56)");
    update("veff03", 1.0, "r3");
    restore_stack_pos(pos);

    // T4b
    printf(" [ 9/22] V * T1 * T{01}2\n");
    reorder("t1-v", "r2", "21");
    reorder("s2-v12", "r1", "1342");
    mult("r2", "vvhh", "r3", 1);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "156234");
    perm("r5", "(1/23|6/45)");
    update("veff03", 1.0, "r5");
    restore_stack_pos(pos);

    // T4c
    printf(" [10/22] V * T{01}1 * T{01}2\n");
    reorder("vphv", "r1", "1432");
    reorder("s2-v12", "r2", "1342");
    mult("s1c", "r1", "r3", 1);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "154236");
    perm("r5", "(123|6/45)");
    update("veff03", -1.0, "r5");
    restore_stack_pos(pos);

    // T4d
    printf(" [11/22] V * T1 * T{02}2\n");
    reorder("t1-v", "r1", "21");
    reorder("vphv", "r2", "1423");
    reorder("x2-v2", "r3", "1243");
    mult("r1", "r2", "r4", 1);
    mult("r4", "r3", "r5", 1);
    reorder("r5", "r6", "245136");
    perm("r6", "(1/23|456)");
    update("veff03", -1.0, "r6");
    restore_stack_pos(pos);

    // T7a
    printf(" [12/22] V * T1 * T{01}1 * T{02}2\n");
    reorder("pphv", "r1", "4231");
    reorder("t1-v", "r2", "21");
    reorder("x2-v2", "r3", "1243");
    mult("s1c", "r1", "r4", 1);
    mult("r2", "r4", "r5", 1);
    mult("r5", "r3", "r6", 1);
    reorder("r6", "r7", "245136");
    perm("r7", "(1/23|456)");
    update("veff03", -1.0, "r7");
    restore_stack_pos(pos);

    // T7b
    printf(" [13/22] V * T1 * T{01}1 * T{01}2\n");
    reorder("pvhh", "r1", "2431");
    reorder("t1-v", "r2", "21");
    reorder("s2-v12", "s2_21", "2143");
    set_order("s2_21", "1234");
    reorder("s2_21", "r3", "2341");
    mult("s1c", "r1", "r4", 1);
    mult("r2", "r4", "r5", 1);
    mult("r5", "r3", "r6", 1);
    reorder("r6", "r7", "234156");
    perm("r7", "(123|4/56)");
    update("veff03", 1.0, "r7");
    restore_stack_pos(pos);

    // T7c
    printf(" [14/22] V * T{01}1^2 * T{01}2\n");
    reorder("ppvh", "r1", "3412");
    reorder("s2-v12", "s2_21", "2143");
    set_order("s2_21", "1234");
    reorder("s2_21", "r2", "2341");
    mult("s1c", "r1", "r3", 1);
    mult("s1c", "r3", "r4", 1);
    mult("r4", "r2", "r5", 1);
    reorder("r5", "r6", "124356");
    perm("r6", "(3/12|4/56)");
    update("veff03", -1.0, "r6");
    restore_stack_pos(pos);

    // T7d
    printf(" [15/22] V * T1^2 * T{02}2\n");
    reorder("t1-v", "r1", "21");
    reorder("x2-v2", "r2", "1243");
    mult("r1", "vphh", "r3", 1);
    mult("r1", "r3", "r4", 1);
    mult("r4", "r2", "r5", 1);
    reorder("r5", "r6", "345126");
    perm("r6", "(1/23|6/45)");
    update("veff03", 1.0, "r6");
    restore_stack_pos(pos);

    // T8a
    // last reorder & perm like T7d -> combine ?
    printf(" [16/22] V * T1 * T{01}2 * T{02}2\n");
    reorder("pphh", "r1", "4231");
    reorder("s2-v12", "r2", "3412");
    reorder("x2-v2", "r3", "1243");
    mult("r1", "t1c", "r4", 2);
    mult("r3", "r4", "r5", 1);
    mult("r2", "r5", "r6", 1);
    reorder("r6", "r7", "345126");
    perm("r7", "(1/23|6/45)");
    update("veff03", -1.0, "r7");
    restore_stack_pos(pos);

    // T8b
    printf(" [17/22] V * T{01}1 * T{01}2^2\n");
    reorder("pphh", "r1", "4312");
    reorder("s2-v1", "r2", "1324");
    reorder("s2-v12", "s2_21", "2143");
    set_order("s2_21", "1234");
    reorder("s2_21", "r3", "2341");
    mult("s1c", "r1", "r4", 1);
    mult("r2", "r4", "r5", 2);
    mult("r5", "r3", "r6", 1);
    reorder("r6", "r7", "134256");
    perm("r7", "(123|4/56)");
    update("veff03", -1.0, "r7");
    restore_stack_pos(pos);

    // T8c
    printf(" [18/22] V * T1 * T{01}2 * T{02}2\n");
    reorder("pphh", "r1", "2314");
    reorder("t1-v", "r2", "21");
    reorder("s2-v1", "r3", "1324");
    reorder("x2-v2", "r4", "1243");
    mult("r2", "r1", "r5", 1);
    mult("r3", "r5", "r6", 2);
    mult("r6", "r4", "r7", 1);
    reorder("r7", "r8", "145236");
    perm("r8", "(1/23|456)");
    update("veff03", -1.0, "r8");
    restore_stack_pos(pos);

    // T8d
    printf(" [19/22] V * T2 * T{01}1 * T{02}2\n");
    reorder("pphh", "r1", "2134");
    reorder("t2-v12", "r2", "3412");
    reorder("x2-v2", "r3", "1243");
    mult("r2", "r1", "r4", 2);
    mult("s1c", "r4", "r5", 1);
    mult("r5", "r3", "r6", 1);
    reorder("r6", "r7", "145236");
    perm("r7", "(1/23|6/45)");
    update("veff03", 0.5, "r7");
    restore_stack_pos(pos);

    // T8e
    printf(" [20/22] V * T1 * T{01}2 * T{02}2\n");
    reorder("pphh", "r1", "4312");
    reorder("t1-v", "r2", "21");
    reorder("s2-v12", "s2_21", "2143");
    set_order("s2_21", "1234");
    reorder("s2_21", "r3", "2341");
    mult("x2c", "r1", "r4", 2);
    mult("r2", "r4", "r5", 1);
    mult("r5", "r3", "r6", 1);
    reorder("r6", "r7", "234156");
    perm("r7", "(3/12|4/56)");
    update("veff03", 0.5, "r7");
    restore_stack_pos(pos);

    // T10a
    printf(" [21/22] V * T1 * T{01}1^2 * T{02}1\n");
    reorder("pphh", "r1", "4123");
    reorder("t1-v", "r2", "21");
    reorder("s2-v12", "s2_21", "2143");
    set_order("s2_21", "1234");
    reorder("s2_21", "r3", "2341");
    mult("r2", "r1", "r4", 1);
    mult("s1c", "r4", "r5", 1);
    mult("s1c", "r5", "r6", 1);
    mult("r6", "r3", "r7", 1);
    reorder("r7", "r8", "124356");
    perm("r8", "(3/12|4/56)");
    update("veff03", 1.0, "r8");
    restore_stack_pos(pos);

    // T10b
    printf(" [22/22] V * T1^2 * T{01}1 * T{02}2\n");
    reorder("pphh", "r1", "2134");
    reorder("t1-v", "r2", "21");
    reorder("x2-v2", "r3", "1243");
    mult("r2", "r1", "r4", 1);
    mult("r2", "r4", "r5", 1);
    mult("s1c", "r5", "r6", 1);
    mult("r6", "r3", "r7", 1);
    reorder("r7", "r8", "145236");
    perm("r8", "(1/23|6/45)");
    update("veff03", 1.0, "r8");
    restore_stack_pos(pos);

    timer_stop("03-CCSD-Heff");
}


void heff03_triples()
{
    double t1, t2;

    // CCSDT-1: estimate S3, then estimate contribution to Heff

    printf(" +T(3) correction to the Heff(0h3p) operator:\n");

    timer_new_entry("03-S3", "0h3p -- Triples equations");
    timer_start("03-S3");
    t1 = abs_time();

    printf(" estimating T{0h3p}_3 triples amplitudes ...\n");
    tmplt("z3nw", "pppppp", "111000", "123456", IS_PERM_UNIQUE);
    dg_stack_pos_t pos = get_stack_pos();

    // T1a
    printf(" [ 1/ 8] V * T{02}2\n");
    reorder("pvpp", "r1", "2341");
    mult("x2c", "r1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|4/56)");
    update("z3nw", 1.0, "r3");
    restore_stack_pos(pos);

    // T1b
    printf(" [ 2/ 8] V * T{01}2\n");
    reorder("s2c", "r1", "3412");
    reorder("vvhp", "r2", "1243");
    mult("r1", "r2", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(1/23|6/45)");
    update("z3nw", -1.0, "r4");
    restore_stack_pos(pos);

    //print_ampl_vs_denom("z3nw", "z3c_eps.dat");

    if (cc_opts->cc_model == CC_MODEL_CCSDT_1B_PRIME) {
        // folded: - T{02}2 * V
        // undressed Veff: FS-CCSDT-1'
        // TODO: if 'IS_PERM_UNIQUE' in this diagram, the code is broken (segfault)
        printf(" [ 3/ 8] - T{02}2 * V\n");
        reorder("x2c", "r1", "2341");
        mult("vvvv", "r1", "r2", 1);
        reorder("r2", "r3", "124356");
        tmplt("r4", "pppppp", "111000", "123456", NOT_PERM_UNIQUE);
        expand_diagram("r3", "r4");
        perm("r4", "(3/12|4/56)");
        update("z3nw", -1.0, "r4");
        restore_stack_pos(pos);
    }
    else {
        // folded: - T{02}2 * Veff{02}
        // dressed Veff: FS-CCSD+T(3), FS-CCSD+T(4), FS-CCSDT-1
        // TODO: if 'IS_PERM_UNIQUE' in this diagram, the code is broken (segfault)
        printf(" [ 3/ 8] - T{02}2 * Veff{02}\n");
        reorder("x2c", "r1", "2341");
        mult("veff02", "r1", "r2", 1);
        reorder("r2", "r3", "124356");
        tmplt("r4", "pppppp", "111000", "123456", NOT_PERM_UNIQUE);
        expand_diagram("r3", "r4");
        perm("r4", "(3/12|4/56)");
        update("z3nw", -1.0, "r4");
        restore_stack_pos(pos);
    }

    closed("z3nw", "z3nw_cl");
    diveps("z3nw");

    timer_stop("03-S3");
    t2 = abs_time();
    printf(" estimated in %.3f seconds\n", t2 - t1);


    // Triples contribution to Heff(0h3p)
    printf(" Evalution of triples correction +T(3)\n");
    timer_new_entry("03-T(3)", "0h3p -- Triples contribution to Heff");
    timer_start("03-T(3)");
    t1 = abs_time();

    tmplt("+t(3)", "pppppp", "111111", "123456", NOT_PERM_UNIQUE);

    restrict_valence("z3nw", "z3nw-v12", "111110", 0);
    restrict_valence("z3nw", "z3nw-v1", "111100", 0);
    restrict_valence("x3c", "x3-v123", "110111", 0);
    restrict_valence("x3c", "x3-v12", "110110", 0);
    restrict_valence("s3c", "s3-v123", "100111", 0);

    dg_stack_pos_t pos2 = get_stack_pos();

    // T2a
    printf(" [ 4/ 8] f * T{03}3\n");
    reorder("x3-v123", "r1", "456123");
    mult("r1", "vh", "r2", 1);
    reorder("r2", "r3", "456123");
    perm("r3", "(3/12)");
    update("+t(3)", -1.0, "r3");
    restore_stack_pos(pos2);

    // T2b
    printf(" [ 5/ 8] f * T{02}3\n");
    reorder("pv", "ppr_", "21");
    mult("z3nw-v12", "ppr_", "r1", 1);
    perm("r1", "(6/45)");
    update("+t(3)", 1.0, "r1");
    restore_stack_pos(pos2);

    // T2c
    printf(" [ 6/ 8] V * T{03}3\n");
    mult("z3nw-v1", "ppvvr", "r1", 2);
    perm("r1", "(4/56)");
    update("+t(3)", 0.5, "r1");
    restore_stack_pos(pos2);

    // T2d
    printf(" [ 7/ 8] V * T{01}3\n");
    reorder("s3-v123", "r1", "456123");
    mult("r1", "vvhh", "r2", 2);
    reorder("r2", "r3", "456123");
    perm("r3", "(1/23)");
    update("+t(3)", 0.5, "r3");
    restore_stack_pos(pos2);

    // T2e
    printf(" [ 8/ 8] V * T{02}3\n");
    reorder("x3-v12", "r1", "124536");
    reorder("pvhv", "r2", "2431");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "125346");
    perm("r4", "(3/12|6/45)");
    update("+t(3)", 1.0, "r4");
    restore_stack_pos(pos2);

    // add +T(3) correction to the effective interaction operator
    update("veff03", 1.0, "+t(3)");

    timer_stop("03-T(3)");
    t2 = abs_time();
    printf(" +T(3) correction done in %.3f seconds\n", t2 - t1);
}


/**
 * Triples amplitude equations.
 * This code is invoked only if cc_model >= CCSDT
 */
void construct_triples_0h3p()
{
    timer_new_entry("03-T3", "0h3p -- Triples equations (T3)");
    timer_start("03-T3");

    //tmplt("z3nw", "pppppp", "111000", "123456");
    copy("z3_0", "z3nw");

    // constant terms
    //diag_T3_0h3p_k_ij_a_bc();
    //diag_T3_0h3p_k_ij_abc();
    //diag_T3_0h3p_i_jk_c_ab();
    //diag_T3_0h3p_ijk_c_ab();

    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        diag_T3_0h3p_c_ab();
        diag_T3_0h3p_a_bc();
        diag_T3_0h3p_abc();
    }

    // constant terms
    //diag_T3_0h3p_k_ij_c_ab();
    //diag_T3_0h3p_i_jk();
    //diag_T3_0h3p_ijk();
    //diag_T3_0h3p_k_ij();

    timer_stop("03-T3");
}


/*
 * Final T3 permutation: (1/23|6/45) aka (i/jk|c/ab)
 * Diagrams included:
 * (i/jk|c/ab): T1b T3a T3d T4b T5d T7c T8a T8e T10a
 */
void diag_T3_0h3p_i_jk_c_ab()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "ppph", "1100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("vvhp", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3a
        reorder("ph", "phr_", "21");
        reorder("x2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3c
        /*reorder("t2c", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        perm("r4", "(12)");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);*/

        // T3d
        reorder("t2c", "r1", "3412");
        reorder("pphp", "r2", "4312");
        mult("x2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4b
        reorder("t1c", "r2", "21");
        mult("vvhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        /*mult("t1c", "phhp", "r1", 1);
        perm("r1", "(12)");
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);*/

        // T7b
        /*reorder("t1c", "r1", "21");
        reorder("hphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        perm("r5", "(12)");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);*/

        // T7c
        reorder("pphp", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("x2c", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        /*reorder("pphh", "r1", "3142");
        reorder("t2c", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        perm("r4", "(12)");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);*/

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("t1c", "r3", "21");
        mult("x2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        // T10a
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("s1c", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    if (cc_opts->cc_model > CC_MODEL_CCSDT_3) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        reorder("s3c", "t5d_r1", "146235");
        reorder("pphh", "t5d_r2", "2341");
        reorder("x2c", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        diagram_stack_erase("t5d_r1");
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        diagram_stack_erase("t5d_r4");
        reorder("t5d_r5", "t5d_r6", "156234");
        diagram_stack_erase("t5d_r5");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(1/23|6/45)");
    update("z3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123|6/45) aka (ijk|c/ab)
 * Diagrams included:
 * (ijk|c/ab):  T3c T4c T7b T8b
 */
void diag_T3_0h3p_ijk_c_ab()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "ppph", "1100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3c
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        reorder("pvhp", "r0", "2431");
        mult("s1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("t1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(123|6/45)");
    update("z3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|4/56) aka (k/ij|a/bc)
 * Diagrams included:
 * (k/ij|a/bc): T1a T3e T4a T5e T7d T8d T10b
 */
void diag_T3_0h3p_k_ij_a_bc()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pppp", "1000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("pvpp", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b
        /*reorder("ppph", "r1", "1342");
        reorder("t2c", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        perm("r4", "(34)");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);*/

        // T3e
        reorder("t2c", "r1", "3412");
        mult("r1", "pvhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4a
        mult("ppppr", "s1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T4d
        /*reorder("phhp", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        perm("r3", "(34)");
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);*/

        // T7a
        /*reorder("t1c", "r1", "21");
        reorder("ppph", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        perm("r5", "(34)");
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);*/

        // T7d
        reorder("t1c", "r1", "21");
        mult("r1", "pvhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8c
        /*reorder("t1c", "r1", "21");
        reorder("t2c", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        perm("r6", "(34)");
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);*/

        // T8d
        reorder("t2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("s1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("t1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("x2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    if (cc_opts->cc_model > CC_MODEL_CCSDT_3) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        reorder("x3c", "t5e_r1", "124653");  // interchange 2 <-> 3, then reorder 134562
        reorder("pphh", "t5e_r2", "4123");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "t5e_r3", "2341");
        mult("t5e_r1", "t5e_r2", "t5e_r4", 3);
        diagram_stack_erase("t5e_r1");
        mult("t5e_r4", "t5e_r3", "t5e_r5", 1);
        diagram_stack_erase("t5e_r3");
        reorder("t5e_r5", "t5e_r6", "124356");
        diagram_stack_erase("t5e_r5");
        update("r3", -0.5, "t5e_r6");
        restore_stack_pos(pos2);
    }

    perm("r3", "(3/12|4/56)");
    update("z3_0", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|456) aka (k/ij|abc)
 * Diagrams included:
 * (k/ij|abc):  T3b T4d T7a T8c
 */
void diag_T3_0h3p_k_ij_abc()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pppp", "1000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b
        reorder("ppph", "r1", "1342");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4d
        reorder("pvhp", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("t1c", "r1", "21");
        reorder("ppph", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("s2c", "s2c_21", "2143");
        reorder("t1c", "r1", "21");
        reorder("s2c_21", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    mult("x2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|456)");
    update("z3_0", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (6/45) aka (c/ab)
 * Diagrams included:
 * T2a T5c T6b T6e T9b
 */
void diag_T3_0h3p_c_ab()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }

    tmplt("i1", "pp", "00", "21", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2a
        reorder("pp", "ppr_", "21");
        update("i1", 1.0, "ppr_");
        restore_stack_pos(pos2);

        // T5c
        reorder("t2c", "r1", "3412");
        mult("r1", "pphh", "r2", 3);
        update("i1", -0.5, "r2");
        restore_stack_pos(pos2);

        // T6b
        reorder("t1c", "r1", "21");
        mult("r1", "ph", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T6e
        reorder("ppph", "r1", "3142");
        mult("r1", "t1c", "r2", 2);
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T9b
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "1342");
        mult("r2", "t1c", "r3", 2);
        mult("r1", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    mult("z3c", "i1", "r1", 1);
    perm("r1", "(6/45)");
    update("z3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (4/56) aka (a/bc)
 * Diagrams included:
 * (a/bc): T2c T5g T9d
 */
void diag_T3_0h3p_a_bc()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }

    tmplt("i1", "pppp", "0000", "3412", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2c
        update("i1", 0.5, "ppppr");

        // T5g
        reorder("t2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        update("i1", 0.25, "r2");
        restore_stack_pos(pos2);

        // T9d
        reorder("t1c", "r1", "21");
        reorder("t1c", "r2", "21");
        mult("r1", "pphh", "r3", 1);
        mult("r2", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);

        // T6g
        /*reorder("t1c", "r1", "21");
        mult("ppph", "r1", "r2", 1);
        perm("r2", "(34)");
        reorder("r2", "r3", "3412");
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);*/
    }
    mult("z3c", "i1", "r1", 2);
    perm("r1", "(4/56)");
    update("z3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (456) aka (abc)
 * Diagrams included:
 * (abc):  T6g
 */
void diag_T3_0h3p_abc()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }

    tmplt("i1", "pppp", "0000", "3412", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T6g
        reorder("t1c", "r1", "21");
        mult("ppph", "r1", "r2", 1);
        reorder("r2", "r3", "3412");
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);
    }
    mult("z3c", "i1", "r1", 2);
    perm("r1", "(456)");
    update("z3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|6/45) aka (k/ij|c/ab)
 * Diagrams included:
 * T2e T5a T6c T6d T9e
 */
void diag_T3_0h3p_k_ij_c_ab()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }

    tmplt("i1", "pphp", "1000", "2431", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2e
        reorder("pvhp", "r1", "2431");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T5a
        reorder("pphh", "r2", "3142");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r3", "2413");
        mult("r3", "r2", "r4", 2);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T6c
        reorder("pphp", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T6d
        reorder("pvhh", "r1", "2314");
        reorder("t1c", "r2", "21");
        mult("r2", "r1", "r3", 1);
        reorder("r3", "r4", "2134");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T9e
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "3124");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }
    reorder("x3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125346");
    diagram_stack_erase("r3");
    perm("r4", "(3/12|6/45)");
    update("z3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123) aka (ijk)
 * Diagrams included:
 * (ijk):  T6h
 */
void diag_T3_0h3p_ijk()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }

    tmplt("i1", "pphh", "1100", "1234", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T6h
        reorder("vphh", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    reorder("s3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(123)");
    update("z3_0", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (1/23) aka (i/jk)
 * Diagrams included:
 * (i/jk): T2d T5f T9c
 */
void diag_T3_0h3p_i_jk()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }

    tmplt("i1", "pphh", "1100", "1234", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2d
        update("i1", 0.5, "vvhh");

        // T5f
        reorder("pphh", "r2", "3412");
        mult("x2c", "r2", "r3", 2);
        update("i1", 0.25, "r3");
        restore_stack_pos(pos2);

        // T9c
        reorder("pphh", "r2", "3412");
        mult("s1c", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);

        // T6h
        /*reorder("hphh", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        perm("r4", "(12)");
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);*/
    }

    reorder("s3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(1/23)");
    update("z3_0", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12) aka (k/ij)
 * Diagrams included:
 * T2b T5b T6a T6f T9a
 */
void diag_T3_0h3p_k_ij()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }

    tmplt("i1", "ph", "10", "12", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2b
        update("i1", -1.0, "vh");

        // T5b
        reorder("pphh", "r2", "3412");
        mult("s2c", "r2", "r3", 3);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T6a
        reorder("ph", "r2", "21");
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T6f
        reorder("vphh", "r2", "1342");
        mult("r2", "t1c", "r3", 2);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T9a
        reorder("pphh", "r2", "3142");
        mult("r2", "t1c", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    reorder("x3c", "r1", "456123");
    mult("r1", "i1", "r2", 1);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(3/12)");
    update("z3_0", 1.0, "r3");
    restore_stack_pos(pos);
}


/**
 * Folded diagrams for the 0h3p sector:
 * F1: T{03}_3 x Veff{03}
 * F2: T{02}_2 x Veff{02}
 * F3: T{01}_1 T{02}_2 x Veff{02}
 * F4: T{03}_3 x Veff{02}
 * F5: T{01}_1 x Veff{03}
 * F6: T{02}_2 x Veff{03}
 * F7: T{01}_1 S_{01}_1 x Veff{03}
 * F8: T{01}_1 S_{01}_1 S_{01}_1 x Veff{03}
 * F9: T{01}_1 T{02}_2 x Veff{03}
 * F10: T{03}_3 x Veff{03}
 *
 * Diagrams F1, F4, F10 are omitted in the CCSDT-2 and CCSDT-3 models
 * (since T3 operator in amplitude equations is prohibited in these models).
 */
void construct_folded_0h3p()
{
    dg_stack_pos_t pos = get_stack_pos();
    timer_new_entry("folded_Z3", "Folded contribution to T{0h3p}_3");
    timer_start("folded_Z3");

    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        // F1
        reorder("z3c", "r1", "234561");
        mult("veff01", "r1", "r2", 1);
        perm("r2", "(1/23)");
        update("z3nw", -1.0, "r2");
        restore_stack_pos(pos);
    }

    // F2
    reorder("x2c", "r1", "2341");
    mult("veff02", "r1", "r2", 1);
    reorder("r2", "r3", "124356");
    tmplt("r4", "pppppp", "111000", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r3", "r4");
    perm("r4", "(3/12|4/56)");
    update("z3nw", -1.0, "r4");
    restore_stack_pos(pos);

    // F3
    reorder("x2c", "r1", "2341");
    reorder("s1c", "r2", "21");
    reorder("veff02", "r3", "1243");
    mult("r3", "r2", "r4", 1);
    reorder("r4", "r5", "1243");
    mult("r5", "r1", "r6", 1);
    reorder("r6", "r7", "124356");
    diagram_stack_erase("r6");
    perm("r7", "(3/12|4/56)");
    update("z3nw", -1.0, "r7");
    restore_stack_pos(pos);

    // F4
    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        reorder("z3c", "r1", "345612");
        mult("veff02", "r1", "r2", 2);
        perm("r2", "(3/12)");
        update("z3nw", -0.5, "r2");
        restore_stack_pos(pos);
    }

    // F5
    reorder("s1c", "r1", "21");
    mult("veff03", "r1", "r2", 1);
    tmplt("r3", "pppppp", "111000", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r2", "r3");
    perm("r3", "(6/45)");
    update("z3nw", -1.0, "r3");
    restore_stack_pos(pos);

    // F6
    reorder("x2c", "r1", "3412");
    mult("veff03", "r1", "r2", 2);
    tmplt("r3", "pppppp", "111000", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r2", "r3");
    perm("r3", "(4/56)");
    update("z3nw", -0.5, "r3");
    restore_stack_pos(pos);

    // F7
    reorder("s1c", "r1", "21");
    mult("veff03", "r1", "r2", 1);
    reorder("r2", "r3", "123465");
    mult("r3", "r1", "r4", 1);
    reorder("r4", "r5", "123465");
    tmplt("r6", "pppppp", "111000", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r5", "r6");
    diagram_stack_erase("r5");
    perm("r6", "(4/56)");
    update("z3nw", -1.0, "r6");
    restore_stack_pos(pos);

    // F8
    reorder("s1c", "r1", "21");
    mult("veff03", "r1", "r2", 1);
    reorder("r2", "r3", "123645");
    mult("r3", "r1", "r4", 1);
    reorder("r4", "r5", "123465");
    mult("r5", "r1", "r6", 1);
    reorder("r6", "r7", "123654");
    update("z3nw", -1.0, "r7");
    restore_stack_pos(pos);

    // F9
    reorder("s1c", "r1", "21");
    reorder("x2c", "r2", "3412");
    mult("veff03", "r2", "r3", 2);
    reorder("r3", "r4", "123564");
    mult("r4", "r1", "r5", 1);
    reorder("r5", "r6", "123645");
    diagram_stack_erase("r5");
    perm("r6", "(4/56)");
    update("z3nw", -0.5, "r6");
    restore_stack_pos(pos);

    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        // F10
        reorder("z3c", "r1", "456123");
        mult("veff03", "r1", "r2", 3);
        diagram_stack_erase("r1");
        update("z3nw", -0.125, "r2");
        restore_stack_pos(pos);
    }

    timer_stop("folded_Z3");
}
