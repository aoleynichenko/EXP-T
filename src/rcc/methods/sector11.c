/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2023 The EXP-T developers.
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
 * with one hole and one particle (Fock space sector 1h1p)
 * Calculates energies of excitations from the reference state (EE-FS-MRCC)
 *
 * Symbolic names used for diagrams:
 * e1c   current approximation to T^{1h1p}_1 amplitudes
 * e2c   current approximation to T^{1h1p}_2 amplitudes
 * e3c   current approximation to T^{1h1p}_3 amplitudes
 * e1nw  next (new) approximation to T^{1h1p}_1 amplitudes
 * e2nw  next (new) approximation to T^{1h1p}_2 amplitudes
 * e3nw  next (new) approximation to T^{1h1p}_3 amplitudes
 *
 * For the T^{1h1p}_2 operator we exclude amplitudes with four valence indices.
 *
 * E1 operator is not used for calculation of excitation energies, however,
 * it is of extreme importance for properties, thus it is calculated separately.
 *
 * This code was obtained by modification of the CC(0,0) program.
 */

#include "methods.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "ccutils.h"
#include "diis.h"
#include "engine.h"
#include "datamodel.h"
#include "intham_imms.h"
#include "heff.h"
#include "options.h"
#include "sort.h"
#include "utils.h"

void sort_integrals_1h1p();

void init_amplitudes_1h1p();

void const_terms_1h1p();

void const_terms_1h1p_triples(int pt_order);

void construct_singles_1h1p();

void construct_doubles_1h1p();

void construct_triples_1h1p(int pt_order);

void construct_folded_1h1p();

void folded_diagrams_singles_1h1p();

void folded_diagrams_doubles_1h1p();

void folded_diagrams_triples_1h1p(int pt_order);

void t3_contrib_to_singles_1h1p(int pt_order);

void t3_contrib_to_doubles_1h1p(int pt_order);

void diagram_1h1p_D2a_1_const();

void diagram_1h1p_D2a_2();

void diagram_1h1p_D2b_1();

void diagram_1h1p_D2b_2_const();

void diagram_1h1p_D2c_const();

void diagram_1h1p_D2d_const();

void diagram_1h1p_D2e_1();

void diagram_1h1p_D2e_2_const();

void diagram_1h1p_D2e_3_const();

void diagram_1h1p_D2e_4_const();

void diagram_1h1p_D3a_const();

void diagram_1h1p_D3b_1_const();

void diagram_1h1p_D3b_2();

void diagram_1h1p_D3c_1();

void diagram_1h1p_D3c_2_const();

void diagram_1h1p_D3d_1_const();

void diagram_T3_1h1p_D3d_2();

void diagram_1h1p_D4a_1_const();

void diagram_1h1p_D4a_2_const();

void diagram_1h1p_D4b_1_const();

void diagram_1h1p_D4b_2_const();

void diagram_1h1p_D5a_1();

void diagram_1h1p_D5a_2_const();

void diagram_1h1p_D5b_1();

void diagram_1h1p_D5b_2_const();

void diagram_1h1p_D5c_1_const();

void diagram_1h1p_D5c_2_const();

void diagram_1h1p_D5c_3_const();

void diagram_1h1p_D5c_4();

void diagram_1h1p_D5d_1_const();

void diagram_1h1p_D5d_2_const();

void diagram_1h1p_D5d_3_const();

void diagram_1h1p_D5d_4();

void diagram_1h1p_D5e_1_const();

void diagram_1h1p_D5e_2_const();

void diagram_1h1p_D5f_1_const();

void diagram_1h1p_D5f_2_const();

void diagram_1h1p_D5g_1();

void diagram_1h1p_D5g_2_const();

void diagram_1h1p_D5h_1();

void diagram_1h1p_D5h_2_const();

void diagram_1h1p_D6a_const();

void diagram_1h1p_D6b_const();

void diagram_1h1p_D6c_1_const();

void diagram_1h1p_D6c_2_const();

void diagram_1h1p_D6c_3_const();

void diagram_1h1p_D6c_4_const();

void diagram_1h1p_D7a_const();

void diagram_1h1p_D7b_const();

void diagram_1h1p_D7c_1_const();

void diagram_1h1p_D7c_2_const();

void diagram_1h1p_D7c_3_const();

void diagram_1h1p_D7c_4();

void diagram_1h1p_D7d_1_const();

void diagram_1h1p_D7d_2();

void diagram_1h1p_D7e_1_const();

void diagram_1h1p_D7e_2();

void diagram_1h1p_D8a_1_const();

void diagram_1h1p_D8a_2_const();

void diagram_1h1p_D8b_1_const();

void diagram_1h1p_D8b_2_const();

void diagram_1h1p_D9_const();


/**
 * Solves amplitude equations in Fock space sector 1h1p
 */
int sector_1h1p(cc_options_t *opts)
{
    int triples = triples_enabled();

    print_sector_banner(1, 1);
    opts->curr_sector_h = 1;
    opts->curr_sector_p = 1;

    // setup shifts for the IH-like technique by A. V. Zaitsevskii
    if (intham_imms_in_sector(1, 1)) {
        intham_imms_setup(1, 1);
    }

    sort_integrals_1h1p();
    const_terms_1h1p();
    init_amplitudes_1h1p();
    intruder_state_analysis(NULL, "e2c", triples ? "e3c" : NULL);

    /*
     * main loop
     */
    int success = solve_amplitude_equations(
            1, 1,
            "e1c", "e1nw",
            "e2c", "e2nw",
            triples ? "e3c" : NULL, triples ? "e3nw" : NULL,
            "veff11",
            construct_singles_1h1p,
            construct_doubles_1h1p,
            construct_triples_1h1p,
            construct_folded_1h1p
    );
    if (success == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    /*
     * analyze cluster operators and flush them to disk
     */
    print_cluster_operator_analysis(1, 1, "e1c", "e2c", triples ? "e3c" : NULL);
    save_cluster_amplitudes(1, 1, "e1c", "e2c", triples ? "e3c" : NULL, "veff11");
    if (opts->do_flush_amplitudes_txt) {
        save_cluster_amplitudes_formatted(1, 1, "e1c", "e2c", triples ? "e3c" : NULL);
    }

    /*
     * construct and diagonalize effective Hamiltonian
     * analyze its eigenvectors & eigenvalues
     */
    heff_analysis(1, 1, "veff01", "veff10", "veff11");

    /*
     * Ionization potential and electron affinity
     */
    double ip = cc_opts->ground_energy_1h0p;
    double ea = -cc_opts->ground_energy_0h1p;

    print_ionization_potential("Ionization potential 0h0p -> 1h0p", ip);
    print_ionization_potential("Electron affinity    0h0p -> 0h1p", ea);

    /*
     * calculate quasi-natural orbitals and model-space estimates of properties
     */
    calculate_wavefunction_norms_and_overlaps(1, 1);
    calculate_model_space_natural_spinors(1, 1);
    calculate_model_space_properties(1, 1);
    calculate_properties_direct(1, 1);

    return EXIT_SUCCESS;
}


/**
 * Sorting of integrals used in the CC(1h1p) amplitude equations.
 */
void sort_integrals_1h1p()
{
    // sorting -- interaction diagrams with valence lines
    request_sorting("vg", "ph", "11", "12");
    request_sorting("vhpg", "phph", "1001", "1234");
    request_sorting("pppgr", "ppph", "0001", "3412");
    request_sorting("hvhh", "hphh", "0100", "1234");
    request_sorting("phph", "phph", "0000", "1234");
    request_sorting("vphg", "pphh", "1001", "1234");
    request_sorting("pvgh", "pphh", "0110", "1234"); // can be reordered [2143] from 'vphg'
    request_sorting("pvhg", "pphh", "0101", "2431");
    request_sorting("phpg", "phph", "0001", "1234");
    request_sorting("pvgp", "pphp", "0110", "1234");
    request_sorting("vhhg", "phhh", "1001", "1234");
    //request_sorting("vvhh", "pphh", "1100", "1234");
    //request_sorting("vvhp", "pphp", "1100", "1234");
    if (cc_opts->cc_model >= CC_MODEL_CCSD_T3) {
        request_sorting("vhhp", "phhp", "1000", "1234");
        //request_sorting("vphh", "pphh", "1000", "1234");
    }
    perform_sorting();

    // prepare amplitudes from previous sectors
    reorder("s1c", "s1r", "21");
    reorder("s2c", "s2r", "3412");
    reorder("h1c", "h1r", "21");
    reorder("h2c", "h2r", "3412");

    if (cc_opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }
}


/**
 * Constructs initial guess to E2 cluster amplitudes (1h1p sector) or
 * read them from disk (if required)
 */
void init_amplitudes_1h1p()
{
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;
    int calc_t1 = 1, calc_t2 = 1, calc_t3 = 1, calc_veff = 1;

    printf("\n Initialization of T{1h1p}_2 amplitudes ...\n");

    if (cc_opts->reuse_amplitudes[1][1]) {
        printf(" Trying to read amplitudes and effective interaction (sector 1h1p) from disk ...\n");
        // Singles
        if (diagram_read("e1c.dg") != NULL) {
            printf(" T{1h1p}_1 amplitudes successfully read from disk\n");
            calc_t1 = 0;
        }
        else {
            printf(" T{1h1p}_1 amplitudes will be calculated\n");
        }
        // Doubles
        if (diagram_read("e2c.dg") != NULL) {
            printf(" T{1h1p}_2 amplitudes successfully read from disk\n");
            calc_t2 = 0;
        }
        else {
            printf(" T{1h1p}_2 amplitudes will be calculated\n");
        }
        // Effective interaction
        if (diagram_read("veff11.dg") != NULL) {
            printf(" Heff{1h1p} diagram successfully read from disk\n");
            calc_veff = 0;
        }
        else {
            printf(" Heff{1h1p} diagram will be calculated\n");
        }
        if (triples) {
            if (diagram_read("e3c.dg") != NULL) {
                printf(" T{1h1p}_3 amplitudes successfully read from disk\n");
                calc_t3 = 0;
            }
            else {
                printf(" T{1h1p}_3 amplitudes will be calculated\n");
            }
        }
    }

    // init amplitudes if needed
    if (calc_t1) {
        copy("vg", "e1c");
        diveps("e1c");
        printf(" Calculated: T{1h1p}_1\n");
    }
    if (calc_t2 && calc_veff) {
        copy("vhpg", "e2c");
        closed("e2c", "veff11");
        diveps("e2c");
        printf(" Calculated: T{1h1p}_2 and Heff{1h1p}\n");
    }
    else if (calc_t2) {
        copy("vhpg", "e2c");
        closed("e2c", "veff11_reserve");
        diveps("e2c");
        printf(" Calculated: T{1h1p}_2\n");
    }
    else if (calc_veff) {
        copy("vhpg", "veff11_open");
        closed("veff11_open", "veff11");
        printf(" Calculated: Heff{1h1p}\n");
    }
    if (triples && calc_t3) {
        tmplt("e3c", "phhphp", "100010", "123456", IS_PERM_UNIQUE);
        printf(" Calculated: T{1h1p}_3\n");
    }

    if (cc_opts->cc_model == CC_MODEL_CCS) {
        clear("e2c");
    }
    if (cc_opts->cc_model == CC_MODEL_CCD) {
        clear("e1c");
    }

    printf(" done\n\n");
}


/**
 * In sector 0h2p amplitude equations there are diagrams which are independent
 * on the X{01}_2 amplitudes and can be evaluated only once -- at
 * start of the calculation. They include only f, V, T{00}_1, T{00}_2, T{01}_1,
 * T{01}_2 oper-s.
 */
void const_terms_1h1p()
{
    timer_new_entry("const_s11", "Constant part of 1h1p amplitudes");
    timer_start("const_s11");

    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        tmplt("e3_0", "phhphp", "100010", "123456", IS_PERM_UNIQUE);
    }

    copy("vhpg", "e2_0");

    // T2
    diagram_1h1p_D2a_1_const();
    //diagram_1h1p_D2a_2();
    //diagram_1h1p_D2b_1();
    diagram_1h1p_D2b_2_const();
    diagram_1h1p_D2c_const();
    diagram_1h1p_D2d_const();
    //diagram_1h1p_D2e_1();
    diagram_1h1p_D2e_2_const();
    diagram_1h1p_D2e_3_const();
    diagram_1h1p_D2e_4_const();

    // T2^2
    diagram_1h1p_D3a_const();
    diagram_1h1p_D3b_1_const();
    //diagram_1h1p_D3b_2();
    //diagram_1h1p_D3c_1();
    diagram_1h1p_D3c_2_const();
    diagram_1h1p_D3d_1_const();
    //diagram_T3_1h1p_D3d_2();

    // T1
    diagram_1h1p_D4a_1_const();
    diagram_1h1p_D4a_2_const();
    diagram_1h1p_D4b_1_const();
    diagram_1h1p_D4b_2_const();

    // T1 T2
    //diagram_1h1p_D5a_1();
    diagram_1h1p_D5a_2_const();
    //diagram_1h1p_D5b_1();
    diagram_1h1p_D5b_2_const();
    diagram_1h1p_D5c_1_const();
    diagram_1h1p_D5c_2_const();
    diagram_1h1p_D5c_3_const();
    //diagram_1h1p_D5c_4();
    diagram_1h1p_D5d_1_const();
    diagram_1h1p_D5d_2_const();
    diagram_1h1p_D5d_3_const();
    //diagram_1h1p_D5d_4();
    diagram_1h1p_D5e_1_const();
    diagram_1h1p_D5e_2_const();
    diagram_1h1p_D5f_1_const();
    diagram_1h1p_D5f_2_const();
    //diagram_1h1p_D5g_1();
    diagram_1h1p_D5g_2_const();
    //diagram_1h1p_D5h_1();
    diagram_1h1p_D5h_2_const();

    // T1^2
    diagram_1h1p_D6a_const();
    diagram_1h1p_D6b_const();
    diagram_1h1p_D6c_1_const();
    diagram_1h1p_D6c_2_const();
    diagram_1h1p_D6c_3_const();
    diagram_1h1p_D6c_4_const();

    // T1^2 T2
    diagram_1h1p_D7a_const();
    diagram_1h1p_D7b_const();
    diagram_1h1p_D7c_1_const();
    diagram_1h1p_D7c_2_const();
    diagram_1h1p_D7c_3_const();
    //diagram_1h1p_D7c_4();
    diagram_1h1p_D7d_1_const();
    //diagram_1h1p_D7d_2();
    diagram_1h1p_D7e_1_const();
    //diagram_1h1p_D7e_2();

    // T1^3
    diagram_1h1p_D8a_1_const();
    diagram_1h1p_D8a_2_const();
    diagram_1h1p_D8b_1_const();
    diagram_1h1p_D8b_2_const();

    // T1^4
    diagram_1h1p_D9_const();

    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        const_terms_1h1p_triples(PT_INF);
    }

    timer_stop("const_s11");
}


/**
 * Singles equations for the 1h1p sector.
 */
void construct_singles_1h1p()
{
    timer_new_entry("11-E1", "1h1p -- Singles equations (Exc1)");
    timer_start("11-E1");

    // S1
    copy("vg", "e1nw");

    dg_stack_pos_t pos = get_stack_pos();

    // S2a
    reorder("e2c", "e2c_m", "1243");
    reorder("e2c_m", "r1", "1324");
    reorder("ph", "phr", "21");
    mult("r1", "phr", "r2", 2);
    update("e1nw", -1.0, "r2");
    restore_stack_pos(pos);

    // S2b
    reorder("pphg", "r1", "4321");
    mult("s2c", "r1", "r3", 3);
    update("e1nw", 0.5, "r3");
    restore_stack_pos(pos);

    // S2c
    reorder("pvhh", "r1", "2341");
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "r2", "4123");
    mult("r1", "r2", "r3", 3);
    update("e1nw", -0.5, "r3");
    restore_stack_pos(pos);

    // S3a
    reorder("pg", "pgr", "21");
    mult("s1c", "pgr", "r1", 1);
    update("e1nw", 1.0, "r1");
    restore_stack_pos(pos);

    // S3b
    reorder("h1c", "h1_r", "21");
    mult("vh", "h1_r", "r1", 1);
    update("e1nw", -1.0, "r1");
    restore_stack_pos(pos);

    // S3c
    mult("pvhg", "t1c", "r2", 2);
    update("e1nw", 1.0, "r2");
    restore_stack_pos(pos);

    // S4a
    reorder("h2c", "h2_r", "3412");
    mult("h2_r", "pphh", "r1", 3);
    mult("s1c", "r1", "r2", 1);
    update("e1nw", -0.5, "r2");
    restore_stack_pos(pos);

    // S4b
    reorder("h1c", "h1_r", "21");
    reorder("pphh", "vr1", "3412");
    mult("s2c", "vr1", "r1", 3);
    mult("r1", "h1_r", "r2", 1);
    update("e1nw", -0.5, "r2");
    restore_stack_pos(pos);

    // S4c
    reorder("t1c", "t1_r", "21");
    reorder("pphh", "r1", "2413");
    mult("r1", "t1_r", "r2", 2);
    reorder("r2", "r3", "21");
    reorder("e2c", "e2c_m_21", "2134");
    reorder("e2c_m_21", "r4", "2413");
    mult("r4", "r3", "r5", 2);
    update("e1nw", -1.0, "r5");
    restore_stack_pos(pos);

    // S5a
    reorder("h1c", "h1_r", "21");
    reorder("ph", "phr", "21");
    mult("s1c", "phr", "r2", 1);
    mult("r2", "h1_r", "r3", 1);
    update("e1nw", -1.0, "r3");
    restore_stack_pos(pos);

    // S5b
    reorder("pphg", "r1", "4231");
    mult("r1", "t1c", "r2", 2);
    mult("s1c", "r2", "r3", 1);
    update("e1nw", 1.0, "r3");
    restore_stack_pos(pos);

    // S5c
    reorder("h1c", "h1_r", "21");
    reorder("pvhh", "r1", "2431");
    mult("r1", "t1c", "r2", 2);
    mult("r2", "h1_r", "r3", 1);
    update("e1nw", -1.0, "r3");
    restore_stack_pos(pos);

    // S6
    reorder("h1c", "h1_r", "21");
    reorder("t1c", "t1_r", "21");
    reorder("pphh", "r1", "2413");
    mult("r1", "t1_r", "r2", 2);
    reorder("r2", "r3", "21");
    mult("s1c", "r3", "r4", 1);
    mult("r4", "h1_r", "r5", 1);
    update("e1nw", -1.0, "r5");
    restore_stack_pos(pos);

    // Triples contribution to Singles
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_contrib_to_singles_1h1p(PT_INF);
    }

    timer_stop("11-E1");
}


void diagram_1h1p_D2a_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pg", "pgr", "21");
    mult("s2c", "pgr", "r1", 1);
    update("e2_0", 1.0, "r1");

    restore_stack_pos(pos);
}


void diagram_1h1p_D2a_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pp", "ppr_", "21");
    reorder("e2c", "e2_21", "2143");
    mult("e2_21", "ppr_", "r1", 1);
    reorder("r1", "r2", "2143");
    update("e2nw", 1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h1p_D2b_1()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("e2c", "e2r_", "3412");
    mult("e2r_", "hh", "r1", 1);
    reorder("r1", "r2", "3412");
    update("e2nw", -1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h1p_D2b_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    // D2b-2
    reorder("h2c", "h2r_", "3412");
    mult("h2r_", "vh", "r1", 1);
    reorder("r1", "r2", "3412");
    reorder("r2", "r3", "2143");
    update("e2_0", -1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D2c_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    mult("pppgr", "s2c", "r1", 2);
    reorder("r1", "r2", "3412");
    update("e2_0", 0.5, "r2");

    restore_stack_pos(pos);
}


void diagram_1h1p_D2d_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("h2c", "r1", "3412");
    mult("hvhh", "r1", "r2", 2);
    reorder("r2", "r3", "2143");
    update("e2_0", 0.5, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D2e_1()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("e2c", "r1", "1423");
    reorder("phph", "r2", "2341");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1342");
    update("e2nw", -1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h1p_D2e_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("t2c", "r1", "1423");
    reorder("pvgh", "r2", "2341");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1342");
    reorder("r4", "r5", "2143");
    update("e2_0", -1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h1p_D2e_3_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("s2c", "r1", "1324");
    reorder("phhg", "r0", "2431");
    mult("r1", "r0", "r3", 2);
    reorder("r3", "r4", "1324");
    update("e2_0", 1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h1p_D2e_4_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("h2c", "r1", "1324");
    reorder("pvhp", "r0", "2431");
    mult("r1", "r0", "r3", 2);
    reorder("r3", "r4", "1324");
    reorder("r4", "r5", "2143");
    update("e2_0", 1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h1p_D3a_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "h2r_", "3412");
    reorder("pphh", "v1", "3412");
    mult("s2c", "v1", "r1", 2);
    mult("r1", "h2r_", "r2", 2);
    update("e2_0", 0.25, "r2");

    restore_stack_pos(pos);
}


void diagram_1h1p_D3b_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("s2c", "s2_r1", "1324");
    reorder("h2c", "h2_r1", "1324");
    reorder("pphh", "r2", "4231");
    mult("s2_r1", "r2", "r3", 2);
    mult("r3", "h2_r1", "r4", 2);
    reorder("r4", "r5", "1324");
    update("e2_0", 1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h1p_D3b_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("e2c", "e2c_m", "1243");
    reorder("e2c_m", "e2_r1", "1324");
    reorder("t2c", "t2_r1", "1324");
    reorder("pphh", "r2", "4231");
    mult("e2_r1", "r2", "r3", 2);
    mult("r3", "t2_r1", "r4", 2);
    reorder("r4", "r5", "1324");
    reorder("r5", "r6", "1243");
    update("e2nw", 1.0, "r6");

    restore_stack_pos(pos);
}


void diagram_1h1p_D3c_1()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 3);
    reorder("e2c", "e2_21", "2143");
    reorder("e2_21", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    reorder("r3", "r4", "2143");
    update("e2nw", -0.5, "r4");

    restore_stack_pos(pos);
}


void diagram_1h1p_D3c_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v1", "3412");
    mult("s2c", "v1", "r1", 3);
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    update("e2_0", -0.5, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D3d_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("h2c", "h2r_", "3412");
    mult("h2r_", "pphh", "r1", 3);
    mult("s2c", "r1", "r2", 1);
    update("e2_0", -0.5, "r2");

    restore_stack_pos(pos);
}


void diagram_T3_1h1p_D3d_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("t2c", "t2r_", "3412");
    mult("t2r_", "pphh", "r1", 3);
    reorder("e2c", "e2_21", "2143");
    mult("e2_21", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("e2nw", -0.5, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D4a_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phpg", "r1", "2341");
    mult("s1c", "r1", "r2", 1);
    update("e2_0", 1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h1p_D4a_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pvgp", "r1", "2341");
    mult("t1c", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("e2_0", 1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D4b_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("hvhp", "r1", "1243");
    reorder("h1c", "h1_21", "21");
    mult("r1", "h1_21", "r2", 1);
    reorder("r2", "r3", "2134");
    update("e2_0", -1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D4b_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();


    reorder("vhhg", "r1", "1243");
    reorder("t1c", "t1_21", "21");
    mult("r1", "t1_21", "r2", 1);
    reorder("r2", "r3", "1243");
    update("e2_0", -1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5a_1()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("ph", "r1", "21");
    reorder("e2c", "e2_21", "3412");
    mult("t1c", "r1", "r2", 1);
    mult("e2_21", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("e2nw", -1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5a_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("ph", "r1", "21");
    reorder("h2c", "h2_21", "3412");
    mult("s1c", "r1", "r2", 1);
    mult("h2_21", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    reorder("r4", "r5", "2143");
    update("e2_0", -1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5b_1()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("t1c", "t1_21", "21");
    mult("t1_21", "ph", "r1", 1);
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("e2nw", -1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5b_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("h1c", "h1_21", "21");
    mult("h1_21", "ph", "r1", 1);
    mult("s2c", "r1", "r2", 1);
    update("e2_0", -1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5c_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphg", "r1", "4312");
    mult("s1c", "r1", "r2", 1);
    reorder("t2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    reorder("r5", "r6", "2134");
    update("e2_0", -1.0, "r6");   // minus

    restore_stack_pos(pos);
}


void diagram_1h1p_D5c_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphp", "r1", "4312");
    mult("s1c", "r1", "r2", 1);
    reorder("h2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    reorder("r5", "r6", "2143");
    update("e2_0", 1.0, "r6");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5c_3_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphg", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    reorder("s2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    update("e2_0", 1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5c_4()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("e2c", "e2c_m", "1243");
    reorder("pphp", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    reorder("e2c_m", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    reorder("r5", "r6", "1243");
    update("e2nw", 1.0, "r6");   // minus

    restore_stack_pos(pos);
}


void diagram_1h1p_D5d_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pvhh", "r1", "2314");
    reorder("h1c", "h1_21", "21");
    mult("h1_21", "r1", "i1", 1);
    reorder("t2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    reorder("r3", "r4", "2134");
    update("e2_0", 1.0, "r4");   // minus

    restore_stack_pos(pos);
}


void diagram_1h1p_D5d_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pvhh", "r1", "2314");
    reorder("t1c", "t1_21", "21");
    mult("t1_21", "r1", "i1", 1);
    reorder("h2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    reorder("r3", "r4", "2143");
    update("e2_0", -1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5d_3_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phhh", "r1", "2314");
    reorder("h1c", "h1_21", "21");
    mult("h1_21", "r1", "i1", 1);
    reorder("s2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    update("e2_0", -1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5d_4()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phhh", "r1", "2314");
    reorder("t1c", "t1_21", "21");
    mult("t1_21", "r1", "i1", 1);
    reorder("e2c", "e2c_m", "1243");
    reorder("e2c_m", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    reorder("r3", "r4", "1243");
    update("e2nw", -1.0, "r4");   // minus

    restore_stack_pos(pos);
}


void diagram_1h1p_D5e_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphp", "r1", "4312");
    reorder("h1c", "h1_21", "21");
    mult("s2c", "r1", "i1", 2);
    mult("i1", "h1_21", "r2", 1);
    update("e2_0", 0.5, "r2");   // minus

    restore_stack_pos(pos);
}


void diagram_1h1p_D5e_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphg", "r1", "4312");
    reorder("t1c", "t1_21", "21");
    mult("s2c", "r1", "i1", 2);
    mult("i1", "t1_21", "r2", 1);
    reorder("r2", "r3", "1243");
    update("e2_0", -0.5, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5f_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phhh", "r1", "2341");
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "h2r", "3412");
    mult("s1c", "r1", "i1", 1);
    mult("i1", "h2r", "r2", 2);
    update("e2_0", 0.5, "r2");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5f_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pvhh", "r1", "2341");
    reorder("h2c", "h2r", "3412");
    mult("t1c", "r1", "i1", 1);
    mult("i1", "h2r", "r2", 2);
    reorder("r2", "r3", "2143");
    update("e2_0", 0.5, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5g_1()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphp", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("e2nw", 1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5g_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphg", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("s2c", "i1", "r2", 1);
    update("e2_0", 1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5h_1()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("e2c", "e2c_21", "2143");
    reorder("e2c_21", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    reorder("r3", "r4", "2143");
    update("e2nw", -1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h1p_D5h_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pvhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    update("e2_0", -1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D6a_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    mult("pppgr", "t1c", "i1", 1);
    reorder("i1", "r1", "4123");
    mult("s1c", "r1", "r2", 1);
    update("e2_0", 1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h1p_D6b_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("t1c", "t1_21", "21");
    reorder("h1c", "h1_21", "21");
    mult("h1_21", "vhhh", "i1", 1);
    mult("t1_21", "i1", "i2", 1);
    reorder("i2", "r1", "3412");
    update("e2_0", 1.0, "r1");

    restore_stack_pos(pos);
}


void diagram_1h1p_D6c_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phph", "r1", "2341");
    reorder("h1c", "r2", "21");
    mult("s1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    update("e2_0", -1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h1p_D6c_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pvgh", "r1", "2341");
    reorder("t1c", "r2", "21");
    mult("t1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    reorder("r4", "r5", "2143");
    update("e2_0", -1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h1p_D6c_3_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phhg", "r0", "2431");
    mult("s1c", "r0", "i1", 1);
    reorder("t1c", "t1_21", "21");
    mult("i1", "t1_21", "i2", 1);
    reorder("i2", "r2", "1243");
    update("e2_0", -1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h1p_D6c_4_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pvhp", "r0", "2431");
    mult("t1c", "r0", "i1", 1);
    reorder("h1c", "h1_21", "21");
    mult("i1", "h1_21", "i2", 1);
    reorder("i2", "r2", "1243");
    reorder("r2", "r3", "2143");
    update("e2_0", -1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D7a_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "vr1", "3412");
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "h2_21_r", "3412");
    mult("t1c", "vr1", "r1", 1);
    mult("s1c", "r1", "r2", 1);
    mult("r2", "h2_21_r", "r3", 2);
    update("e2_0", 0.5, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D7b_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "vr1", "3412");
    reorder("t1c", "t1_r", "21");
    reorder("h1c", "h1_r", "21");
    mult("s2c", "vr1", "r1", 2);
    mult("h1_r", "r1", "r2", 1);
    mult("t1_r", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("e2_0", 0.5, "r4");

    restore_stack_pos(pos);
}


void diagram_1h1p_D7c_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "r1", "3142");
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "r2", "2413");
    reorder("t1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("s1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    reorder("r6", "r7", "1243");
    update("e2_0", -1.0, "r7");

    restore_stack_pos(pos);
}


void diagram_1h1p_D7c_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "r1", "3142");
    reorder("s2c", "s2_21", "2143");
    reorder("s2_21", "r2", "2413");
    reorder("h1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("t1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    reorder("r6", "r7", "1243");
    reorder("r7", "r8", "2143");
    update("e2_0", -1.0, "r8");

    restore_stack_pos(pos);
}


void diagram_1h1p_D7c_3_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "r1", "3142");
    reorder("t2c", "r2", "2413");
    reorder("h1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("s1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    update("e2_0", 1.0, "r6");   // minus

    restore_stack_pos(pos);
}


void diagram_1h1p_D7c_4()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("e2c", "e2c_m_21", "2134");
    reorder("pphh", "r1", "3142");
    reorder("e2c_m_21", "r2", "2413");
    reorder("t1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("t1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    reorder("r6", "r7", "2143");
    update("e2nw", -1.0, "r7");   // minus

    restore_stack_pos(pos);
}


void diagram_1h1p_D7d_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v_", "2431");
    mult("v_", "t1c", "r1", 2);
    reorder("r1", "r2", "21");
    mult("s1c", "r2", "r3", 1);
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "h2_21_r", "2341");
    mult("r3", "h2_21_r", "r4", 1);
    update("e2_0", -1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h1p_D7d_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v_", "2431");
    mult("v_", "t1c", "r1", 2);
    reorder("r1", "r2", "21");
    mult("t1c", "r2", "r3", 1);
    reorder("e2c", "e2c_21", "2143");
    reorder("e2c_21", "e2c_21_r", "2341");
    mult("r3", "e2c_21_r", "r4", 1);
    reorder("r4", "r5", "2143");
    update("e2nw", -1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h1p_D7e_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v_", "2431");
    reorder("h1c", "h1_21", "21");
    mult("v_", "t1c", "r1", 2);
    mult("h1_21", "r1", "r2", 1);
    reorder("s2c", "s2_21", "2143");
    reorder("s2_21", "s2_21_r", "1243");
    mult("s2_21_r", "r2", "r3", 1);
    reorder("r3", "r4", "2134");
    update("e2_0", -1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h1p_D7e_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v_", "2431");
    reorder("t1c", "t1_21", "21");
    mult("v_", "t1c", "r1", 2);
    mult("t1_21", "r1", "r2", 1);
    reorder("e2c", "e2r_", "1243");
    mult("e2r_", "r2", "r3", 1);
    reorder("r3", "r4", "1243");
    update("e2nw", -1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h1p_D8a_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphg", "r1", "4312");
    reorder("t1c", "t1_21", "21");
    mult("t1c", "r1", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    mult("r3", "t1_21", "r4", 1);
    reorder("r4", "r5", "1243");
    update("e2_0", -1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h1p_D8a_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphp", "r1", "4312");
    reorder("h1c", "h1_21", "21");
    mult("s1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    mult("r3", "h1_21", "r4", 1);
    reorder("r4", "r5", "2134");
    update("e2_0", -1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h1p_D8b_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phhh", "r1", "2341");
    reorder("t1c", "t1_21", "21");
    reorder("h1c", "h1_21", "21");
    mult("s1c", "r1", "i1", 1);
    mult("h1_21", "i1", "i2", 1);
    mult("t1_21", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    update("e2_0", 1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h1p_D8b_2_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pvhh", "r1", "2341");
    reorder("t1c", "t1_21", "21");
    reorder("h1c", "h1_21", "21");
    mult("t1c", "r1", "i1", 1);
    mult("t1_21", "i1", "i2", 1);
    mult("h1_21", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    reorder("r3", "r4", "2143");
    update("e2_0", 1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h1p_D9_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("t1c", "t1_21", "21");
    reorder("h1c", "h1_21", "21");
    reorder("pphh", "vr1", "3412");
    mult("t1c", "vr1", "r1", 1);
    mult("s1c", "r1", "r2", 1);
    mult("h1_21", "r2", "r3", 1);
    mult("t1_21", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    update("e2_0", 1.0, "r5");

    restore_stack_pos(pos);
}


/**
 * Constructs the next approximation to T{11}_2 amplitudes ("Doubles").
 * Only T{11}_2/3 dependent diagrams will be calculated; for independent terms
 * see const_terms_1h1p().
 */
void construct_doubles_1h1p()
{
    timer_new_entry("11-S2", "1h1p -- Doubles equations (T{1h1p}_2)");
    timer_start("11-S2");

    // diagram D1
    copy("e2_0", "e2nw");

    // T2
    //diagram_1h1p_D2a_1_const();
    diagram_1h1p_D2a_2();
    diagram_1h1p_D2b_1();
    //diagram_1h1p_D2b_2_const();
    //diagram_1h1p_D2c_const();
    //diagram_1h1p_D2d_const();
    diagram_1h1p_D2e_1();
    //diagram_1h1p_D2e_2_const();
    //diagram_1h1p_D2e_3_const();
    //diagram_1h1p_D2e_4_const();

    // T2^2
    //diagram_1h1p_D3a_const();
    //diagram_1h1p_D3b_1_const();
    diagram_1h1p_D3b_2();
    diagram_1h1p_D3c_1();
    //diagram_1h1p_D3c_2_const();
    //diagram_1h1p_D3d_1_const();
    diagram_T3_1h1p_D3d_2();

    // T1
    //diagram_1h1p_D4a_1_const();
    //diagram_1h1p_D4a_2_const();
    //diagram_1h1p_D4b_1_const();
    //diagram_1h1p_D4b_2_const();

    // T1 T2
    diagram_1h1p_D5a_1();
    //diagram_1h1p_D5a_2_const();
    diagram_1h1p_D5b_1();
    //diagram_1h1p_D5b_2_const();
    //diagram_1h1p_D5c_1_const();
    //diagram_1h1p_D5c_2_const();
    //diagram_1h1p_D5c_3_const();
    diagram_1h1p_D5c_4();
    //diagram_1h1p_D5d_1_const();
    //diagram_1h1p_D5d_2_const();
    //diagram_1h1p_D5d_3_const();
    diagram_1h1p_D5d_4();
    //diagram_1h1p_D5e_1_const();
    //diagram_1h1p_D5e_2_const();
    //diagram_1h1p_D5f_1_const();
    //diagram_1h1p_D5f_2_const();
    diagram_1h1p_D5g_1();
    //diagram_1h1p_D5g_2_const();
    diagram_1h1p_D5h_1();
    //diagram_1h1p_D5h_2_const();

    // T1^2
    //diagram_1h1p_D6a_const();
    //diagram_1h1p_D6b_const();
    //diagram_1h1p_D6c_1_const();
    //diagram_1h1p_D6c_2_const();
    //diagram_1h1p_D6c_3_const();
    //diagram_1h1p_D6c_4_const();

    // T1^2 T2
    //diagram_1h1p_D7a_const();
    //diagram_1h1p_D7b_const();
    //diagram_1h1p_D7c_1_const();
    //diagram_1h1p_D7c_2_const();
    //diagram_1h1p_D7c_3_const();
    diagram_1h1p_D7c_4();
    //diagram_1h1p_D7d_1_const();
    diagram_1h1p_D7d_2();
    //diagram_1h1p_D7e_1_const();
    diagram_1h1p_D7e_2();

    // T1^3
    //diagram_1h1p_D8a_1_const();
    //diagram_1h1p_D8a_2_const();
    //diagram_1h1p_D8b_1_const();
    //diagram_1h1p_D8b_2_const();

    // T1^4
    //diagram_1h1p_D9_const();

    // Triples contribution to Doubles
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_contrib_to_doubles_1h1p(PT_INF);
    }

    timer_stop("11-S2");
}


/**
 * Evaluates folded diagrams in the (1,1) sector.
 */
void construct_folded_1h1p()
{
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;

    folded_diagrams_singles_1h1p();
    folded_diagrams_doubles_1h1p();
    if (triples) {
        folded_diagrams_triples_1h1p(PT_INF);
    }
}


/**
 * Folded diagrams in the 1h1p sector.
 */
void folded_diagrams_singles_1h1p()
{
    dg_stack_pos_t pos = get_stack_pos();

    // Singles:
    // SF1
    reorder("e1c", "r1", "21");
    mult("veff01", "r1", "r2", 1);
    update("e1nw", -1.0, "r2");
    restore_stack_pos(pos);

    // SF2
    reorder("veff10", "r1", "21");
    mult("e1c", "r1", "r2", 1);
    update("e1nw", +1.0, "r2");
    restore_stack_pos(pos);

    // SF3
    reorder("veff11", "r1", "1432");
    mult("r1", "e1c", "r2", 2);
    update("e1nw", +1.0, "r2");
    restore_stack_pos(pos);
}


void folded_diagrams_doubles_1h1p()
{
    dg_stack_pos_t pos = get_stack_pos();

    // Doubles:
    // F1
    reorder("e2c", "r1", "2314");
    reorder("veff11", "r2", "1432");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "3124");
    update("e2nw", +1.0, "r4");
    restore_stack_pos(pos);

    // F2
    reorder("s1c", "r1", "21");
    reorder("veff11", "r2", "4123");
    mult("r1", "r2", "r3", 1);
    mult("h1c", "r3", "r4", 1);
    reorder("r4", "r5", "4123");
    update("e2nw", +1.0, "r5");
    restore_stack_pos(pos);

    // F3
    reorder("e2c", "r1", "2341");
    mult("veff01", "r1", "r2", 1);
    update("e2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // F4
    reorder("veff10", "r1", "21");
    mult("e2c", "r1", "r2", 1);
    update("e2nw", +1.0, "r2");
    restore_stack_pos(pos);

    // F5
    reorder("veff11", "r1", "1243");
    reorder("s1c", "r2", "21");
    mult("r1", "r2", "r3", 1);
    reorder("r3", "r4", "1243");
    tmplt("r5", "phph", "1001", "1234", NOT_PERM_UNIQUE);
    expand_diagram("r4", "r5");
    update("e2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // F6
    reorder("veff11", "r1", "1342");
    mult("r1", "h1c", "r2", 1);
    reorder("r2", "r3", "1423");
    tmplt("r4", "phph", "1001", "1234", NOT_PERM_UNIQUE);
    expand_diagram("r3", "r4");
    update("e2nw", +1.0, "r4");
    restore_stack_pos(pos);
}


