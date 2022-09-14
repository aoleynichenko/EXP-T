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
 * with one hole (Fock space sector 1h0p)
 *
 * Symbolic names used for diagrams:
 * h1c   current approximation to T^{1h0p}_1 amplitudes
 * h2c   current approximation to T^{1h0p}_2 amplitudes
 * h3c   current approximation to T^{1h0p}_3 amplitudes
 * h1nw  next (new) approximation to T^{1h0p}_1 amplitudes
 * h2nw  next (new) approximation to T^{1h0p}_2 amplitudes
 * h3nw  next (new) approximation to T^{1h0p}_3 amplitudes
 *
 * This code was obtained by modification of the CC(0,0) program:
 * all we need is to flip particle creation lines to turn them to valence hole
 * annihilation lines. See Kaldor, J. Comp. Chem. V. 8, P.448 (1987) for details.
 *
 * 2019-2022 Alexander Oleynichenko
 */

#include "methods.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ccutils.h"
#include "diis.h"
#include "engine.h"
#include "datamodel.h"
#include "heff.h"
#include "linalg.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"

void sort_integrals_1h0p();

void init_amplitudes_1h0p();

void const_terms_1h0p();

void construct_singles_1h0p();

void construct_doubles_1h0p();

void construct_triples_1h0p();

void construct_folded_1h0p();

void t3_contrib_to_folded_1h0p();

void t3_contrib_to_singles_1h0p();

void t3_contrib_to_doubles_1h0p();

void diagram_1h0p_D1_const();

void diagram_1h0p_D2a_1_const();

void diagram_1h0p_D2a_2();

void diagram_1h0p_D2b();

void diagram_1h0p_D2c_const();

void diagram_1h0p_D2d();

void diagram_1h0p_D2e_1_const();

void diagram_1h0p_D2e_2();

void diagram_1h0p_D3a();

void diagram_1h0p_D3b();

void diagram_1h0p_D3c();

void diagram_1h0p_D3d_1();

void diagram_1h0p_D3d_2();

void diagram_1h0p_D4a_const();

void diagram_1h0p_D4b_1_const();

void diagram_1h0p_D4b_2();

void diagram_1h0p_D5a();

void diagram_1h0p_D5b_1();

void diagram_1h0p_D5b_2();

void diagram_1h0p_D5c_1_const();

void diagram_1h0p_D5c_2();

void diagram_1h0p_D5d_1();

void diagram_1h0p_D5d_2();

void diagram_1h0p_D5e_1_const();

void diagram_1h0p_D5e_2();

void diagram_1h0p_D5f();

void diagram_1h0p_D5g_1_const();

void diagram_1h0p_D5g_2();

void diagram_1h0p_D5h();

void diagram_1h0p_D6a_const();

void diagram_1h0p_D6b();

void diagram_1h0p_D6c_1_const();

void diagram_1h0p_D6c_2();

void diagram_1h0p_D7a();

void diagram_1h0p_D7b();

void diagram_1h0p_D7c_1();

void diagram_1h0p_D7c_2();

void diagram_1h0p_D7d();

void diagram_1h0p_D7e_1();

void diagram_1h0p_D7e_2();

void diagram_1h0p_D8a_1_const();

void diagram_1h0p_D8a_2();

void diagram_1h0p_D8b();

void diagram_1h0p_D9();


/**
 * Solves amplitude equations and find correlation energy for ground state
 * (Fock space sector 1h0p)
 */
int sector_1h0p(cc_options_t *opts)
{
    int triples = triples_enabled();

    print_sector_banner(1, 0);
    opts->curr_sector_h = 1;
    opts->curr_sector_p = 0;

    // setup shifts for the IH-like technique by A. V. Zaitsevskii
    if (intham_imms_in_sector(1, 0)) {
        intham_imms_setup(1, 0);
    }

    sort_integrals_1h0p();
    const_terms_1h0p();
    init_amplitudes_1h0p();

    /*
     * main loop
     */
    int success = solve_amplitude_equations(
            1, 0,
            "h1c", "h1nw",
            "h2c", "h2nw",
            triples ? "h3c" : NULL, triples ? "h3nw" : NULL,
            "veff10",
            construct_singles_1h0p,
            construct_doubles_1h0p,
            construct_triples_1h0p,
            construct_folded_1h0p
    );
    if (success == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    /*
     * analyze cluster amplitudes and flush them to disk
     */
    print_cluster_operator_analysis(1, 0, "h1c", "h2c", triples ? "h3c" : NULL);
    save_cluster_amplitudes(1, 0, "h1c", "h2c", triples ? "h3c" : NULL, "veff10");
    if (opts->do_flush_amplitudes_txt) {
        save_cluster_amplitudes_formatted(1, 0, "h1c", "h2c", triples ? "h3c" : NULL);
    }

    /*
     * construct and diagonalize effective Hamiltonian
     * analyze its eigenvectors & eigenvalues
     */
    heff_analysis(1, 0, "veff10");

    /*
     * Ionization potentials
     */
    double ion_pot = cc_opts->ground_energy_1h0p;
    print_ionization_potential("Ionization potential 0h0p -> 1h0p", ion_pot);

    /*
     * calculate quasi-natural orbitals and model-space estimates of properties
     */
    calculate_model_space_properties(1, 0);
    calculate_model_space_natural_spinors(1, 0);

    return EXIT_SUCCESS;
}


/**
 * Sorting of integrals used in the CC(1h0p) amplitude equations.
 */
void sort_integrals_1h0p()
{
    // prepare one-electron integrals
    request_sorting("hg", "hh", "01", "12");
    request_sorting("pg", "ph", "01", "12");

    // prepare two-electron integrals
    request_sorting("pphg", "pphh", "0001", "1234");
    request_sorting("phhg", "phhh", "0001", "1234");
    request_sorting("hhgp", "hhhp", "0010", "1234");
    request_sorting("phgp", "phhp", "0010", "1234");
    request_sorting("hhhg", "hhhh", "0001", "1234");
    request_sorting("ppgp", "pphp", "0010", "3412");

    // special for triples
    if (cc_opts->cc_model >= CC_MODEL_CCSD_T3) {
        request_sorting("ppgh", "pphh", "0010", "1234");
    }
    perform_sorting();

    // prepare cluster amplitudes from the 0h0p sector
    reorder("t1c", "t1r", "21");
    reorder("t2c", "t2r", "3412");
    reorder("pp", "ppr", "21");

    if (cc_opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }
}


/**
 * Constructs initial guess to cluster amplitudes (1h0p sector) or
 * read them from disk (if required)
 */
void init_amplitudes_1h0p()
{
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;
    int calc_s1 = 1, calc_s2 = 1, calc_s3 = 1, calc_veff = 1;

    printf("\n Initialization of T{1h0p}_1 and T{1h0p}_2 amplitudes ...\n");

    if (cc_opts->reuse_amplitudes[1][0]) {
        printf(" Trying to read amplitudes from disk ...\n");
        // Singles
        if (diagram_read("h1c.dg") != NULL) {
            printf(" T{1h0p}_1 amplitudes successfully read from disk\n");
            calc_s1 = 0;
        }
        else {
            printf(" T{1h0p}_1 amplitudes will be calculated\n");
        }
        // Doubles
        if (diagram_read("h2c.dg") != NULL) {
            printf(" T{1h0p}_2 amplitudes successfully read from disk\n");
            calc_s2 = 0;
        }
        else {
            printf(" T{1h0p}_2 amplitudes will be calculated\n");
        }
        // Effective interaction
        if (diagram_read("veff10.dg") != NULL) {
            printf(" Heff{10} diagram successfully read from disk\n");
            calc_veff = 0;
        }
        else {
            printf(" Heff{10} diagram will be calculated\n");
        }
        // maybe Triples
        if (triples) {
            if (diagram_read("h3c.dg") != NULL) {
                printf(" T{1h0p}_3 amplitudes successfully read from disk\n");
                calc_s3 = 0;
            }
            else {
                printf(" T{1h0p}_3 amplitudes will be calculated\n");
            }
        }
    }

    // calculate initial approximation to T{1h0p} amplitudes (if needed)
    if (calc_s1) {
        copy("h1_0", "h1c");
        closed("h1c", "veff10_reserve");
        diveps("h1c");
    }
    if (calc_s2) {
        copy("h2_0", "h2c");
        diveps("h2c");
    }
    if (calc_veff) {
        copy("h1_0", "veff10_open");
        closed("veff10_open", "veff10");
    }
    if (triples && calc_s3) {
        tmplt("h3c", "hhhhpp", "000100", "123456", IS_PERM_UNIQUE);
    }
    printf("\n");

    if (cc_opts->cc_model == CC_MODEL_CCS) {
        clear("h2c");
    }
    else if (cc_opts->cc_model == CC_MODEL_CCD) {
        clear("h1c");
    }
}


/**
 * In sector 1h0p amplitude equations there are diagrams which are independent
 * on the T{1h0p}_1 and T{1h0p}_2 amplitudes and can be evaluated only once -- at
 * start of the calculation. They include only f, V, T{00}_1 and T{00}_2 oper-s.
 */
void const_terms_1h0p()
{
    printf(" Construction of T^(1h0p)-independent contributions to the IPCCSD-equations ...\n");

    timer_new_entry("const_s10", "Constant part of 1h0p amplitudes");
    timer_start("const_s10");

    tmplt("h1_0", "hh", "01", "12", NOT_PERM_UNIQUE);
    tmplt("h2_0", "hhhp", "0010", "1234", IS_PERM_UNIQUE);

    diagram_1h0p_D1_const();

    diagram_1h0p_D2a_1_const();
    //diagram_1h0p_D2a_2();
    //diagram_1h0p_D2b();
    diagram_1h0p_D2c_const();
    //diagram_1h0p_D2d();
    diagram_1h0p_D2e_1_const();
    //diagram_1h0p_D2e_2();

    //diagram_1h0p_D3a();
    //diagram_1h0p_D3b();
    //diagram_1h0p_D3c();
    //diagram_1h0p_D3d_1();
    //diagram_1h0p_D3d_2();

    diagram_1h0p_D4a_const();
    diagram_1h0p_D4b_1_const();
    //diagram_1h0p_D4b_2();

    //diagram_1h0p_D5a();
    //diagram_1h0p_D5b_1();
    //diagram_1h0p_D5b_2();
    diagram_1h0p_D5c_1_const();
    //diagram_1h0p_D5c_2();
    //diagram_1h0p_D5d_1();
    //diagram_1h0p_D5d_2();
    diagram_1h0p_D5e_1_const();
    //diagram_1h0p_D5e_2();
    //diagram_1h0p_D5f();
    diagram_1h0p_D5g_1_const();
    //diagram_1h0p_D5g_2();
    //diagram_1h0p_D5h();

    diagram_1h0p_D6a_const();
    //diagram_1h0p_D6b();
    diagram_1h0p_D6c_1_const();
    //diagram_1h0p_D6c_2();

    //diagram_1h0p_D7a();
    //diagram_1h0p_D7b();
    //diagram_1h0p_D7c_1();
    //diagram_1h0p_D7c_2();
    //diagram_1h0p_D7d();
    //diagram_1h0p_D7e_1();
    //diagram_1h0p_D7e_2();

    diagram_1h0p_D8a_1_const();
    //diagram_1h0p_D8a_2();
    //diagram_1h0p_D8b();

    //diagram_1h0p_D9();

    timer_stop("const_s10");
}


/**
 * Constructs the next approximation to T{1h0p}_1 amplitudes ("Singles").
 */
void construct_singles_1h0p()
{
    // prepare amplitudes
    reorder("h2c", "h2cr", "3412");
    reorder("h1c", "h1cr", "21");

    copy("h1_0", "h1nw");
    dg_stack_pos_t pos = get_stack_pos();

    // S1
    copy("hg", "h1nw");

    // S2a
    reorder("h2cr", "r1", "1324");
    mult("r1", "ph", "r2", 2);
    reorder("r2", "r3", "21");
    update("h1nw", 1.0, "r3");
    restore_stack_pos(pos);

    // S2b
    reorder("pphg", "r1", "4321");
    reorder("t2c", "r2", "2143");
    mult("r2", "r1", "r3", 3);
    update("h1nw", 0.5, "r3");
    restore_stack_pos(pos);

    // S2c
    reorder("h2c", "r0", "2143");
    reorder("r0", "r2", "4123");
    reorder("phhh", "r1", "2341");
    mult("r1", "r2", "r3", 3);
    update("h1nw", -0.5, "r3");
    restore_stack_pos(pos);

    // S3a
    reorder("pg", "pgr", "21");
    mult("t1c", "pgr", "r1", 1);
    update("h1nw", 1.0, "r1");
    restore_stack_pos(pos);

    // S3b
    mult("hh", "h1cr", "r1", 1);
    update("h1nw", -1.0, "r1");
    restore_stack_pos(pos);

    // S3c
    reorder("phhg", "r1", "2431");
    mult("r1", "t1c", "r2", 2);
    update("h1nw", 1.0, "r2");
    restore_stack_pos(pos);

    // S4a
    mult("h2cr", "pphh", "r1", 3);
    mult("t1c", "r1", "r2", 1);
    update("h1nw", -0.5, "r2");
    restore_stack_pos(pos);

    // S4b
    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 3);
    mult("r1", "h1cr", "r2", 1);
    update("h1nw", -0.5, "r2");
    restore_stack_pos(pos);

    // S4c
    reorder("pphh", "r1", "2413");
    mult("r1", "t1r", "r2", 2);
    reorder("r2", "r3", "2134");
    reorder("h2c", "r4", "1324");
    mult("r4", "r3", "r5", 2);
    update("h1nw", 1.0, "r5");
    restore_stack_pos(pos);

    // S5a
    reorder("ph", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    mult("r2", "h1cr", "r3", 1);
    update("h1nw", -1.0, "r3");
    restore_stack_pos(pos);

    // S5b
    reorder("pphg", "r1", "4231");
    mult("r1", "t1c", "r2", 2);
    mult("t1c", "r2", "r3", 1);
    update("h1nw", 1.0, "r3");
    restore_stack_pos(pos);

    // S5c
    reorder("phhh", "r1", "2431");
    mult("r1", "t1c", "r2", 2);
    mult("r2", "h1cr", "r3", 1);
    update("h1nw", -1.0, "r3");
    restore_stack_pos(pos);

    // S6c
    reorder("pphh", "r1", "2413");
    mult("r1", "t1r", "r2", 2);
    reorder("r2", "r3", "21");
    mult("t1c", "r3", "r4", 1);
    mult("r4", "h1cr", "r5", 1);
    update("h1nw", -1.0, "r5");
    restore_stack_pos(pos);

    // Triples contribution to Singles; for CCSDT-1a, ... only
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_contrib_to_singles_1h0p();
    }
}


/**
 * Constructs the next approximation to T{1h0p}_2 amplitudes ("Doubles").
 */
void construct_doubles_1h0p()
{
    copy("h2_0", "h2nw");

    //diagram_1h0p_D1_const();

    //diagram_1h0p_D2a_1_const();
    diagram_1h0p_D2a_2();
    diagram_1h0p_D2b();
    //diagram_1h0p_D2c_const();
    diagram_1h0p_D2d();
    //diagram_1h0p_D2e_1_const();
    diagram_1h0p_D2e_2();

    diagram_1h0p_D3a();
    diagram_1h0p_D3b();
    diagram_1h0p_D3c();
    diagram_1h0p_D3d_1();
    diagram_1h0p_D3d_2();

    //diagram_1h0p_D4a_const();
    //diagram_1h0p_D4b_1_const();
    diagram_1h0p_D4b_2();

    diagram_1h0p_D5a();
    diagram_1h0p_D5b_1();
    diagram_1h0p_D5b_2();
    //diagram_1h0p_D5c_1_const();
    diagram_1h0p_D5c_2();
    diagram_1h0p_D5d_1();
    diagram_1h0p_D5d_2();
    //diagram_1h0p_D5e_1_const();
    diagram_1h0p_D5e_2();
    diagram_1h0p_D5f();
    //diagram_1h0p_D5g_1_const();
    diagram_1h0p_D5g_2();
    diagram_1h0p_D5h();

    //diagram_1h0p_D6a_const();
    diagram_1h0p_D6b();
    //diagram_1h0p_D6c_1_const();
    diagram_1h0p_D6c_2();

    diagram_1h0p_D7a();
    diagram_1h0p_D7b();
    diagram_1h0p_D7c_1();
    diagram_1h0p_D7c_2();
    diagram_1h0p_D7d();
    diagram_1h0p_D7e_1();
    diagram_1h0p_D7e_2();

    //diagram_1h0p_D8a_1_const();
    diagram_1h0p_D8a_2();
    diagram_1h0p_D8b();

    diagram_1h0p_D9();

    // Triples contribution to Doubles
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_contrib_to_doubles_1h0p();
    }
}


void diagram_1h0p_D1_const()
{
    copy("hhgp", "h2_0");
}


void diagram_1h0p_D2a_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pg", "pgr", "21");
    mult("t2c", "pgr", "r1", 1);
    reorder("r1", "r2", "2143");
    update("h2_0", 1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h0p_D2a_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    mult("h2c", "ppr", "r1", 1);
    update("h2nw", 1.0, "r1");

    restore_stack_pos(pos);
}


void diagram_1h0p_D2b()
{
    dg_stack_pos_t pos = get_stack_pos();

    // D2b
    mult("h2cr", "hh", "r1", 1);
    reorder("r1", "r2", "3412");
    perm("r2", "(12)");
    update("h2nw", -1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h0p_D2c_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    // D2c
    mult("ppgp", "t2c", "r1", 2);
    reorder("r1", "r2", "3412");
    update("h2_0", 0.5, "r2");

    restore_stack_pos(pos);
}


void diagram_1h0p_D2d()
{
    dg_stack_pos_t pos = get_stack_pos();

    // D2d
    mult("hhhh", "h2cr", "r1", 2);
    update("h2nw", 0.5, "r1");

    restore_stack_pos(pos);
}


void diagram_1h0p_D2e_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("t2c", "r1", "1324");
    reorder("phhg", "r2", "2431");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1324");
    reorder("r4", "r5", "2143");
    perm("r5", "(12)");
    update("h2_0", 1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h0p_D2e_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("h2c", "r1", "1324");
    mult("r1", "phhp", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(12)");
    update("h2nw", 1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h0p_D3a()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 2);
    mult("r1", "h2cr", "r2", 2);
    update("h2nw", 0.25, "r2");

    restore_stack_pos(pos);
}


void diagram_1h0p_D3b()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("h2c", "r1", "1324");
    reorder("t2c", "r2", "2413");
    reorder("pphh", "v1", "4231");
    mult("r1", "v1", "i1", 2);
    mult("i1", "r2", "i2", 2);
    reorder("i2", "r3", "1324");
    perm("r3", "(12)");
    update("h2nw", 1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h0p_D3c()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 3);
    reorder("h2c", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    perm("r3", "(12)");
    update("h2nw", -0.5, "r3");

    restore_stack_pos(pos);
}


void diagram_1h0p_D3d_1()
{
    dg_stack_pos_t pos = get_stack_pos();

    mult("t2r", "pphh", "r1", 3);
    mult("h2c", "r1", "r2", 1);
    update("h2nw", -0.5, "r2");

    restore_stack_pos(pos);
}


void diagram_1h0p_D3d_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    mult("h2cr", "pphh", "r1", 3);
    mult("t2c", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("h2nw", -0.5, "r3");

    restore_stack_pos(pos);
}


void diagram_1h0p_D4a_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phgp", "r1", "2341");
    mult("t1c", "r1", "r2", 1);
    perm("r2", "(12)");
    update("h2_0", 1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h0p_D4b_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("hhhg", "r1", "1243");
    mult("r1", "t1r", "r2", 1);
    reorder("r2", "r3", "1243");
    reorder("r3", "r4", "2143");
    update("h2_0", -1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h0p_D4b_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("hhhp", "r1", "1243");
    mult("r1", "h1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    update("h2nw", -1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5a()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("ph", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    mult("h2cr", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    perm("r4", "(12)");
    update("h2nw", -1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5b_1()
{
    dg_stack_pos_t pos = get_stack_pos();

    mult("t1r", "ph", "r1", 1);
    mult("h2c", "r1", "r2", 1);
    update("h2nw", -1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5b_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    mult("h1cr", "ph", "r1", 1);
    mult("t2c", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("h2nw", -1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5c_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphg", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    reorder("t2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    reorder("r5", "r6", "2143");
    perm("r6", "(12)");
    update("h2_0", 1.0, "r6");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5c_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphp", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    reorder("h2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12)");
    update("h2nw", 1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5d_1()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phhh", "r1", "2314");
    mult("t1r", "r1", "i1", 1);
    reorder("h2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    perm("r3", "(12)");
    update("h2nw", -1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5d_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phhh", "r1", "2314");
    mult("h1cr", "r1", "i1", 1);
    reorder("t2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    reorder("r3", "r4", "2143");
    perm("r4", "(12)");
    update("h2nw", -1.0, "r4");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5e_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphg", "r1", "4312");
    mult("t2c", "r1", "i1", 2);
    mult("i1", "t1r", "r2", 1);
    reorder("r2", "r3", "1243");
    reorder("r3", "r4", "2143");
    update("h2_0", -0.5, "r4");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5e_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphp", "r1", "4312");
    mult("t2c", "r1", "i1", 2);
    mult("i1", "h1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    update("h2nw", -0.5, "r3");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5f()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phhh", "r1", "2341");
    mult("t1c", "r1", "i1", 1);
    mult("i1", "h2cr", "r2", 2);
    perm("r2", "(12)");
    update("h2nw", 0.5, "r2");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5g_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphg", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("h2_0", 1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5g_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphp", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("h2c", "i1", "r2", 1);
    update("h2nw", 1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h0p_D5h()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("h2c", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    perm("r3", "(12)");
    update("h2nw", -1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h0p_D6a_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    mult("ppgp", "t1c", "i1", 1);
    reorder("i1", "r1", "4123");
    mult("t1c", "r1", "r2", 1);
    update("h2_0", 1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h0p_D6b()
{
    dg_stack_pos_t pos = get_stack_pos();

    mult("t1r", "hhhh", "i1", 1);
    mult("h1cr", "i1", "i2", 1);
    reorder("i2", "r1", "3412");
    update("h2nw", 1.0, "r1");

    restore_stack_pos(pos);
}


void diagram_1h0p_D6c_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phhg", "r1", "2431");
    mult("t1c", "r1", "i1", 1);
    mult("i1", "t1r", "i2", 1);
    reorder("i2", "r2", "1243");
    reorder("r2", "r3", "2143");
    perm("r3", "(12)");
    update("h2_0", -1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h0p_D6c_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    mult("t1c", "phhp", "i1", 1);
    mult("i1", "h1cr", "i2", 1);
    reorder("i2", "r2", "1243");
    perm("r2", "(12)");
    update("h2nw", -1.0, "r2");

    restore_stack_pos(pos);
}


void diagram_1h0p_D7a()
{
    dg_stack_pos_t pos = get_stack_pos();

    // D7a
    reorder("pphh", "v1  ", "3412");
    mult("t1c", "v1  ", "r1  ", 1);
    mult("t1c", "r1  ", "r2  ", 1);
    mult("r2  ", "h2cr", "r3  ", 2);
    update("h2nw", 0.5, "r3  ");

    restore_stack_pos(pos);
}


void diagram_1h0p_D7b()
{
    dg_stack_pos_t pos = get_stack_pos();

    // D7b
    reorder("pphh", "v1  ", "3412");
    mult("t2c", "v1  ", "r1  ", 2);
    mult("t1r", "r1  ", "r2  ", 1);
    mult("h1cr", "r2  ", "r3  ", 1);
    reorder("r3  ", "r4  ", "3412");
    update("h2nw", 0.5, "r4  ");

    restore_stack_pos(pos);
}


void diagram_1h0p_D7c_1()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v   ", "3241");
    mult("t1c", "v   ", "r1  ", 1);
    reorder("t2c", "t2x ", "2431");
    mult("t2x ", "r1  ", "r2  ", 2);
    mult("r2  ", "h1cr", "r3  ", 1);
    reorder("r3  ", "r4  ", "3142");
    perm("r4  ", "(12)");
    update("h2nw", -1.0, "r4  ");

    restore_stack_pos(pos);
}


void diagram_1h0p_D7c_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v   ", "3241");
    mult("t1c", "v   ", "r1  ", 1);
    reorder("h2c", "r0  ", "2143");
    reorder("r0  ", "t2x ", "2431");
    mult("t2x ", "r1  ", "r2  ", 2);
    mult("r2  ", "t1r", "r3  ", 1);
    reorder("r3  ", "r4", "3142");
    reorder("r4", "r5", "2143");
    perm("r5", "(12)");
    update("h2nw", -1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h0p_D7d()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v   ", "2431");
    mult("v   ", "t1c", "r1  ", 2);
    reorder("r1  ", "r2  ", "21");
    mult("t1c", "r2  ", "r3  ", 1);
    reorder("h2c", "h2r ", "2341");
    mult("r3  ", "h2r ", "r4  ", 1);
    perm("r4  ", "(12)");
    update("h2nw", -1.0, "r4  ");

    restore_stack_pos(pos);
}


void diagram_1h0p_D7e_1()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v   ", "2431");
    mult("v   ", "t1c", "r1  ", 2);
    mult("h1cr", "r1  ", "r2  ", 1);
    reorder("t2c", "t2x ", "1243");
    mult("t2x ", "r2  ", "r3  ", 1);
    reorder("r3  ", "r4  ", "1243");
    update("h2nw", -1.0, "r4  ");

    restore_stack_pos(pos);
}


void diagram_1h0p_D7e_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v   ", "2431");
    mult("v   ", "t1c", "r1  ", 2);
    mult("t1r", "r1  ", "r2  ", 1);
    reorder("h2c", "r0  ", "2143");
    reorder("r0  ", "t2x ", "1243");
    mult("t2x ", "r2  ", "r3  ", 1);
    reorder("r3  ", "r4", "1243");
    reorder("r4", "r5", "2143");
    update("h2nw", -1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h0p_D8a_1_const()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphg", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    mult("r3", "t1r", "r4", 1);
    reorder("r4", "r5", "1243");
    reorder("r5", "r6", "2143");
    update("h2_0", -1.0, "r6");

    restore_stack_pos(pos);
}


void diagram_1h0p_D8a_2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphp", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    mult("r3", "h1cr", "r4", 1);
    reorder("r4", "r5", "1243");
    update("h2nw", -1.0, "r5");

    restore_stack_pos(pos);
}


void diagram_1h0p_D8b()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("phhh", "r1", "2341");
    mult("t1c", "r1", "i1", 1);
    mult("t1r", "i1", "i2", 1);
    mult("h1cr", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update("h2nw", 1.0, "r3");

    restore_stack_pos(pos);
}


void diagram_1h0p_D9()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("pphh", "v1  ", "3412");
    mult("t1c", "v1  ", "r1  ", 1);
    mult("t1c", "r1  ", "r2  ", 1);
    mult("t1r", "r2  ", "r3  ", 1);
    mult("h1cr", "r3  ", "r4  ", 1);
    reorder("r4  ", "r5  ", "3412");
    update("h2nw", 1.0, "r5  ");

    restore_stack_pos(pos);
}


/**
 * Evaluates folded diagrams in the 1h0p sector.
 */
void construct_folded_1h0p()
{
    dg_stack_pos_t pos = get_stack_pos();

    // singles
    reorder("veff10", "r1", "21");
    mult("h1c", "r1", "r2", 1);
    update("h1nw", 1.0, "r2");
    restore_stack_pos(pos);

    // doubles
    reorder("h2c", "r1", "1243");
    reorder("veff10", "r2", "21");
    mult("r1", "r2", "i1", 1);
    reorder("i1", "r3", "1243");
    update("h2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // triples
    // this folded diagram appears in PT3 and hence is excluded in the CCSDT-n models
    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        t3_contrib_to_folded_1h0p();
    }
}
