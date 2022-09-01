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
 * with two particles (Fock space sector 0h2p)
 *
 * Symbolic names used for diagrams:
 * x2c   current approximation to T^{0h2p}_2 amplitudes
 * x3c   current approximation to T^{0h2p}_3 amplitudes
 * x2nw  next (new) approximation to T^{0h2p}_2 amplitudes
 * x3nw  next (new) approximation to T^{0h2p}_3 amplitudes
 *
 * For the T^{0h2p}_2 operator we exclude amplitudes with four valence indices.
 *
 * This code was obtained by modification of the CC(0,0) program:
 * all we need is to flip hole creation lines to turn them to valence particle
 * annihilation lines. See Kaldor, J. Comp. Chem. V. 8, P.448 (1987) for details.
 * Only two-particle diagrams are required.
 *
 * 2019-2022 Alexander Oleynichenko
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ccutils.h"
#include "datamodel.h"
#include "diis.h"
#include "engine.h"
#include "heff.h"
#include "methods.h"
#include "options.h"
#include "sort.h"
#include "utils.h"

void sort_integrals_0h2p();

void init_amplitudes_0h2p();

void const_terms_0h2p();

void t3corr_0h2p();

void sector_0h2p_ccsd_t3();

void construct_doubles_0h2p();

void construct_triples_0h2p(int pt_order);

void construct_folded_0h2p();

void t3_const_contrib_to_triples_0h2p(int pt_order);

void t3_contrib_to_doubles_0h2p(int pt_order);

void t3_contrib_to_folded_0h2p(int pt_order);

void sector_0h2p_ccsd_t4();


/**
 * Solves amplitude equations in Fock space sector 0h2p
 */
int sector_0h2p(cc_options_t *opts)
{
    print_sector_banner(0, 2);

    opts->curr_sector_h = 0;
    opts->curr_sector_p = 2;

    int triples = triples_enabled();

    // setup shifts for the IH-like technique by A. V. Zaitsevskii
    if (intham_imms_in_sector(0, 2)) {
        intham_imms_setup(0, 2);
    }

    sort_integrals_0h2p();
    const_terms_0h2p();
    init_amplitudes_0h2p();
    intruder_state_analysis(NULL, "x2c", triples ? "x3c" : NULL);

    /*
     * main loop
     */
    int success = solve_amplitude_equations(
            0, 2,
            NULL, NULL,
            "x2c", "x2nw",
            triples ? "x3c" : NULL, triples ? "x3nw" : NULL,
            "veff02",
            NULL,
            construct_doubles_0h2p,
            construct_triples_0h2p,
            construct_folded_0h2p
    );
    if (success == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    /*
     * analyze converged amplitudes and save them to disk
     */
    print_cluster_operator_analysis(0, 2, NULL, "x2c", triples ? "x3c" : NULL);
    save_cluster_amplitudes(0, 2, NULL, "x2c", triples ? "x3c" : NULL, "veff02");
    if (opts->do_flush_amplitudes_txt) {
        save_cluster_amplitudes_formatted(0, 2, NULL, "x2c", triples ? "x3c" : NULL);
    }

    /*
     * construct and diagonalize effective Hamiltonian,
     * analyze its eigenvectors & eigenvalues
     */
    heff_analysis(0, 2, "veff01", "veff02");
    model_space_properties_and_natural_orbitals(0, 2);

    /*
     * Ionization potentials
     */
    double ion_pot_1 = cc_opts->ground_energy_0h1p - cc_opts->ground_energy_0h2p;
    double ion_pot_2 = -cc_opts->ground_energy_0h1p;

    print_ionization_potential("Ionization potential 0h2p -> 0h1p", ion_pot_1);
    print_ionization_potential("Ionization potential 0h1p -> 0h0p", ion_pot_2);

    /*
     * perturbative correction to the effective interaction.
     * Heff will be re-constructed and diagonalized again with corrections added
     */
    if (cc_opts->cc_model == CC_MODEL_CCSD_T3) {
        sector_0h2p_ccsd_t3();
    }
    else if (cc_opts->cc_model == CC_MODEL_CCSD_T4) {
        sector_0h2p_ccsd_t4();
    }

    return EXIT_SUCCESS;
}


/**
 * Sorting of integrals used in the CC(0h2p) amplitude equations.
 */
void sort_integrals_0h2p()
{
    // sorting -- interaction diagrams with valence lines
    request_sorting("vvpp", "pppp", "1100", "1234");
    request_sorting("vvhh", "pphh", "1100", "1234");
    request_sorting("vvhp", "pphp", "1100", "1234");
    if (cc_opts->cc_model >= CC_MODEL_CCSD_T3) {
        request_sorting("vhhp", "phhp", "1000", "1234");
    }
    perform_sorting();

    // prepare cluster amplitudes from previous sectors
    reorder("s1c", "s1r", "21");
    reorder("s2c", "s2r", "3412");

    if (cc_opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }
}


/**
 * Constructs initial guess to X2 cluster amplitudes (0h2p sector) or
 * read them from disk (if required)
 */
void init_amplitudes_0h2p()
{
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;
    int calc_t2 = 1, calc_t3 = 1, calc_veff = 1;

    printf("\n Initialization of T{0h2p}_2 amplitudes ...\n");

    if (cc_opts->reuse_amplitudes[0][2]) {
        printf(" Trying to read amplitudes (sector 0h2p) from disk ...\n");
        if (diagram_read("x2c.dg") != NULL) {
            printf(" T{0h2p}_2 amplitudes successfully read from disk\n");
            calc_t2 = 0;
        }
        else {
            printf(" T{0h2p}_2 amplitudes will be calculated\n");
        }
        if (triples) {
            if (diagram_read("x3c.dg") != NULL) {
                printf(" T{0h2p}_3 amplitudes successfully read from disk\n");
                calc_t3 = 0;
            }
            else {
                printf(" T{0h2p}_3 amplitudes will be calculated\n");
            }
        }
    }

    // init amplitudes if needed
    if (calc_t2 && calc_veff) {
        copy("vvpp", "x2c");
        closed("x2c", "veff02");
        diveps("x2c");
        printf(" Calculated: T{0h2p}_2 and Heff{0h2p}\n");
    }
    else if (calc_t2) {
        copy("vvpp", "x2c");
        closed("x2c", "veff02_reserve");
        diveps("x2c");
        printf(" Calculated: T{0h2p}_2\n");
    }
    else if (calc_veff) {
        copy("vvpp", "veff02_open");
        closed("veff02_open", "veff02");
        printf(" Calculated: Heff{0h2p}\n");
    }
    if (triples && calc_t3) {
        copy("x3_0", "x3c");
        diveps("x3c");
        printf(" Calculated: T{0h2p}_3\n");
    }

    if (cc_opts->cc_model == CC_MODEL_CCS) {
        clear("x2c");
    }

    printf(" done\n\n");
}


/**
 * In sector 0h2p amplitude equations there are diagrams which are independent
 * on the X{01}_2 amplitudes and can be evaluated only once -- at
 * start of the calculation. They include only f, V, T{00}_1, T{00}_2, T{01}_1,
 * T{01}_2 oper-s.
 */
void const_terms_0h2p()
{
    timer_new_entry("const_s02", "Constant part of 0h2p amplitudes");
    timer_start("const_s02");

    printf(" Construction of T^(0h2p)-independent contributions to the FSCC-equations ...\n");

    copy("vvpp", "x2_0");
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        tmplt("x3_0", "pphppp", "110000", "123456", IS_PERM_UNIQUE);
    }
    dg_stack_pos_t pos = get_stack_pos();

    // D2b
    mult("s2r", "vh", "r1", 1);
    reorder("r1", "r2", "3412");
    perm("r2", "(12)");
    update("x2_0", -1.0, "r2");
    restore_stack_pos(pos);

    // D2d
    reorder("t2c", "r1", "3412");
    mult("vvhh", "r1", "r2", 2);
    update("x2_0", 0.5, "r2");
    restore_stack_pos(pos);

    // D2e
    reorder("s2c", "r1", "1324");
    reorder("pvhp", "pvhp_r", "2431");
    mult("r1", "pvhp_r", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(12|34)");
    update("x2_0", 1.0, "r4");
    restore_stack_pos(pos);

    // D3b
    reorder("s2c", "r1", "1324");
    reorder("s2c", "r1-2", "1324");
    reorder("pphh", "r2", "4231");
    mult("r1", "r2", "r3", 2);
    mult("r3", "r1-2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12)");
    update("x2_0", 1.0, "r5");
    restore_stack_pos(pos);

    // D3c
    reorder("pphh", "v1", "3412");
    mult("s2c", "v1", "r1", 3);
    mult("s2r", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update("x2_0", -0.5, "r3");
    restore_stack_pos(pos);

    // D4a
    reorder("pvpp", "r1", "2341");
    mult("s1c", "r1", "r2", 1);
    perm("r2", "(12)");
    update("x2_0", 1.0, "r2");
    restore_stack_pos(pos);

    // D4b -- OK
    reorder("vvhp", "r1", "1243");
    mult("r1", "t1r", "r2", 1);
    reorder("r2", "r3", "1243");
    perm("r3", "(34)");
    update("x2_0", -1.0, "r3");
    restore_stack_pos(pos);

    // D5a
    reorder("ph", "r1", "21");
    mult("s1c", "r1", "r2", 1);
    mult("s2r", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    perm("r4", "(12)");
    update("x2_0", -1.0, "r4");
    restore_stack_pos(pos);

    // D5c
    reorder("pphp", "r1", "4312");
    mult("s1c", "r1", "r2", 1);
    reorder("s2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12|34)");
    update("x2_0", 1.0, "r5");
    restore_stack_pos(pos);

    // D5d
    reorder("pvhh", "r1", "2314");
    mult("t1r", "r1", "i1", 1);
    reorder("s2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    perm("r3", "(12|34)");
    update("x2_0", -1.0, "r3");
    restore_stack_pos(pos);

    // D5f
    reorder("pvhh", "r1", "2341");
    mult("s1c", "r1", "i1", 1);
    mult("i1", "t2r", "r2", 2);
    perm("r2", "(12)");
    update("x2_0", 0.5, "r2");
    restore_stack_pos(pos);

    // D5h
    reorder("pvhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("s2c", "s2-1", "2143");
    reorder("s2-1", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    reorder("r3", "r4", "2143");
    perm("r4", "(12)");
    update("x2_0", -1.0, "r4");
    restore_stack_pos(pos);

    // D6a
    timer_start("mult_pppp");
    mult("ppppr", "s1c", "i1", 1);
    timer_stop("mult_pppp");
    reorder("i1", "r1", "4123");
    mult("s1c", "r1", "r2", 1);
    update("x2_0", 1.0, "r2");
    restore_stack_pos(pos);

    // D6b
    mult("t1r", "vvhh", "i1", 1);
    mult("t1r", "i1", "i2", 1);
    reorder("i2", "r1", "3412");
    update("x2_0", 1.0, "r1");
    restore_stack_pos(pos);

    // D6c
    reorder("pvhp", "r1", "2431");
    mult("s1c", "r1", "i1", 1);
    mult("i1", "t1r", "i2", 1);
    reorder("i2", "r2", "1243");
    perm("r2", "(12|34)");
    update("x2_0", -1.0, "r2");
    restore_stack_pos(pos);

    // D7a
    reorder("pphh", "vr1", "3412");
    mult("s1c", "vr1", "r1", 1);
    mult("s1c", "r1", "r2", 1);
    mult("r2", "t2r", "r3", 2);
    update("x2_0", 0.5, "r3");
    restore_stack_pos(pos);

    // D7c
    reorder("pphh", "v", "3142");
    mult("s1c", "v", "r2", 1);
    mult("t1r", "r2", "r3", 1);
    reorder("s2c", "r4", "1324");
    mult("r4", "r3", "r5", 2);
    reorder("r5", "r6", "1423");
    perm("r6", "(12|34)");
    update("x2_0", -1.0, "r6");
    restore_stack_pos(pos);

    // D7d
    reorder("pphh", "v_", "2431");
    mult("v_", "t1c", "r1", 2);
    reorder("r1", "r2", "21");
    mult("s1c", "r2", "r3", 1);
    reorder("s2c", "s2i0", "2143");
    reorder("s2i0", "s2i", "2341");
    mult("r3", "s2i", "r4", 1);
    perm("r4", "(12)");
    update("x2_0", -1.0, "r4");
    restore_stack_pos(pos);

    // D8a
    reorder("pphp", "r1", "4312");
    mult("s1c", "r1", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    mult("r3", "t1r", "r4", 1);
    reorder("r4", "r5", "1243");
    perm("r5", "(34)");
    update("x2_0", -1.0, "r5");
    restore_stack_pos(pos);

    // D8b
    reorder("pvhh", "r1", "2341");
    mult("s1c", "r1", "i1", 1);
    mult("t1r", "i1", "i2", 1);
    mult("t1r", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update("x2_0", 1.0, "r3");
    restore_stack_pos(pos);

    // D9
    reorder("pphh", "vr1", "3412");
    mult("s1c", "vr1", "r1", 1);
    mult("s1c", "r1", "r2", 1);
    mult("t1r", "r2", "r3", 1);
    mult("t1r", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    update("x2_0", 1.0, "r5");
    restore_stack_pos(pos);

    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_const_contrib_to_triples_0h2p(PT_INF);
    }

    timer_stop("const_s02");
}


/**
 * Constructs the next approximation to T{02}_2 amplitudes ("Doubles").
 * Only T{02}_2/3 dependent diagrams will be calculated; for independent terms
 * see const_terms_0h2p().
 */
void construct_doubles_0h2p()
{
    // prepare diagrams
    reorder("x2c", "x2cr", "3412");

    // D1
    copy("x2_0", "x2nw");

    dg_stack_pos_t pos = get_stack_pos();

    // D2a
    reorder("pp", "ppr_", "21");
    mult("x2c", "ppr_", "r1", 1);
    perm("r1", "(34)");
    update("x2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // D2c
    timer_start("mult_pppp");
    mult("ppppr", "x2c", "r1", 2);
    timer_stop("mult_pppp");
    reorder("r1", "r2", "3412");
    update("x2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // D3a
    reorder("pphh", "v1", "3412");
    mult("x2c", "v1", "r1", 2);
    mult("r1", "t2r", "r2", 2);
    update("x2nw", 0.25, "r2");
    restore_stack_pos(pos);

    // D3d
    mult("t2r", "pphh", "r1", 3);
    mult("x2c", "r1", "r2", 1);
    perm("r2", "(34)");
    update("x2nw", -0.5, "r2");
    restore_stack_pos(pos);

    // D5b
    mult("t1r", "ph", "r1", 1);
    mult("x2c", "r1", "r2", 1);
    perm("r2", "(34)");
    update("x2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D5e
    reorder("pphp", "r1", "4312");
    mult("x2c", "r1", "i1", 2);
    mult("i1", "t1r", "r2", 1);
    reorder("r2", "r3", "1243");
    perm("r3", "(34)");
    update("x2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // D5g
    reorder("pphp", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("x2c", "i1", "r2", 1);
    perm("r2", "(34)");
    update("x2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D7b
    reorder("pphh", "vr1", "3412");
    mult("x2c", "vr1", "r1", 2);
    mult("t1r", "r1", "r2", 1);
    mult("t1r", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("x2nw", 0.5, "r4");
    restore_stack_pos(pos);

    // D7e
    reorder("pphh", "v_", "2431");
    mult("v_", "t1c", "r1", 2);
    mult("t1r", "r1", "r2", 1);
    reorder("x2c", "x2i", "1243");
    mult("x2i", "r2", "r3", 1);
    reorder("r3", "r4", "1243");
    perm("r4", "(34)");
    update("x2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // Triples contribution to Doubles
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_contrib_to_doubles_0h2p(PT_INF);
    }
}


/**
 * Evaluates folded diagrams in the (0,2) sector.
 */
void construct_folded_0h2p()
{
    dg_stack_pos_t pos = get_stack_pos();
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;

    // F1
    mult("veff02", "x2cr", "r1", 2);
    update("x2nw", -0.5, "r1");
    restore_stack_pos(pos);

    // F2
    mult("s1r", "veff02", "r1", 1);
    mult("s1r", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    update("x2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // F3
    mult("x2cr", "veff01", "r1", 1);
    reorder("r1", "r2", "3412");
    perm("r2", "(12)");
    update("x2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // F4
    mult("veff02", "s1r", "r1", 1);
    tmplt("r2", "pppp", "1100", "1234", NOT_PERM_UNIQUE);
    // unfolding
    expand_diagram("r1", "r2");
    perm("r2", "(34)");
    // permutation (34)
    update("x2nw", -1.0, "r2");
    restore_stack_pos(pos);

    if (triples) {
        t3_contrib_to_folded_0h2p(PT_INF);
    }
}
