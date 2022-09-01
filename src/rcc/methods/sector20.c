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

/**
 * Multireference Fock-space coupled cluster method for model spaces
 * with two holes (Fock space sector 2h0p)
 *
 * Symbolic names used for diagrams:
 * h2c   current approximation to T^{2h0p}_2 amplitudes
 * h3c   current approximation to T^{2h0p}_3 amplitudes
 * h2nw  next (new) approximation to T^{2h0p}_2 amplitudes
 * h3nw  next (new) approximation to T^{2h0p}_3 amplitudes
 *
 * For the T^{2h0p}_2 operator we exclude amplitudes with four valence indices.
 *
 * This code was obtained by modification of the CC(0,0) program:
 * all we need is to flip particle creation lines to turn them to valence hole
 * annihilation lines. See Kaldor, J. Comp. Chem. V. 8, P.448 (1987) for details.
 *
 * 2019-2022 Alexander Oleynichenko
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
#include "heff.h"
#include "options.h"
#include "sort.h"
#include "utils.h"

void sort_integrals_2h0p();

void init_amplitudes_2h0p();

void const_terms_2h0p();

void construct_doubles_2h0p();

void construct_triples_2h0p(int pt_order);

void construct_folded_2h0p();

void t3_contrib_to_folded_2h0p(int pt_order);

void t3_contrib_to_doubles_2h0p(int pt_order);


/**
 * Solves amplitude equations in Fock space sector 2h0p
 */
int sector_2h0p(cc_options_t *opts)
{
    int triples = triples_enabled();

    print_sector_banner(2, 0);
    opts->curr_sector_h = 2;
    opts->curr_sector_p = 0;

    // setup shifts for the IH-like technique by A. V. Zaitsevskii
    if (intham_imms_in_sector(2, 0)) {
        intham_imms_setup(2, 0);
    }

    sort_integrals_2h0p();
    const_terms_2h0p();
    init_amplitudes_2h0p();
    intruder_state_analysis(NULL, "g2c", triples ? "g3c" : NULL);

    /*
     * main loop
     */
    int success = solve_amplitude_equations(
            2, 0,
            NULL, NULL,
            "g2c", "g2nw",
            triples ? "g3c" : NULL, triples ? "g3nw" : NULL,
            "veff20",
            NULL,
            construct_doubles_2h0p,
            construct_triples_2h0p,
            construct_folded_2h0p
    );
    if (success == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    /*
     * analyze cluster operators and save them to disk
     */
    print_cluster_operator_analysis(2, 0, NULL, "g2c", triples ? "g3c" : NULL);
    save_cluster_amplitudes(2, 0, NULL, "g2c", triples ? "g3c" : NULL, "veff20");
    if (opts->do_flush_amplitudes_txt) {
        save_cluster_amplitudes_formatted(2, 0, NULL, "g2c", triples ? "g3c" : NULL);
    }

    /*
     * construct and diagonalize effective Hamiltonian
     * analyze its eigenvectors & eigenvalues
     */
    heff_analysis(2, 0, "veff10", "veff20");
    model_space_properties_and_natural_orbitals(2, 0);

    /*
     * Ionization potentials
     */
    double ion_pot_1 = cc_opts->ground_energy_1h0p;
    double ion_pot_2 = cc_opts->ground_energy_2h0p - cc_opts->ground_energy_1h0p;

    print_ionization_potential("Ionization potential 0h0p -> 1h0p", ion_pot_1);
    print_ionization_potential("Ionization potential 1h0p -> 2h0p", ion_pot_2);

    return EXIT_SUCCESS;
}


/**
 * Sorting of integrals used in the CC(2h0p) amplitude equations.
 */
void sort_integrals_2h0p()
{
    // sorting -- interaction diagrams with valence lines
    request_sorting("hhgg", "hhhh", "0011", "1234");
    request_sorting("ppgg", "pphh", "0011", "1234");
    request_sorting("phgg", "phhh", "0011", "1234");
    perform_sorting();

    // prepare amplitudes from previous sectors
    reorder("h1c", "h1r", "21");
    reorder("h2c", "h2r", "3412");

    if (cc_opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }
}


/**
 * Constructs initial guess to G2 cluster amplitudes (2h0p sector) or
 * read them from disk (if required)
 */
void init_amplitudes_2h0p()
{
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;
    int calc_t2 = 1, calc_t3 = 1, calc_veff = 1;

    printf("\n Initialization of T{2h0p}_2 amplitudes ...\n");

    if (cc_opts->reuse_amplitudes[2][0]) {
        printf(" Trying to read amplitudes (sector 2h0p) from disk ...\n");
        if (diagram_read("g2c.dg") != NULL) {
            printf(" T{2h0p}_2 amplitudes successfully read from disk\n");
            calc_t2 = 0;
        }
        else {
            printf(" T{2h0p}_2 amplitudes will be calculated\n");
        }
        if (triples) {
            if (diagram_read("g3c.dg") != NULL) {
                printf(" T{2h0p}_3 amplitudes successfully read from disk\n");
                calc_t3 = 0;
            }
            else {
                printf(" T{2h0p}_3 amplitudes will be calculated\n");
            }
        }
    }

    // init amplitudes if needed
    if (calc_t2 && calc_veff) {
        copy("hhgg", "g2c");
        closed("g2c", "veff20");
        diveps("g2c");
        printf(" Calculated: T{2h0p}_2 and Heff{2h0p}\n");
    }
    else if (calc_t2) {
        copy("hhgg", "g2c");
        closed("g2c", "veff20_reserve");
        diveps("g2c");
        printf(" Calculated: T{2h0p}_2\n");
    }
    else if (calc_veff) {
        copy("hhgg", "veff20_open");
        closed("veff20_open", "veff20");
        printf(" Calculated: Heff{2h0p}\n");
    }
    if (triples && calc_t3) {
        tmplt("g3c", "hhhhhp", "000110", "123456", IS_PERM_UNIQUE);
        printf(" Calculated: T{2h0p}_3\n");
    }

    if (cc_opts->cc_model == CC_MODEL_CCS) {
        clear("g2c");
    }

    printf(" done\n\n");
}


void const_terms_2h0p()
{
}


/**
 * Constructs the next approximation to T{20}_2 amplitudes ("Doubles").
 */
void construct_doubles_2h0p()
{
    // prepare amplitudes
    reorder("g2c", "g2cr", "3412");

    // D1
    copy("hhgg", "g2nw");

    dg_stack_pos_t pos = get_stack_pos();

    // D2a
    reorder("pg", "pgr_", "21");
    mult("h2c", "pgr_", "r1", 1);
    perm("r1", "(34)");
    update("g2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // D2b
    reorder("g2c", "g2cr_", "3412");
    mult("g2cr_", "hh", "r1", 1);
    reorder("r1", "r2", "3412");
    perm("r2", "(12)");
    update("g2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D2c
    reorder("ppgg", "ppggr", "3412");
    mult("ppggr", "t2c", "r1", 2);
    reorder("r1", "r2", "3412");
    update("g2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // D2d
    reorder("g2c", "r1", "3412");
    mult("hhhh", "r1", "r2", 2);
    update("g2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // D2e
    reorder("h2c", "r1", "1324");
    reorder("phhg", "r0", "2431");
    mult("r1", "r0", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(12|34)");
    update("g2nw", 1.0, "r4");
    restore_stack_pos(pos);

    // D3a
    reorder("g2c", "g2cr_", "3412");
    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 2);
    mult("r1", "g2cr_", "r2", 2);
    update("g2nw", 0.25, "r2");
    restore_stack_pos(pos);

    // D3b
    reorder("h2c", "r1", "1324");
    reorder("pphh", "r2", "4231");
    mult("r1", "r2", "r3", 2);
    mult("r3", "r1", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12)");
    update("g2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // D3c
    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 3);
    reorder("g2c", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    perm("r3", "(12)");
    update("g2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // D3d
    reorder("h2c", "h2r_", "3412");
    mult("h2r_", "pphh", "r1", 3);
    mult("h2c", "r1", "r2", 1);
    perm("r2", "(34)");
    update("g2nw", -0.5, "r2");
    restore_stack_pos(pos);

    // D4a
    reorder("phgg", "r1", "2341");
    mult("t1c", "r1", "r2", 1);
    perm("r2", "(12)");
    update("g2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D4b
    reorder("hhhg", "r1", "1243");
    reorder("h1c", "h1r_", "21");
    mult("r1", "h1r_", "r2", 1);
    reorder("r2", "r3", "1243");
    perm("r3", "(34)");
    update("g2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D5a
    reorder("ph", "r1", "21");
    reorder("g2c", "g2cr_", "3412");
    mult("t1c", "r1", "r2", 1);
    mult("g2cr_", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    perm("r4", "(12)");
    update("g2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D5b
    reorder("h1c", "h1cr_", "21");
    mult("h1cr_", "ph", "r1", 1);
    mult("h2c", "r1", "r2", 1);
    perm("r2", "(34)");
    update("g2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D5c
    reorder("pphg", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    reorder("h2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12|34)");
    update("g2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // D5d
    reorder("phhh", "r1", "2314");
    reorder("h1c", "h1cr_", "21");
    mult("h1cr_", "r1", "i1", 1);
    reorder("h2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    perm("r3", "(12|34)");
    update("g2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D5e
    reorder("pphg", "r1", "4312");
    reorder("h1c", "h1cr_", "21");
    mult("t2c", "r1", "i1", 2);
    mult("i1", "h1cr_", "r2", 1);
    reorder("r2", "r3", "1243");
    perm("r3", "(34)");
    update("g2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // D5f
    reorder("phhh", "r1", "2341");
    reorder("g2c", "g2cr_", "3412");
    mult("t1c", "r1", "i1", 1);
    mult("i1", "g2cr_", "r2", 2);
    perm("r2", "(12)");
    update("g2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // D5g
    reorder("pphg", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("h2c", "i1", "r2", 1);
    perm("r2", "(34)");
    update("g2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D5h
    reorder("phhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("g2c", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    perm("r3", "(12)");
    update("g2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D6a
    reorder("ppgg", "ppggr", "3412");
    mult("ppggr", "t1c", "i1", 1);
    reorder("i1", "r1", "4123");
    mult("t1c", "r1", "r2", 1);
    update("g2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D6b
    reorder("h1c", "h1cr_", "21");
    mult("h1cr_", "hhhh", "i1", 1);
    mult("h1cr_", "i1", "i2", 1);
    reorder("i2", "r1", "3412");
    update("g2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // D6c
    reorder("phhg", "r1", "2413");
    mult("h1r", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    reorder("r3", "r4", "1324");
    perm("r4", "(12|34)");
    update("g2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D7a
    reorder("pphh", "vr1", "3412");
    reorder("g2c", "g2cr_", "3412");
    mult("t1c", "vr1", "r1", 1);
    mult("t1c", "r1", "r2", 1);
    mult("r2", "g2cr_", "r3", 2);
    update("g2nw", 0.5, "r3");
    restore_stack_pos(pos);

    // D7b
    reorder("pphh", "vr1", "3412");
    reorder("h1c", "h1r", "21");
    mult("t2c", "vr1", "r1", 2);
    mult("h1r", "r1", "r2", 1);
    mult("h1r", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("g2nw", 0.5, "r4");
    restore_stack_pos(pos);

    // D7c
    // TODO: too many reorderings
    reorder("pphh", "v_", "3241");
    mult("t1c", "v_", "r1", 1);
    reorder("h2c", "r0", "2143");
    reorder("r0", "h2r_", "2431");
    mult("h2r_", "r1", "r2", 2);
    mult("r2", "h1r", "r3", 1);
    reorder("r3", "r4", "3142");
    reorder("r4", "r5", "2143");
    perm("r5", "(12|34)");
    update("g2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // D7d
    reorder("pphh", "v_", "2431");
    mult("v_", "t1c", "r1", 2);
    reorder("r1", "r2", "21");
    mult("t1c", "r2", "r3", 1);
    reorder("g2c", "g2r_", "2341");
    mult("r3", "g2r_", "r4", 1);
    perm("r4", "(12)");
    update("g2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D7e
    reorder("pphh", "v_", "2431");
    mult("v_", "t1c", "r1", 2);
    mult("h1r", "r1", "r2", 1);
    reorder("h2c", "r0", "2143");
    reorder("r0", "h2r_", "1243");
    mult("h2r_", "r2", "r3", 1);
    reorder("r3", "r4", "1243");
    perm("r4", "(34)");
    update("g2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D8a
    reorder("pphg", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    mult("r3", "h1r", "r4", 1);
    reorder("r4", "r5", "1243");
    perm("r5", "(34)");
    update("g2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // D8b
    reorder("phhh", "r1", "2341");
    mult("t1c", "r1", "i1", 1);
    mult("h1r", "i1", "i2", 1);
    mult("h1r", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update("g2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // D9
    reorder("pphh", "vr1", "3412");
    mult("t1c", "vr1", "r1", 1);
    mult("t1c", "r1", "r2", 1);
    mult("h1r", "r2", "r3", 1);
    mult("h1r", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    update("g2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // Triples contribution to Doubles
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_contrib_to_doubles_2h0p(PT_INF);
    }
}


/**
 * Evaluates folded diagrams in the (2,0) sector.
 */
void construct_folded_2h0p()
{
    dg_stack_pos_t pos = get_stack_pos();

    // F1
    reorder("veff20", "r1", "3412");
    mult("g2c", "r1", "r2", 2);
    update("g2nw", -0.5, "r2");
    restore_stack_pos(pos);

    // F2
    reorder("veff20", "r1", "3412");
    mult("h1c", "r1", "r2", 1);
    mult("h1c", "r2", "r3", 1);
    update("g2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // F3
    reorder("veff10", "r1", "21");
    mult("g2c", "r1", "r2", 1);
    perm("r2", "(34)");
    update("g2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // F4
    reorder("veff20", "r1", "3412");
    mult("r1", "h1c", "r2", 1);
    reorder("r2", "r2r", "3412");
    tmplt("r3", "hhhh", "0011", "1234", NOT_PERM_UNIQUE);
    // unfolding
    expand_diagram("r2r", "r3");
    perm("r3", "(12)");
    update("g2nw", 1.0, "r3");
    restore_stack_pos(pos);

    if (triples_enabled()) {
        t3_contrib_to_folded_2h0p(PT_INF);
    }
}
