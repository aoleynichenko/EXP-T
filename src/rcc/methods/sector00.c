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
 * Single-reference coupled cluster method (Fock space sector 0h0p).
 *
 * Symbolic names used for diagrams:
 * t1c   current approximation to T1 amplitudes
 * t2c   current approximation to T2 amplitudes
 * t3c   current approximation to T3 amplitudes
 * t1nw  next (new) approximation to T1 amplitudes
 * t2nw  next (new) approximation to T2 amplitudes
 * t3nw  next (new) approximation to T3 amplitudes
 *
 * CC diagrams are listed in the book:
 *   I. Shavitt, R. J. Bartlett,
 *   Many-Body Methods in Chemistry and Physics. MBPT and Coupled-Cluster Theory.
 *   Cambridge University Press, 2010
 * (the numbering of diagrams is also adopted from this book)
 */

#include "methods.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "cc_properies.h"
#include "ccutils.h"
#include "engine.h"
#include "datamodel.h"
#include "diis.h"
#include "heff.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"


/*
 * declarations of functions used in this file
 */
void sort_integrals_0h0p();

void construct_singles_0h0p();

void construct_doubles_0h0p();

void construct_triples_0h0p();

void sector_0h0p_perturbative_triples_correction();

double mp2_energy();

double cc_energy();

void print_cc_energy(double escf, double ecorr);

void init_amplitudes_0h0p();

void t3_contrib_to_singles_0h0p(int pt_order);

void t3_contrib_to_doubles_0h0p(int pt_order);

int sector_0h0p_lambda_equations();

void sector_0h0p_analytic_density_matrix_lambda();

void sector_0h0p_analytic_density_matrix_expectation(int pt_order);


/**
 * Solves amplitude equations and find correlation energy for ground state
 * (Fock space sector 0h0p)
 */
int sector_0h0p(cc_options_t *opts)
{
    int triples = triples_enabled();

    print_sector_banner(0, 0);
    opts->curr_sector_h = 0;
    opts->curr_sector_p = 0;

    // preliminary integral sorting
    sort_integrals_0h0p();

    // initialize T1 and T2 amplitudes (diagrams t1c and t2c)
    init_amplitudes_0h0p();

    /*
     * main loop
     */
    int success = solve_amplitude_equations(
            0, 0,
            "t1c", "t1nw",
            "t2c", "t2nw",
            triples ? "t3c" : NULL, triples ? "t3nw" : NULL, NULL,
            construct_singles_0h0p,
            construct_doubles_0h0p,
            construct_triples_0h0p,
            NULL
    );
    if (success == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    /*
     * Total CC energy
     */
    double ecorr = cc_energy();
    print_cc_energy(opts->escf, ecorr);
    write_formatted_heff_0h0p(opts->escf + ecorr);

    /*
     * analyze cluster amplitudes and save them to disk
     */
    print_cluster_operator_analysis(0, 0, "t1c", "t2c", triples ? "t3c" : NULL);
    save_cluster_amplitudes(0, 0, "t1c", "t2c", triples ? "t3c" : NULL, NULL);
    if (opts->do_flush_amplitudes_txt) {
        save_cluster_amplitudes_formatted(0, 0, "t1c", "t2c", triples ? "t3c" : NULL);
    }

    /*
     * print T1 diagnostic
     */
    {
        int nelec = get_num_electrons();
        double t1d = t1_diagnostic("t1c", nelec);
        printf(" T1 diagnostic = %.8f\n", t1d);
        printf(" More on T1 diagnostic: [T.J.Lee,P.R.Taylor,IJQC 23,199(1989)]\n\n");
    }

    /*
     * [T] and (T) perturbative corrections
     */
    if (opts->cc_model == CC_MODEL_CCSD_T3 || opts->cc_model == CC_MODEL_CCSD_T4) {
        sector_0h0p_perturbative_triples_correction();
    }

    if (triples) {
        diagram_stack_erase("t3nw");
    }

    /*
     * density matrix & expectation values of property operators
     */
    if (cc_opts->calc_density_0h0p == CC_DENSITY_MATRIX_LAMBDA) {
        int exit_code = sector_0h0p_lambda_equations();
        if (exit_code == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
        sector_0h0p_analytic_density_matrix_lambda();
    }
    else if (cc_opts->calc_density_0h0p == CC_DENSITY_MATRIX_EXPECTATION) {
        sector_0h0p_analytic_density_matrix_expectation(PT_2);
    }

    sector_0h0p_calculate_properties();

    return EXIT_SUCCESS;
}


/**
 * Sorting of integrals used in the CC(0h0p) amplitude equations.
 */
void sort_integrals_0h0p()
{
    int triples = (cc_opts->cc_model < CC_MODEL_CCSDT_1A) ? 0 : 1;

    // prepare one-electron integrals
    request_sorting("hh", "hh", "00", "12");
    request_sorting("hp", "hp", "00", "12");
    request_sorting("ph", "ph", "00", "12");
    request_sorting("pp", "pp", "00", "12");

    // prepare two-electron integrals
    request_sorting("hhpp", "hhpp", "0000", "1234");
    request_sorting("pphh", "pphh", "0000", "1234");
    request_sorting("hhhh", "hhhh", "0000", "1234");
    request_sorting("phhp", "phhp", "0000", "2431");
    request_sorting("ppppr", "pppp", "0000", "3412");
    request_sorting("pphp", "pphp", "0000", "1234");
    request_sorting("phpp", "phpp", "0000", "1234");
    request_sorting("phhh", "phhh", "0000", "1234");
    request_sorting("hhhp", "hhhp", "0000", "1234");
    request_sorting("hphh", "hphh", "0000", "1234");  // for Fock matrix only
    request_sorting("hphp", "hphp", "0000", "1234");  // for Fock matrix only
    if (triples || cc_opts->cc_model == CC_MODEL_CCSD_T3 || cc_opts->cc_model == CC_MODEL_CCSD_T4) {
        request_sorting("ppph", "ppph", "0000", "1234");
        request_sorting("hhph", "hhph", "0000", "1234");
    }

    // for lambda equations
    if (cc_opts->calc_density_0h0p == CC_DENSITY_MATRIX_LAMBDA) {
        request_sorting("hpph", "hpph", "0000", "1234");
        request_sorting("pppp", "pppp", "0000", "1234");
    }

    perform_sorting();

    // hphp -> phph
    reorder("hphp", "phph", "2143");
    set_order("phph", "1234");

    reorder("pphh", "v", "3412");

    if (cc_opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }
}


/**
 * Constructs initial guess (MP2) to cluster amplitudes
 * (or read them from disk).
 */
void init_amplitudes_0h0p()
{
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;
    int calc_t1 = 1, calc_t2 = 1, calc_t3 = 1;

    printf("\n Initial guess\n");
    printf(" -------------\n");

    if (cc_opts->reuse_amplitudes[0][0]) {
        printf(" Trying to read amplitudes from disk ...\n");
        if (diagram_read("t1c.dg") != NULL) {
            printf(" T1 amplitudes successfully read from disk\n");
            calc_t1 = 0;
        }
        else {
            printf(" T1 amplitudes will be calculated\n");
        }
        if (diagram_read("t2c.dg") != NULL) {
            printf(" T2 amplitudes successfully read from disk\n");
            calc_t2 = 0;
        }
        else {
            printf(" T2 amplitudes will be calculated\n");
        }
        if (triples) {
            if (diagram_read("t3c.dg") != NULL) {
                printf(" T3 amplitudes successfully read from disk\n");
                calc_t3 = 0;
            }
            else {
                printf(" T3 amplitudes will be calculated\n");
            }
        }
        printf("\n");
    }

    // init amplitudes if needed
    if (calc_t1) {
        copy("hp", "t1c");
        diveps("t1c");
    }
    if (calc_t2) {
        copy("hhpp", "t2c");
        diveps("t2c");
    }
    if (triples && calc_t3) {
        tmplt("t3c", "hhhppp", "000000", "123456", IS_PERM_UNIQUE);
    }

    // calculate MP2 energy
    if (calc_t1 && calc_t2) {
        double emp2 = mp2_energy();
        printf(" MP2 correlation energy = %20.12f\n", emp2);
        printf("       Total MP2 energy = %20.12f\n", cc_opts->escf + emp2);
    }
    printf("\n");

    if (cc_opts->cc_model == CC_MODEL_CCS) {
        clear("t2c");
    }
    else if (cc_opts->cc_model == CC_MODEL_CCD) {
        printf("set T1 amplitudes to zero\n");
        clear("t1c");
    }
}


/**
 * Calculates MP2 energy.
 */
double mp2_energy()
{
    double complex emp2 = 0.0 + 0.0 * I;
    double complex et1 = 0.0 + 0.0 * I;
    double complex et2 = 0.25 * scalar_product("N", "N", "v", "t2c");

    dg_stack_pos_t pos = get_stack_pos();
    reorder("ph", "phr", "21");
    et1 = scalar_product("N", "N", "phr", "t1c");
    restore_stack_pos(pos);

    //if (cc_opts->print_level >= CC_PRINT_HIGH) {
    printf(" Contributions to the MP2 energy:\n");
    printf("   one-particle  %15.10f%15.10f\n", creal(et1), cimag(et1));
    printf("   two-particle  %15.10f%15.10f\n", creal(et2), cimag(et2));
    //}

    emp2 = et1 + et2;

    return creal(emp2);
}


/**
 * Prints banner with correlation, SCF and total energy
 */
void print_cc_energy(double escf, double ecorr)
{
    char total_energy_str[16];
    char cc_model_str[32];

    get_cc_model_name(0, 0, cc_opts->cc_model, cc_model_str);
    strcpy(total_energy_str, " Total ");
    strcat(total_energy_str, cc_model_str);

    printf("\n");
    printf("          SCF reference energy = %20.12f\n", escf);
    printf(" %10s correlation energy = %20.12f\n", cc_model_str, ecorr);
    printf("  %21s energy = %20.12f\n", total_energy_str, escf + ecorr);
    printf("\n");
}


/**
 * Calculates CCSD energy on the current iteration.
 */
double cc_energy()
{
    double complex ecc = 0.0 + 0.0 * I;
    double complex et1 = 0.0 + 0.0 * I;
    double complex et1_2 = 0.0 + 0.0 * I;
    double complex et2 = 0.25 * scalar_product("N", "N", "v", "t2c");

    dg_stack_pos_t pos = get_stack_pos();
    reorder("ph", "phr", "21");
    et1 = scalar_product("N", "N", "phr", "t1c");
    reorder("pphh", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    et1_2 = 0.5 * scalar_product("N", "N", "r2", "t1c");
    restore_stack_pos(pos);

    ecc = et1 + et1_2 + et2;
    cc_opts->eref = cc_opts->escf + creal(ecc);  // save for higher sectors

    return creal(ecc);
}


/**
 * Singles equations.
 */
void construct_singles_0h0p()
{
    timer_new_entry("00-T1", "0h0p -- Singles equations (T1)");
    timer_start("00-T1");

    // S1
    copy("hp", "t1nw");

    dg_stack_pos_t pos = get_stack_pos();

    // S2a
    reorder("t2c", "r1", "1324");
    reorder("ph", "phr", "21");
    mult("r1", "phr", "r2", 2);
    update("t1nw", 1.0, "r2");
    restore_stack_pos(pos);

    // S2b
    reorder("pphp", "r1", "4321");
    reorder("t2c", "r2", "2143");
    mult("r2", "r1", "r3", 3);
    update("t1nw", 0.5, "r3");
    restore_stack_pos(pos);

    // S2c
    reorder("phhh", "r1", "2341");
    reorder("t2c", "r2", "4123");
    mult("r1", "r2", "r3", 3);
    update("t1nw", -0.5, "r3");
    restore_stack_pos(pos);

    // S3a
    reorder("pp", "ppr", "21");
    mult("t1c", "ppr", "r1", 1);
    update("t1nw", 1.0, "r1");
    restore_stack_pos(pos);

    // S3b
    reorder("t1c", "t1cr", "21");
    mult("hh", "t1cr", "r1", 1);
    update("t1nw", -1.0, "r1");
    restore_stack_pos(pos);

    // S3c
    mult("phhp", "t1c", "r2", 2);
    update("t1nw", 1.0, "r2");
    restore_stack_pos(pos);

    // S4a
    reorder("t2c", "t2cr", "3412");
    mult("t2cr", "pphh", "r1", 3);
    mult("t1c", "r1", "r2", 1);
    update("t1nw", -0.5, "r2");
    restore_stack_pos(pos);

    // S4b
    reorder("t1c", "t1cr", "21");
    reorder("pphh", "vr1", "3412");
    mult("t2c", "vr1", "r1", 3);
    mult("r1", "t1cr", "r2", 1);
    update("t1nw", -0.5, "r2");
    restore_stack_pos(pos);

    // S4c
    reorder("t1c", "t1cr", "21");
    reorder("pphh", "r1", "2413");
    mult("r1", "t1cr", "r2", 2);
    reorder("r2", "r3", "21");
    reorder("t2c", "r4", "2413");
    mult("r4", "r3", "r5", 2);
    update("t1nw", 1.0, "r5");
    restore_stack_pos(pos);

    // S5a
    reorder("t1c", "t1cr", "21");
    reorder("ph", "phr", "21");
    mult("t1c", "phr", "r2", 1);
    mult("r2", "t1cr", "r3", 1);
    update("t1nw", -1.0, "r3");
    restore_stack_pos(pos);

    // S5b
    reorder("pphp", "r1", "4231");
    mult("r1", "t1c", "r2", 2);
    mult("t1c", "r2", "r3", 1);
    update("t1nw", 1.0, "r3");
    restore_stack_pos(pos);

    // S5c
    reorder("t1c", "t1cr", "21");
    reorder("phhh", "r1", "2431");
    mult("r1", "t1c", "r2", 2);
    mult("r2", "t1cr", "r3", 1);
    update("t1nw", -1.0, "r3");
    restore_stack_pos(pos);

    // S6
    reorder("t1c", "t1cr", "21");
    reorder("pphh", "r1", "2413");
    mult("r1", "t1cr", "r2", 2);
    reorder("r2", "r3", "21");
    mult("t1c", "r3", "r4", 1);
    mult("r4", "t1cr", "r5", 1);
    update("t1nw", -1.0, "r5");
    restore_stack_pos(pos);

    // Triples contribution to Singles
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_contrib_to_singles_0h0p(PT_INF);
    }

    timer_stop("00-T1");
}


/**
 * Doubles amplitude equations.
 */
void construct_doubles_0h0p()
{
    timer_new_entry("00-T2", "0h0p -- Doubles equations (T2)");
    timer_start("00-T2");

    // D1
    copy("hhpp", "t2nw");

    dg_stack_pos_t pos = get_stack_pos();

    // D2a
    reorder("pp", "ppr_", "21");
    mult("t2c", "ppr_", "r1", 1);
    perm("r1", "(34)");
    update("t2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // D2b
    reorder("t2c", "t2cr_", "3412");
    mult("t2cr_", "hh", "r1", 1);
    reorder("r1", "r2", "3412");
    perm("r2", "(12)");
    update("t2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D2c
    timer_start("mult_pppp");
    mult("ppppr", "t2c", "$r1", 2);
    timer_stop("mult_pppp");
    reorder("$r1", "r2", "3412");
    update("t2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // D2d !
    reorder("t2c", "r1", "3412");
    mult("hhhh", "r1", "r2", 2);
    update("t2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // D2e
    reorder("t2c", "r1", "1324");
    mult("r1", "phhp", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(12|34)");
    update("t2nw", 1.0, "r4");
    restore_stack_pos(pos);

    // D3a !
    reorder("t2c", "t2cr_", "3412");
    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 2);
    mult("r1", "t2cr_", "r2", 2);
    update("t2nw", 0.25, "r2");
    restore_stack_pos(pos);

    // D3b
    reorder("t2c", "r1", "1324");
    reorder("pphh", "r2", "4231");
    mult("r1", "r2", "r3", 2);
    mult("r3", "r1", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12)");
    update("t2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // D3c
    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 3);
    reorder("t2c", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    perm("r3", "(12)");
    update("t2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // D3d
    reorder("t2c", "t2cr_", "3412");
    mult("t2cr_", "pphh", "r1", 3);
    mult("t2c", "r1", "r2", 1);
    perm("r2", "(34)");
    update("t2nw", -0.5, "r2");
    restore_stack_pos(pos);

    // D4a
    reorder("phpp", "r1", "2341");
    mult("t1c", "r1", "r2", 1);
    perm("r2", "(12)");
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D4b
    reorder("hhhp", "r1", "1243");
    reorder("t1c", "t1cr", "21");
    mult("r1", "t1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    perm("r3", "(34)");
    update("t2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D5a
    reorder("ph", "r1", "21");
    reorder("t2c", "t2cr", "3412");
    mult("t1c", "r1", "r2", 1);
    mult("t2cr", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    perm("r4", "(12)");
    update("t2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D5b
    reorder("t1c", "t1cr", "21");
    mult("t1cr", "ph", "r1", 1);
    mult("t2c", "r1", "r2", 1);
    perm("r2", "(34)");
    update("t2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D5c
    reorder("pphp", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    reorder("t2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12|34)");
    update("t2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // D5d
    reorder("phhh", "r1", "2314");
    reorder("t1c", "t1cr", "21");
    mult("t1cr", "r1", "i1", 1);
    reorder("t2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    perm("r3", "(12|34)");
    update("t2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D5e
    reorder("pphp", "r1", "4312");
    reorder("t1c", "t1cr", "21");
    mult("t2c", "r1", "i1", 2);
    mult("i1", "t1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    perm("r3", "(34)");
    update("t2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // D5f
    reorder("phhh", "r1", "2341");
    reorder("t2c", "t2cr", "3412");
    mult("t1c", "r1", "i1", 1);
    mult("i1", "t2cr", "r2", 2);
    perm("r2", "(12)");
    update("t2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // D5g
    reorder("pphp", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("t2c", "i1", "r2", 1);
    perm("r2", "(34)");
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D5h
    reorder("phhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("t2c", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    perm("r3", "(12)");
    update("t2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D6a
    timer_start("mult_pppp");
    mult("ppppr", "t1c", "i1", 1);
    timer_stop("mult_pppp");
    reorder("i1", "r1", "4123");
    mult("t1c", "r1", "r2", 1);
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D6b
    reorder("t1c", "t1cr", "21");
    mult("t1cr", "hhhh", "i1", 1);
    mult("t1cr", "i1", "i2", 1);
    reorder("i2", "r1", "3412");
    update("t2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // D6c
    reorder("phph", "r1", "2341");
    reorder("t1c", "r2", "21");
    mult("t1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    perm("r4", "(12|34)");
    update("t2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D7a
    reorder("pphh", "vr1", "3412");
    reorder("t2c", "t2cr", "3412");
    mult("t1c", "vr1", "r1", 1);
    mult("t1c", "r1", "r2", 1);
    mult("r2", "t2cr", "r3", 2);
    update("t2nw", 0.5, "r3");
    restore_stack_pos(pos);

    // D7b
    reorder("pphh", "vr1", "3412");
    reorder("t1c", "t1cr", "21");
    mult("t2c", "vr1", "r1", 2);
    mult("t1cr", "r1", "r2", 1);
    mult("t1cr", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("t2nw", 0.5, "r4");
    restore_stack_pos(pos);

    // D7c
    reorder("pphh", "r1", "3142");
    reorder("t2c", "r2", "2413");
    reorder("t1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("t1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    reorder("r6", "r7", "1243");
    perm("r7", "(12|34)");
    update("t2nw", -1.0, "r7");
    restore_stack_pos(pos);

    // D7d
    reorder("pphh", "v_", "2431");
    mult("v_", "t1c", "r1", 2);
    reorder("r1", "r2", "21");
    mult("t1c", "r2", "r3", 1);
    reorder("t2c", "t2r", "2341");
    mult("r3", "t2r", "r4", 1);
    perm("r4", "(12)");
    update("t2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D7e
    reorder("pphh", "v_", "2431");
    reorder("t1c", "t1cr", "21");
    mult("v_", "t1c", "r1", 2);
    mult("t1cr", "r1", "r2", 1);
    reorder("t2c", "t2r", "1243");
    mult("t2r", "r2", "r3", 1);
    reorder("r3", "r4", "1243");
    perm("r4", "(34)");
    update("t2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D8a
    reorder("pphp", "r1", "4312");
    reorder("t1c", "t1cr", "21");
    mult("t1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    mult("r3", "t1cr", "r4", 1);
    reorder("r4", "r5", "1243");
    perm("r5", "(34)");
    update("t2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // D8b
    reorder("phhh", "r1", "2341");
    reorder("t1c", "t1cr", "21");
    mult("t1c", "r1", "i1", 1);
    mult("t1cr", "i1", "i2", 1);
    mult("t1cr", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update("t2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // D9
    reorder("t1c", "t1cr", "21");
    reorder("pphh", "vr1", "3412");
    mult("t1c", "vr1", "r1", 1);
    mult("t1c", "r1", "r2", 1);
    mult("t1cr", "r2", "r3", 1);
    mult("t1cr", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    update("t2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // Triples contribution to Doubles
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_contrib_to_doubles_0h0p(PT_INF);
    }

    timer_stop("00-T2");
}
