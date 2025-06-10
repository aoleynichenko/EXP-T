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

#include "engine.h"
#include "options.h"
#include "ccutils.h"
#include "spinors.h"
#include "heff.h"
#include "io.h"

#include <stdio.h>
#include <stdlib.h>

double mp2_energy_goldstone();

double cc_energy_goldstone();

void print_cc_energy(double escf, double ecorr);

void construct_singles_0h0p_goldstone();

void construct_doubles_0h0p_goldstone();

void construct_singles_0h0p_mp3_goldstone();

void construct_doubles_0h0p_mp3_goldstone();

void flush_amplitudes_raw(char *diagram_name, char *path);

void symmetrize_doubles(char *source);


/**
 * Solves amplitude equations and find correlation energy for ground state
 * (Fock space sector 0h0p).
 * This code is based on Goldstone formalism and diagrams.
 */
int sector_0h0p_goldstone(cc_options_t *opts)
{
    printf("\n begin goldstone ccsd\n");

    // initialize T1 and T2 amplitudes (diagrams t1c and t2c)
    copy("hp", "t1c");
    diveps("t1c");
    copy("hhpp", "t2c");
    diveps("t2c");

    flush_amplitudes_raw("t2c", "../t2_mp2.dat");

    /*
     * useful intermediates
     */
    copy("pphh", "pphh_antisym");
    perm("pphh_antisym", "(34)");

    mp2_energy_goldstone();

    /*
     * calculate MP3 amplitudes
     */
    printf("\n");
    printf(" begin MP3 calculation\n");

    construct_singles_0h0p_mp3_goldstone();
    construct_doubles_0h0p_mp3_goldstone();

    copy("t1nw", "t1c");
    diveps("t1c");
    copy("t2nw", "t2c");
    diveps("t2c");

    cc_energy_goldstone();

    flush_amplitudes_raw("t2c", "../t2_mp3.dat");

    printf(" end of MP3 calculation\n");
    printf("\n");

    /*
     * main loop
     */
    int success = solve_amplitude_equations(
        0, 0,
        "t1c", "t1nw",
        "t2c", "t2nw",
        NULL, NULL, NULL,
        construct_singles_0h0p_goldstone,
        construct_doubles_0h0p_goldstone,
        NULL,
        NULL
    );
    if (success == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    /*
     * Total CC energy
     */
    double ecorr = cc_energy_goldstone();
    print_cc_energy(opts->escf, ecorr);
    write_formatted_heff_file_0h0p(opts->escf + ecorr);

    /*
     * analyze performance of the tensor train expansions - if needed
     */
#ifdef TENSOR_TRAIN
    //    tensor_decomposition_benchmark("ppppr", 1e-8);
//    tensor_decomposition_benchmark("t2c", 1e-8);
//    tensor_decomposition_benchmark("t3c", 1e-8);
#endif // TENSOR_TRAIN

    /*
     * analyze cluster amplitudes and save them to disk
     */
    print_cluster_operator_analysis(0, 0, "t1c", "t2c", NULL);
    save_cluster_amplitudes(0, 0, "t1c", "t2c", NULL, NULL);
    if (opts->do_flush_amplitudes_txt) {
        save_cluster_amplitudes_formatted(0, 0, "t1c", "t2c", NULL);
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

    flush_amplitudes_raw("t2c", "../t2_ccsd.dat");

    return EXIT_SUCCESS;
}


void flush_amplitudes_raw(char *diagram_name, char *path)
{
    assert_diagram_exists(diagram_name);
    diagram_t *dg = diagram_stack_find(diagram_name);

    int f = io_open(path, "w");
    io_write(f, dg->blocks[0]->buf, dg->blocks[0]->size * sizeof(double complex));
    io_close(f);
}


/**
 * Calculates MP2 energy using the Goldstone formalism
 * (no antisymmetrization of tensors).
 */
double mp2_energy_goldstone()
{
    double complex emp2 = 0.0 + 0.0 * I;
    double complex et1 = 0.0 + 0.0 * I;

    dg_stack_pos_t pos = get_stack_pos();
    perm("v", "(34)");
    double complex et2 = 0.5 * scalar_product("N", "N", "t2c", "v");
    restore_stack_pos(pos);

    //reorder("ph", "phr", "21");
    //et1 = scalar_product("N", "N", "phr", "t1c");
    //restore_stack_pos(pos);

    //if (cc_opts->print_level >= CC_PRINT_HIGH) {
    printf(" Contributions to the MP2 energy:\n");
    printf("   one-particle  %15.10f%15.10f\n", creal(et1), cimag(et1));
    printf("   two-particle  %15.10f%15.10f\n", creal(et2), cimag(et2));
    //}

    emp2 = et1 + et2;

    return creal(emp2);
}


/**
 * Calculates single-reference CC correlation energy on the current iteration.
 * TODO: optimize
 */
double cc_energy_goldstone()
{
    double complex ecc = 0.0 + 0.0 * I;
    double complex et1 = 0.0 + 0.0 * I;
    double complex et1_2 = 0.0 + 0.0 * I;
    //double complex et2 = 0.25 * scalar_product("N", "N", "v", "t2c");

    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "v", "3412");
    perm("v", "(34)");
    double complex et2 = 0.5 * scalar_product("N", "N", "t2c", "v");

    //printf("begin T1^2\n");
    copy("pphh", "r0");
    perm("r0", "(34)");

    reorder("r0", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    et1_2 = 0.5 * scalar_product("N", "N", "r2", "t1c");

    restore_stack_pos(pos);

    pos = get_stack_pos();
    reorder("ph", "phr", "21");
    et1 = scalar_product("N", "N", "phr", "t1c");
    /*reorder("pphh", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    et1_2 = 0.5 * scalar_product("N", "N", "r2", "t1c");*/
    restore_stack_pos(pos);

    ecc = et1 + et1_2 + et2;
    cc_opts->eref = cc_opts->escf + creal(ecc);  // save for higher sectors
    printf(" correlation energy = %.12f %.12f\n", creal(ecc), cimag(ecc));

    return creal(ecc);
}


/**
 * Singles equations: MP3 only terms
 */
void construct_singles_0h0p_mp3_goldstone()
{
    timer_new_entry("00-MP3-T1", "MP3 -- Singles equations");
    timer_start("00-MP3-T1");

    // S1
    copy("hp", "t1nw");

    dg_stack_pos_t pos = get_stack_pos();

    // S2a
    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");
    reorder("t2c_a", "r1", "1324");
    reorder("ph", "phr", "21");
    mult("r1", "phr", "r2", 2);
    update("t1nw", 1.0, "r2");
    restore_stack_pos(pos);

    // S2b
    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");

    reorder("pphp", "r1", "4321");
    reorder("t2c_a", "r2", "2143");
    mult("r2", "r1", "r3", 3);
    update("t1nw", 1.0, "r3");
    restore_stack_pos(pos);

    // S2c
    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");

    reorder("phhh", "r1", "2341");
    reorder("t2c_a", "r2", "4123");
    mult("r1", "r2", "r3", 3);
    update("t1nw", -1.0, "r3");
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

    // S3c - exchange
    reorder("hphp", "r1", "1432");
    mult("r1", "t1c", "r2", 2);
    update("t1nw", -1.0, "r2");
    restore_stack_pos(pos);

    timer_stop("00-MP3-T1");
}


/**
 * Doubles amplitude equations: MP3 terms only
 */
void construct_doubles_0h0p_mp3_goldstone()
{
    timer_new_entry("00-MP3-T2", "MP3 -- Doubles equations");
    timer_start("00-MP3-T2");

    // D1
    copy("hhpp", "t2nw");

    dg_stack_pos_t pos = get_stack_pos();

    // D2a
    reorder("pp", "ppr_", "21");
    mult("t2c", "ppr_", "r1", 1);
    //perm("r1", "(34)");
    symmetrize_doubles("r1");
    update("t2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // D2b
    reorder("t2c", "t2cr_", "3412");
    mult("t2cr_", "hh", "r1", 1);
    reorder("r1", "r2", "3412");
    //perm("r2", "(12)");
    symmetrize_doubles("r2");
    update("t2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D2c
    timer_start("mult_pppp");
    //tt_enable();
    mult("ppppr", "t2c", "$r1", 2);
    //tt_disable();
    timer_stop("mult_pppp");
    reorder("$r1", "r2", "3412");
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D2d !
    reorder("t2c", "r1", "3412");
    mult("hhhh", "r1", "r2", 2);
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D2e
    copy("t2c", "t2a");
    perm("t2a", "(34)");
    reorder("t2a", "r1", "1324");
    mult("r1", "phhp", "r3", 2);
    reorder("r3", "r4", "1324");
    symmetrize_doubles("r4");
    update("t2nw", 1.0, "r4");
    restore_stack_pos(pos);

    copy("t2c", "t2a");
    perm("t2a", "(34)");
    reorder("t2a", "r1", "1423");
    reorder("phph", "r2", "2341");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1342");
    symmetrize_doubles("r4");
    update("t2nw", -1.0, "r4");
    restore_stack_pos(pos);

    timer_stop("00-MP3-T2");
}


/**
 * Singles equations.
 */
void construct_singles_0h0p_goldstone()
{
    timer_new_entry("00-T1", "0h0p -- Singles equations (T1)");
    timer_start("00-T1");


    // S1
    copy("hp", "t1nw");

    dg_stack_pos_t pos = get_stack_pos();

    // S2a
    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");
    reorder("t2c_a", "r1", "1324");
    reorder("ph", "phr", "21");
    mult("r1", "phr", "r2", 2);
    update("t1nw", 1.0, "r2");
    restore_stack_pos(pos);

    // S2b
    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");

    reorder("pphp", "r1", "4321");
    reorder("t2c_a", "r2", "2143");
    mult("r2", "r1", "r3", 3);
    update("t1nw", 1.0, "r3");
    restore_stack_pos(pos);

    // S2c
    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");

    reorder("phhh", "r1", "2341");
    reorder("t2c_a", "r2", "4123");
    mult("r1", "r2", "r3", 3);
    update("t1nw", -1.0, "r3");
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

    // S3c - exchange
    reorder("hphp", "r1", "1432");
    mult("r1", "t1c", "r2", 2);
    update("t1nw", -1.0, "r2");
    restore_stack_pos(pos);

    // S4a
    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");

    reorder("t2c_a", "t2cr", "3412");
    mult("t2cr", "pphh", "r1", 3);
    mult("t1c", "r1", "r2", 1);
    update("t1nw", -1.0, "r2");
    restore_stack_pos(pos);

    // S4b
    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");

    reorder("t1c", "t1cr", "21");
    reorder("pphh", "vr1", "3412");
    mult("t2c_a", "vr1", "r1", 3);
    mult("r1", "t1cr", "r2", 1);
    update("t1nw", -1.0, "r2");
    restore_stack_pos(pos);

    // S4c
    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");

    copy("pphh", "pphh_a");
    perm("pphh_a", "(34)");

    reorder("t1c", "t1cr", "21");
    reorder("pphh_a", "r1", "2413");
    mult("r1", "t1cr", "r2", 2);
    reorder("r2", "r3", "21");
    reorder("t2c_a", "r4", "2413");
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

    // S5b - exchange
    reorder("pphp", "r1", "4132");
    mult("r1", "t1c", "r2", 2);
    mult("t1c", "r2", "r3", 1);
    update("t1nw", -1.0, "r3");
    restore_stack_pos(pos);

    // S5c
    copy("phhh", "phhh_a");
    perm("phhh_a", "(34)");

    reorder("t1c", "t1cr", "21");
    reorder("phhh_a", "r1", "2431");
    mult("r1", "t1c", "r2", 2);
    mult("r2", "t1cr", "r3", 1);
    update("t1nw", -1.0, "r3");
    restore_stack_pos(pos);

    // S6
    copy("pphh", "pphh_a");
    perm("pphh_a", "(34)");

    reorder("t1c", "t1cr", "21");
    reorder("pphh_a", "r1", "2413");
    mult("r1", "t1cr", "r2", 2);
    reorder("r2", "r3", "21");
    mult("t1c", "r3", "r4", 1);
    mult("r4", "t1cr", "r5", 1);
    update("t1nw", -1.0, "r5");
    restore_stack_pos(pos);

    // Triples contribution to Singles
    /*if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_contrib_to_singles_0h0p(PT_INF);
    }*/

    timer_stop("00-T1");
}


void symmetrize_doubles(char *source)
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder(source, "int_2143", "2143");
    update(source, 1.0, "int_2143");
    restore_stack_pos(pos);
}


/**
 * Doubles amplitude equations.
 */
void construct_doubles_0h0p_goldstone()
{
    timer_new_entry("00-T2", "0h0p -- Doubles equations (T2)");
    timer_start("00-T2");

    cc_energy_goldstone();

    // D1
    copy("hhpp", "t2nw");

    dg_stack_pos_t pos = get_stack_pos();

    // D2a
    reorder("pp", "ppr_", "21");
    mult("t2c", "ppr_", "r1", 1);
    //perm("r1", "(34)");
    symmetrize_doubles("r1");
    update("t2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // D2b
    reorder("t2c", "t2cr_", "3412");
    mult("t2cr_", "hh", "r1", 1);
    reorder("r1", "r2", "3412");
    //perm("r2", "(12)");
    symmetrize_doubles("r2");
    update("t2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D2c
    timer_start("mult_pppp");
    //tt_enable();
    mult("ppppr", "t2c", "$r1", 2);
    //tt_disable();
    timer_stop("mult_pppp");
    reorder("$r1", "r2", "3412");
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D2d !
    reorder("t2c", "r1", "3412");
    mult("hhhh", "r1", "r2", 2);
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D2e
    copy("t2c", "t2a");
    perm("t2a", "(34)");
    reorder("t2a", "r1", "1324");
    mult("r1", "phhp", "r3", 2);
    reorder("r3", "r4", "1324");
    symmetrize_doubles("r4");
    update("t2nw", 1.0, "r4");
    restore_stack_pos(pos);

    copy("t2c", "t2a");
    perm("t2a", "(34)");
    reorder("t2a", "r1", "1423");
    reorder("phph", "r2", "2341");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1342");
    symmetrize_doubles("r4");
    update("t2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D3a !
    reorder("t2c", "t2cr_", "3412");
    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 2);
    mult("r1", "t2cr_", "r2", 2);
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D3b
    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");
    reorder("t2c_a", "r1", "1324");
    reorder("pphh_antisym", "r2", "4231");
    mult("r1", "r2", "r3", 2);
    mult("r3", "r1", "r4", 2);
    reorder("r4", "r5", "1324");
    //perm("r5", "(12)");
    update("t2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // D3c
    reorder("pphh_antisym", "v1", "3412");
    mult("t2c", "v1", "r1", 3);
    reorder("t2c", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    symmetrize_doubles("r3");
    update("t2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D3d
    reorder("t2c", "t2cr_", "3412");
    mult("t2cr_", "pphh_antisym", "r1", 3);
    mult("t2c", "r1", "r2", 1);
    symmetrize_doubles("r2");
    update("t2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D4a
    reorder("phpp", "r1", "2341");
    mult("t1c", "r1", "r2", 1);
    symmetrize_doubles("r2");
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D4b
    reorder("hhhp", "r1", "1243");
    reorder("t1c", "t1cr", "21");
    mult("r1", "t1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    symmetrize_doubles("r3");
    update("t2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D5a
    reorder("ph", "r1", "21");
    reorder("t2c", "t2cr", "3412");
    mult("t1c", "r1", "r2", 1);
    mult("t2cr", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    symmetrize_doubles("r4");
    update("t2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D5b
    reorder("t1c", "t1cr", "21");
    mult("t1cr", "ph", "r1", 1);
    mult("t2c", "r1", "r2", 1);
    symmetrize_doubles("r2");
    update("t2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D5c
    copy("pphp", "pphp_a");
    perm("pphp_a", "(12)");

    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");

    reorder("pphp_a", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    reorder("t2c_a", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    symmetrize_doubles("r5");
    update("t2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // D5d
    copy("phhh", "phhh_a");
    perm("phhh_a", "(34)");

    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");

    reorder("phhh_a", "r1", "2314");
    reorder("t1c", "t1cr", "21");
    mult("t1cr", "r1", "i1", 1);
    reorder("t2c_a", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    symmetrize_doubles("r3");
    update("t2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D5e
    reorder("pphp", "r1", "4312");
    reorder("t1c", "t1cr", "21");
    mult("t2c", "r1", "i1", 2);
    mult("i1", "t1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    symmetrize_doubles("r3");
    update("t2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D5f
    reorder("phhh", "r1", "2341");
    reorder("t2c", "t2cr", "3412");
    mult("t1c", "r1", "i1", 1);
    mult("i1", "t2cr", "r2", 2);
    symmetrize_doubles("r2");
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D5g
    copy("pphp", "pphp_a");
    perm("pphp_a", "(12)");

    reorder("pphp_a", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("t2c", "i1", "r2", 1);
    symmetrize_doubles("r2");
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D5h
    copy("phhh", "phhh_a");
    perm("phhh_a", "(34)");

    reorder("phhh_a", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("t2c", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    symmetrize_doubles("r3");
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
    symmetrize_doubles("r4");
    update("t2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D6c - exchange
    reorder("t1c", "r2", "21");
    mult("t1c", "phhp", "r3", 1);
    mult("r3", "r2", "r4", 1);
    reorder("r4", "r5", "1243");
    symmetrize_doubles("r5");
    update("t2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // D7a
    reorder("pphh", "vr1", "3412");
    reorder("t2c", "t2cr", "3412");
    mult("t1c", "vr1", "r1", 1);
    mult("t1c", "r1", "r2", 1);
    mult("r2", "t2cr", "r3", 2);
    update("t2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // D7b
    reorder("pphh", "vr1", "3412");
    reorder("t1c", "t1cr", "21");
    mult("t2c", "vr1", "r1", 2);
    mult("t1cr", "r1", "r2", 1);
    mult("t1cr", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("t2nw", 1.0, "r4");
    restore_stack_pos(pos);

    // D7c
    copy("pphh", "pphh_a");
    perm("pphh_a", "(34)");

    copy("t2c", "t2c_a");
    perm("t2c_a", "(34)");

    reorder("pphh_a", "r1", "3142");
    reorder("t2c_a", "r2", "2413");
    reorder("t1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("t1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    reorder("r6", "r7", "1243");
    symmetrize_doubles("r7");
    update("t2nw", -1.0, "r7");
    restore_stack_pos(pos);

    // D7d
    copy("pphh", "pphh_a");
    perm("pphh_a", "(34)");

    reorder("pphh_a", "v_", "2431");
    mult("v_", "t1c", "r1", 2);
    reorder("r1", "r2", "21");
    mult("t1c", "r2", "r3", 1);
    reorder("t2c", "t2r", "2341");
    mult("r3", "t2r", "r4", 1);
    symmetrize_doubles("r4");
    update("t2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D7e
    copy("pphh", "pphh_a");
    perm("pphh_a", "(34)");

    reorder("pphh_a", "v_", "2431");
    reorder("t1c", "t1cr", "21");
    mult("v_", "t1c", "r1", 2);
    mult("t1cr", "r1", "r2", 1);
    reorder("t2c", "t2r", "1243");
    mult("t2r", "r2", "r3", 1);
    reorder("r3", "r4", "1243");
    symmetrize_doubles("r4");
    update("t2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D8a
    reorder("pphp", "r1", "4312");
    reorder("t1c", "t1cr", "21");
    mult("t1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    mult("r3", "t1cr", "r4", 1);
    reorder("r4", "r5", "1243");
    symmetrize_doubles("r5");
    update("t2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // D8b
    reorder("phhh", "r1", "2341");
    reorder("t1c", "t1cr", "21");
    mult("t1c", "r1", "i1", 1);
    mult("t1cr", "i1", "i2", 1);
    mult("t1cr", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    symmetrize_doubles("r3");
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

    /*
    // Triples contribution to Doubles
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_contrib_to_doubles_0h0p(PT_INF);
    }*/

    timer_stop("00-T2");
}
