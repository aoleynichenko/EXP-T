/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2021 The EXP-T developers.
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
 * sector01.c
 * ==========
 *
 * Multireference Fock-space coupled cluster method for model spaces
 * with one particle (Fock space sector 0h1p)
 *
 * Cluster Operator: T = T1 + T2 (+ T3) + S1 + S2 (+ S3)
 *
 * This code was obtained by modification of the CC(0,0) program:
 * all we need is to flip hole creation lines to turn them to valence particle
 * annihilation lines. See Kaldor, J. Comp. Chem. V. 8, P.448 (1987) for details.
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

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
#include "options.h"
#include "methods.h"
#include "sort.h"
#include "heff.h"
#include "spinors.h"
#include "utils.h"

void init_amplitudes_0h1p();

void const_terms_0h1p();

void calc_S1();

void calc_S2();

void calc_S3(int pt_order);

void sector_0h1p_ccsd_t3();

void sector_0h1p_ccsd_t4();

void folded_0h1p();

void t3_0h1p_contrib_to_folded(int pt_order);

void t3_0h1p_contrib_to_singles(int pt_order);

void t3_0h1p_contrib_to_doubles(int pt_order);

void t3_0h1p_const_contrib_to_triples(int pt_order);


/*******************************************************************************
 * sector01
 *
 * Solves amplitude equations and find correlation energy for ground state
 * (Fock space sector 0h1p)
 ******************************************************************************/
int sector01(cc_options_t *opts)
{
    double diff1, diff2, diff3;
    int diffmax1_idx[CC_DIAGRAM_MAX_RANK];
    int diffmax2_idx[CC_DIAGRAM_MAX_RANK];
    int diffmax3_idx[CC_DIAGRAM_MAX_RANK];
    int max_t1_idx[2];
    int max_t2_idx[4];
    int max_t3_idx[6];
    double max_t1, max_t2, max_t3;
    int converged;
    int it;
    double t1, t2;
    int triples;

    printf("\n");
    printf("\t\t\t\t*****************\n");
    printf("\t\t\t\t** Sector 0h1p **\n");
    printf("\t\t\t\t*****************\n");
    printf("\n");

    opts->curr_sector_h = 0;
    opts->curr_sector_p = 1;

    triples = (opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;

    // prepare one-electron integrals
    request_sorting("vh", "ph", "10", "12");
    request_sorting("vp", "pp", "10", "12");
    request_sorting("pv", "pp", "01", "12");
    // prepare two-electron integrals
    request_sorting("vhpp", "phpp", "1000", "1234");
    request_sorting("vhhh", "phhh", "1000", "1234");
    request_sorting("pvhp", "pphp", "0100", "1234");
    request_sorting("hvhp", "hphp", "0100", "1234");
    request_sorting("pvhh", "pphh", "0100", "1234");
    request_sorting("pvpp", "pppp", "0100", "1234");
    // special for triples
    if (opts->cc_model >= CC_MODEL_CCSD_T3) {
        request_sorting("vhph", "phph", "1000", "1234");
        request_sorting("hvhh", "hphh", "0100", "1234");
        request_sorting("vppp", "pppp", "1000", "1234");
        request_sorting("vphh", "pphh", "1000", "1234");
    }
    perform_sorting();

    printf("\n Preparing T1 and T2 amplitudes ...\n");
    reorder("t1c", "t1r", "21");
    reorder("t2c", "t2r", "3412");
    reorder("pp", "ppr", "21");

    if (opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }

    printf(" Construction of S^(0,1)-independent contributions to the FSCC-equations ...\n");
    const_terms_0h1p();

    init_amplitudes_0h1p();

    predict_intruders("s2c", 5);
    if (triples) {
        predict_intruders("s3c", 5);
    }

    printf(" Solution of amplitude equations (sector 0h1p)\t\t");
    print_asctime();
    if (triples == 0) {
        printf(" ---------------------------------------------------------------------------------------\n");
        printf(" it.       diffmax(S1)       diffmax(S2)     max(S1)     max(S2)    t,sec       mem,Gb\n");
        printf(" ---------------------------------------------------------------------------------------\n");
    }
    else {
        printf(" ---------------------------------------------------------------------------------------------------------------------\n");
        printf(" it.       diffmax(S1)       diffmax(S2)       diffmax(S3)     max(S1)     max(S2)     max(S3)    t,sec       mem,Gb\n");
        printf(" ---------------------------------------------------------------------------------------------------------------------\n");
    }
    diis_queue_t *diis_queue = new_diis_queue(1, 1, triples && opts->diis_triples);
    converged = 0;
    t1 = abs_time();
    for (it = 1; it <= opts->maxiter; it++) {
        double it_t1, it_t2;
        it_t1 = abs_time();

        if (opts->skip_sector[0][1]) {
            converged = 1;
            break;
        }

        reorder("s2c", "s2cr", "3412");

        calc_S1();
        calc_S2();
#ifdef VERSION_DEVEL
        if (triples) {
            calc_S3(PT_INF);
        }
#endif

        closed("s1nw", "veff01");

        folded_0h1p();

        diveps("s1nw");
        diveps("s2nw");
        if (triples) {
            diveps("s3nw");
        }

#ifdef VERSION_DEVEL
        apply_selections(0, 1, "s1nw");
        apply_selections(0, 1, "s2nw");
        if (triples) {
            apply_selections(0, 1, "s3nw");
        }
#endif

        if (cc_opts->cc_model == CC_MODEL_CCS) {
            clear("s2nw");
        }
        else if (cc_opts->cc_model == CC_MODEL_CCD) {
            clear("s1nw");
        }

        diffmax("s1c", "s1nw", &diff1, diffmax1_idx);
        diffmax("s2c", "s2nw", &diff2, diffmax2_idx);
        findmax("s1nw", &max_t1, max_t1_idx);
        findmax("s2nw", &max_t2, max_t2_idx);
        if (triples) {
            diffmax("s3c", "s3nw", &diff3, diffmax3_idx);
            findmax("s3nw", &max_t3, max_t3_idx);
        }

        // CCSD
        if (triples == 0) {
            printf(" %3d%18.12f%18.12f%12.6f%12.6f", it, diff1, diff2, max_t1, max_t2);
            if (fabs(diff1) < opts->conv && fabs(diff2) < opts->conv) {
                converged = 1;
                goto end_iter;
            }
        }
            // CCSDT
        else {
            printf(" %3d%18.12f%18.12f%18.12f%12.6f%12.6f%12.6f", it, diff1, diff2, diff3, max_t1, max_t2, max_t3);
            if (fabs(diff1) < opts->conv && fabs(diff2) < opts->conv && fabs(diff3) < opts->conv) {
                converged = 1;
                goto end_iter;
            }
        }

        if (check_divergence("s1nw") == 1 ||
            check_divergence("s2nw") == 1) {
            printf("\n\t\t\tDIVERGED\n");
            return EXIT_FAILURE;
        }

        if (opts->diis_enabled) {
            diis_put(diis_queue, "s1nw", "s1c", "s2nw", "s2c", "s3nw", "s3c", it);
            if (it >= 2) {
                diis_truncate(diis_queue, opts->diis_dim);
                diis_extrapolate(diis_queue, "s1nw", "s2nw", "s3nw");
            }
        }

        // damping
        damping(0, 1, "s1c", "s1nw", it);
        damping(0, 1, "s2c", "s2nw", it);
        if (triples) {
            damping(0, 1, "s3c", "s3nw", it);
        }

        copy("s1nw", "s1c");
        copy("s2nw", "s2c");
        if (triples) {
            copy("s3nw", "s3c");
        }

        end_iter:

        // print time spent for iteration
        it_t2 = abs_time();
        double iter_time = it_t2 - it_t1;
        printf("%9.1f", iter_time);

        // flush amplitudes to disk if needed
        // (each opts->do_flush_iter iterations)
        if (opts->do_flush_iter && it % opts->do_flush_iter == 0) {
            diagram_write(diagram_stack_find("s1c"), "s1c.dg");
            diagram_write(diagram_stack_find("s2c"), "s2c.dg");
            diagram_write(diagram_stack_find("veff01"), "veff01.dg");
            if (triples) {
                diagram_write(diagram_stack_find("s3c"), "s3c.dg");
            }
        }

        // print memory statistics if the print level is high enough
        double curr_usage = cc_get_current_memory_usage() / (1024.0 * 1024.0 * 1024.0);
        double peak_usage = cc_get_peak_memory_usage() / (1024.0 * 1024.0 * 1024.0);
        printf("%8.2f/%.2f\n", curr_usage, peak_usage);

        if (converged) {
            break;
        }
    }
    t2 = abs_time();
    if (triples == 0) {
        printf(" ---------------------------------------------------------------------------------------\n");
    }
    else {
        printf(" ---------------------------------------------------------------------------------------------------------------------\n");
    }
    if (converged == 0) {
        printf("\tnot converged!\n");
        return EXIT_FAILURE;
    }
    else {
        printf("\tconverged in %d iterations\n", it);
    }

    printf(" average time per iteration = %.3f sec\n", (t2 - t1) / it);

    // extract max T1 and T2 amplitudes
    printf(" (absolute values)\n");
    print_max_ampl("Max S{01}_1 amplitude (s{01}_ia)", "s1c");
    print_max_ampl("Max S{01}_2 amplitude (s{01}_ijab)", "s2c");
    if (triples) {
        print_max_ampl("Max S{01}_3 amplitude (s{01}_ijkabc)", "s3c");
    }
    printf("\n");

    // norms of cluster operators
    print_ampl_norm("Norm |S{01}_1|", "s1c");
    print_ampl_norm("Norm |S{01}_2|", "s2c");
    if (triples) {
        print_ampl_norm("Norm |S{01}_3|", "s3c");
    }
    printf("\n");

    delete_diis_queue(diis_queue);

    // flush amplitudes and effective interaction to disk
    diagram_write(diagram_stack_find("s1c"), "s1c.dg");
    diagram_write(diagram_stack_find("s2c"), "s2c.dg");
    if (triples) {
        diagram_write(diagram_stack_find("s3c"), "s3c.dg");
    }
    diagram_write(diagram_stack_find("veff01"), "veff01.dg");

    /*print_ampl_vs_denom("s1c", "s1c_eps.dat");
    print_ampl_vs_denom("s2c", "s2c_eps.dat");
    if (triples) {
        print_ampl_vs_denom("s3c", "s3c_eps.dat");
    }*/

    // construct and diagonalize effective Hamiltonian
    // analyze its eigenvectors & eigenvalues
    diag_heff(0, 1, "veff01");

    // perturbative correction to the effective interaction
    // Heff will be re-constructed and diagonalized again with corrections added
#ifdef VERSION_DEVEL
    if (cc_opts->cc_model == CC_MODEL_CCSD_T3) {
        sector_0h1p_ccsd_t3();
    }
    else if (cc_opts->cc_model == CC_MODEL_CCSD_T4) {
        sector_0h1p_ccsd_t4();
    }
#endif

    return EXIT_SUCCESS;
}


/*******************************************************************************
 * init_amplitudes_0h1p
 *
 * Constructs initial guess to S1 and S2 cluster amplitudes (0h1p sector) or
 * read them from disk (if required)
 ******************************************************************************/
void init_amplitudes_0h1p()
{
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;
    int calc_s1 = 1, calc_s2 = 1, calc_s3 = 1, calc_veff = 1;

    printf("\n Initialization of S1 and S2 amplitudes ...\n");

    if (cc_opts->reuse_0h1p) {
        printf(" Trying to read amplitudes from disk ...\n");
        if (diagram_read("s1c.dg") != NULL) {
            printf(" S{01}_1 amplitudes successfully read from disk\n");
            calc_s1 = 0;
        }
        else {
            printf(" S{01}_1 amplitudes will be calculated\n");
        }
        if (diagram_read("s2c.dg") != NULL) {
            printf(" S{01}_2 amplitudes successfully read from disk\n");
            calc_s2 = 0;
        }
        else {
            printf(" S{01}_2 amplitudes will be calculated\n");
        }
        if (diagram_read("veff01.dg") != NULL) {
            printf(" Heff{01} diagram successfully read from disk\n");
            calc_veff = 0;
        }
        else {
            printf(" Heff{01} diagram will be calculated\n");
        }
        if (triples) {
            if (diagram_read("s3c.dg") != NULL) {
                printf(" S{01}_3 amplitudes successfully read from disk\n");
                calc_s3 = 0;
            }
            else {
                printf(" S{01}_3 amplitudes will be calculated\n");
            }
        }
    }

    // calculate initial approximation to S{01} amplitudes (if needed)
    if (calc_s1) {
        copy("s1_0", "s1c");
        closed("s1c", "veff01_reserve");
        diveps("s1c");
    }
    if (calc_s2) {
        copy("s2_0", "s2c");
        diveps("s2c");
    }
    if (calc_veff) {
        copy("s1_0", "veff01_open");
        closed("veff01_open", "veff01");
    }
    if (triples && calc_s3) {
        tmplt("s3c", "phhppp", "100000", "123456", IS_PERM_UNIQUE);
    }
    printf("\n");

    if (cc_opts->cc_model == CC_MODEL_CCS) {
        clear("s2c");
    }
    else if (cc_opts->cc_model == CC_MODEL_CCD) {
        clear("s1c");
    }
}


/*******************************************************************************
 * const_terms_0h1p
 *
 * In sector 0h1p amplitude equations there are diagrams which are independent
 * on the S{01}_1 and S{01}_2 amplitudes and can be evaluated only once -- at
 * start of the calculation. They include only f, V, T{00}_1 and T{00}_2 oper-s.
 ******************************************************************************/
void const_terms_0h1p()
{
    timer_new_entry("const_s01", "Constant part of 0h1p amplitudes");
    timer_start("const_s01");

    tmplt("s1_0", "pp", "10", "12", NOT_PERM_UNIQUE);
    tmplt("s2_0", "phpp", "1000", "1234", IS_PERM_UNIQUE);
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        tmplt("s3_0", "phhppp", "100000", "123456", IS_PERM_UNIQUE);
    }
    dg_stack_pos_t pos = get_stack_pos();
    //check_unique("s2_0");
    //clear_non_unique("s2_0");

    // s1_0 -- singles
    // dgs1 -- CCS
    copy("vp", "s1_0");

    // dgs2c
    reorder("pvhh", "r1", "2341");
    reorder("t2c", "r2", "4123");
    mult("r1", "r2", "r3", 3);
    update("s1_0", -0.5, "r3");
    restore_stack_pos(pos);

    // dgs3b -- CCS, малый вклад
    mult("vh", "t1r", "r1", 1);
    update("s1_0", -1.0, "r1");
    restore_stack_pos(pos);

    // dgs3c -- CCS
    reorder("pvhp", "r1", "2431");
    mult("r1", "t1c", "r2", 2);
    update("s1_0", 1.0, "r2");
    restore_stack_pos(pos);

    // dgs5c -- CCS, малые вклады
    reorder("pvhh", "r1", "2431");
    mult("r1", "t1c", "r2", 2);
    reorder("t1c", "t1r", "21");
    mult("r2", "t1r", "r3", 1);
    update("s1_0", -1.0, "r3");
    restore_stack_pos(pos);

    // s2_0 -- doubles

    // dgd1
    copy("vhpp", "s2_0");
    //check_unique("s2_0");
    //clear_non_unique("s2_0");

    // dgd2b_2
    // PROBLEM IS HERE!
    reorder("t2c", "r1", "2341");
    mult("vh", "r1", "r2", 1);
    //perm("r2", "(12)");
    perm("r2", "(34)");
    update("s2_0", -1.0, "r2");
    restore_stack_pos(pos);

    // dgd2d
    mult("vhhh", "t2r", "r1", 2);
    update("s2_0", 0.5, "r1");
    restore_stack_pos(pos);

    // dgd2e_2
    reorder("t2c", "r1", "1324");
    reorder("pvhp", "r2", "2431");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(34)");
    reorder("r4", "r5", "2143");
    update("s2_0", 1.0, "r5");
    restore_stack_pos(pos);

    // dgd4a_2
    reorder("pvpp", "r1", "2341");
    mult("t1c", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("s2_0", 1.0, "r3");
    restore_stack_pos(pos);

    // dgd4b
    reorder("hvhp", "r1", "1243");
    mult("r1", "t1r", "r2", 1);
    reorder("r2", "r3", "1243");
    reorder("r3", "r4", "2143");
    perm("r4", "(34)");
    update("s2_0", -1.0, "r4");
    restore_stack_pos(pos);

    // dgd5d_2
    reorder("pvhh", "r1", "2314");
    mult("t1r", "r1", "i1", 1);
    reorder("t2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    reorder("r3", "r4", "2143");
    perm("r4", "(34)");
    update("s2_0", -1.0, "r4");
    restore_stack_pos(pos);

    // dgd5f_2
    reorder("pvhh", "r1", "2341");
    mult("t1c", "r1", "i1", 1);
    mult("i1", "t2r", "r2", 2);
    reorder("r2", "r3", "2143");
    update("s2_0", 0.5, "r3");
    restore_stack_pos(pos);

    // dgd5h_1
    reorder("pvhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("t2c", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    update("s2_0", -1.0, "r3");
    restore_stack_pos(pos);

    // dgd6b
    mult("t1r", "vhhh", "i1", 1);
    mult("t1r", "i1", "i2", 1);
    reorder("i2", "r1", "3412");
    update("s2_0", 1.0, "r1");
    restore_stack_pos(pos);

    // dgd6c_2
    reorder("pvhp", "r1", "2431");
    mult("t1c", "r1", "i1", 1);
    mult("i1", "t1r", "i2", 1);
    reorder("i2", "r2", "1243");
    perm("r2", "(34)");
    reorder("r2", "r3", "2143");
    update("s2_0", -1.0, "r3");
    restore_stack_pos(pos);

    // dgd8b_2
    reorder("pvhh", "r1", "2341");
    mult("t1c", "r1", "i1", 1);
    mult("t1r", "i1", "i2", 1);
    mult("t1r", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    reorder("r3", "r4", "2143");
    update("s2_0", 1.0, "r4");
    restore_stack_pos(pos);

    // triples
#ifdef VERSION_DEVEL
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_0h1p_const_contrib_to_triples(PT_INF);
    }
#endif

    timer_stop("const_s01");
}


/*******************************************************************************
 * calc_S1
 *
 * Constructs the next approximation to S{01}_1 amplitudes ("Singles").
 ******************************************************************************/
void calc_S1()
{
    copy("s1_0", "s1nw");
    dg_stack_pos_t pos = get_stack_pos();

    // T2
    // dgs2a
    reorder("s2c", "r1", "1324");
    reorder("ph", "phr", "21");
    mult("r1", "phr", "r2", 2);
    update("s1nw", 1.0, "r2");
    restore_stack_pos(pos);

    reorder("pphp", "r1", "4321");
    mult("s2c", "r1", "r3", 3);
    update("s1nw", 0.5, "r3");
    restore_stack_pos(pos);

    // T1
    // dgs3a -- CCS
    mult("s1c", "ppr", "r1", 1);
    update("s1nw", 1.0, "r1");
    restore_stack_pos(pos);

    // T1*T2
    // dgs4a
    mult("t2r", "pphh", "r1", 3);
    mult("s1c", "r1", "r2", 1);
    update("s1nw", -0.5, "r2");
    restore_stack_pos(pos);
    // dgs4b
    reorder("pphh", "v1", "3412");
    mult("s2c", "v1", "r1", 3);
    mult("r1", "t1r", "r2", 1);
    update("s1nw", -0.5, "r2");
    restore_stack_pos(pos);
    // dgs4c
    reorder("pphh", "r1", "2413");
    mult("r1", "t1r", "r2", 2);
    reorder("r2", "r3", "21");
    reorder("s2c", "s2c2", "2143");
    reorder("s2c2", "r4", "2413");
    mult("r4", "r3", "r5", 2);
    update("s1nw", 1.0, "r5");
    restore_stack_pos(pos);

    // T1^2
    // dgs5a --CCS, оч маленький вклад
    reorder("ph", "r1", "21");
    mult("s1c", "r1", "r2", 1);
    mult("r2", "t1r", "r3", 1);
    update("s1nw", -1.0, "r3");
    restore_stack_pos(pos);
    // dgs5b -- CCS
    reorder("pphp", "r1", "4231");
    mult("r1", "t1c", "r2", 2);
    mult("s1c", "r2", "r3", 1);
    update("s1nw", 1.0, "r3");
    restore_stack_pos(pos);

    // T1^3
    // dgs6, CCS, оч маленький вклад
    reorder("pphh", "r1", "2413");
    mult("r1", "t1r", "r2", 2);
    reorder("r2", "r3", "21");
    mult("s1c", "r3", "r4", 1);
    mult("r4", "t1r", "r5", 1);
    update("s1nw", -1.0, "r5");
    restore_stack_pos(pos);

    // Triples contribution to Singles; for CCSDT-1a, ... only
#ifdef VERSION_DEVEL
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_0h1p_contrib_to_singles(PT_INF);
    }
#endif
}


/*******************************************************************************
 * calc_S2
 *
 * Constructs the next approximation to S{01}_2 amplitudes ("Doubles").
 ******************************************************************************/
void calc_S2()
{
    copy("s2_0", "s2nw");
    dg_stack_pos_t pos = get_stack_pos();

    // dgd2a
    mult("s2c", "ppr", "r1", 1);
    perm("r1", "(34)");
    update("s2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // dgd2b_1
    mult("s2cr", "hh", "r1", 1);
    reorder("r1", "r2", "3412");
    update("s2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // dgd2c
    timer_start("mult_pppp");
    mult("ppppr", "s2c", "r1", 2);
    timer_stop("mult_pppp");
    reorder("r1", "r2", "3412");
    update("s2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // dgd2e_1
    reorder("s2c", "r1", "1324");
    mult("r1", "phhp", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(34)");
    update("s2nw", 1.0, "r4");
    restore_stack_pos(pos);

    // dgd3a
    reorder("pphh", "v1", "3412");
    mult("s2c", "v1", "r1", 2);
    mult("r1", "t2r", "r2", 2);
    update("s2nw", 0.25, "r2");
    restore_stack_pos(pos);

    // dgd3b
    reorder("s2c", "r1", "1324");
    reorder("t2c", "r1-2", "1324");
    reorder("pphh", "r2", "4231");
    mult("r1", "r2", "r3", 2);
    mult("r3", "r1-2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(34)");
    update("s2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // dgd3c_1
    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 3);
    reorder("s2c", "r2", "3412");
    mult("r2", "r1", "r3", 1);
    reorder("r3", "r4", "3412");
    update("s2nw", -0.5, "r4");
    restore_stack_pos(pos);

    // dgd3c_2
    reorder("pphh", "v1", "3412");
    mult("s2c", "v1", "r1", 3);
    reorder("t2c", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    update("s2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // dgd3d
    mult("t2r", "pphh", "r1", 3);
    mult("s2c", "r1", "r2", 1);
    perm("r2", "(34)");
    update("s2nw", -0.5, "r2");
    restore_stack_pos(pos);

    // dgd4a_1
    reorder("phpp", "r1", "2341");
    mult("s1c", "r1", "r2", 1);
    update("s2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // dgd6a
    timer_start("mult_pppp");
    mult("ppppr", "t1c", "i1", 1);
    timer_stop("mult_pppp");
    reorder("i1", "r1", "4123");
    mult("s1c", "r1", "r2", 1);
    update("s2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // dgd6c_1
    mult("s1c", "phhp", "i1", 1);
    mult("i1", "t1r", "i2", 1);
    reorder("i2", "r2", "1243");
    perm("r2", "(34)");
    update("s2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // dgd5a_1
    reorder("ph", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    mult("s2cr", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("s2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // dgd5a_2
    reorder("ph", "r1", "21");
    mult("s1c", "r1", "r2", 1);
    mult("t2r", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    reorder("r4", "r5", "2143");
    update("s2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // dgd5b
    mult("t1r", "ph", "r1", 1);
    mult("s2c", "r1", "r2", 1);
    perm("r2", "(34)");
    update("s2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // dgd5c_1
    reorder("pphp", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    reorder("s2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(34)");
    update("s2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // dgd5c_2
    reorder("pphp", "r1", "4312");
    mult("s1c", "r1", "r2", 1);
    reorder("t2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    reorder("r5", "r6", "2143");
    perm("r6", "(34)");
    update("s2nw", 1.0, "r6");
    restore_stack_pos(pos);

    // dgd5d_1
    reorder("phhh", "r1", "2314");
    mult("t1r", "r1", "i1", 1);
    reorder("s2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    perm("r3", "(34)");
    update("s2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // dgd5e
    reorder("pphp", "r1", "4312");
    mult("s2c", "r1", "i1", 2);
    mult("i1", "t1r", "r2", 1);
    reorder("r2", "r3", "1243");
    perm("r3", "(34)");
    update("s2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // dgd5f_1
    reorder("phhh", "r1", "2341");
    mult("s1c", "r1", "i1", 1);
    mult("i1", "t2r", "r2", 2);
    update("s2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // dgd5g
    reorder("pphp", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("s2c", "i1", "r2", 1);
    perm("r2", "(34)");
    update("s2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // dgd5h_2
    reorder("phhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("s2c", "s2-1", "2143");
    reorder("s2-1", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    reorder("r3", "r4", "2143");
    update("s2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // dgd7a
    reorder("pphh", "v", "3412");
    mult("t1c", "v", "r1", 1);
    mult("s1c", "r1", "r2", 1);
    mult("r2", "t2r", "r3", 2);
    update("s2nw", 0.5, "r3");
    restore_stack_pos(pos);

    // dgd7b
    reorder("pphh", "v", "3412");
    mult("s2c", "v", "r1", 2);
    mult("t1r", "r1", "r2", 1);
    mult("t1r", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("s2nw", 0.5, "r4");
    restore_stack_pos(pos);

    // dgd7c_1
    reorder("pphh", "v", "3241");
    mult("s1c", "v", "r1", 1);
    reorder("t2c", "x", "2431");
    mult("x", "r1", "r2", 2);
    mult("r2", "t1r", "r3", 1);
    reorder("r3", "r4", "3142");
    perm("r4", "(34)");
    update("s2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // dgd7c_2
    reorder("pphh", "v", "3142");
    mult("t1c", "v", "r2", 1);
    mult("t1r", "r2", "r3", 1);
    reorder("s2c", "r4", "1324");
    mult("r4", "r3", "r5", 2);
    reorder("r5", "r6", "1423");
    perm("r6", "(34)");
    update("s2nw", -1.0, "r6");
    restore_stack_pos(pos);

    // dgd7d_1
    reorder("pphh", "r1", "4231");
    mult("r1", "t1c", "r2", 2);
    reorder("t2c", "r3", "2341");
    mult("s1c", "r2", "r4", 1);
    mult("r4", "r3", "r5", 1);
    update("s2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // dgd7d_2
    reorder("pphh", "r1", "4231");
    mult("r1", "t1c", "r2", 2);
    mult("t1c", "r2", "r4", 1);
    reorder("s2c", "s2cx", "2143");
    reorder("s2cx", "r3", "2341");
    mult("r4", "r3", "r5", 1);
    reorder("r5", "r6", "2143");
    update("s2nw", -1.0, "r6");
    restore_stack_pos(pos);

    // dgd7e
    reorder("pphh", "v", "2431");
    mult("v", "t1c", "r1", 2);
    mult("t1r", "r1", "r2", 1);
    reorder("s2c", "s2x", "1243");
    mult("s2x", "r2", "r3", 1);
    reorder("r3", "r4", "1243");
    perm("r4", "(34)");
    update("s2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // dgd8a
    reorder("pphp", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    mult("r3", "t1r", "r4", 1);
    reorder("r4", "r5", "1243");
    perm("r5", "(34)");
    update("s2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // dgd8b_1
    reorder("phhh", "r1", "2341");
    mult("s1c", "r1", "i1", 1);
    mult("t1r", "i1", "i2", 1);
    mult("t1r", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    update("s2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // dgd9
    reorder("pphh", "v", "3412");
    mult("t1c", "v", "r1", 1);
    mult("s1c", "r1", "r2", 1);
    mult("t1r", "r2", "r3", 1);
    mult("t1r", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    update("s2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // Triples contribution to Doubles
#ifdef VERSION_DEVEL
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_0h1p_contrib_to_doubles(PT_INF);
    }
#endif
}


/*******************************************************************************
 * folded_0h1p
 *
 * Evaluates folded diagrams (1 for Singles and 1 for Doubles)
 ******************************************************************************/
void folded_0h1p()
{
    dg_stack_pos_t pos = get_stack_pos();

    /*copy("pp", "vvx");
    closed("vvx", "vv");
    add(-1.0, "vv", 1.0, "veff01", "veff01");*/
    // singles
    reorder("s1c", "s1cr", "21");
    mult("veff01", "s1cr", "r1", 1);
    update("s1nw", -1.0, "r1");
    restore_stack_pos(pos);
    // doubles
    reorder("s2c", "r1", "2341");
    mult("veff01", "r1", "r2", 1);
    update("s2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // triples
    // this folded diagram appears in PT3 and hence is excluded in the CCSDT-n models
#ifdef VERSION_DEVEL
    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        t3_0h1p_contrib_to_folded(PT_INF);
    }
#endif
}