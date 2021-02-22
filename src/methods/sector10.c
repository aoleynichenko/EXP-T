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
 * sector10.c
 * ==========
 *
 * Multireference Fock-space coupled cluster method for model spaces
 * with one hole (Fock space sector 1h0p)
 *
 * Cluster Operator: T = T1 + T2 + S1 + S2
 * (in this program S{10}_1 and S{10}_2 are denoted H1 and H2, respectively)
 * H10 =   |     H02 = |    \  /
 *         |           |     \/
 *        ===          =======
 *         V           V
 *         V           V
 *         |           |
 *
 * This code was obtained by modification of the CC(0,0) program:
 * all we need is to flip particle creation lines to turn them to valence hole
 * annihilation lines. See Kaldor, J. Comp. Chem. V. 8, P.448 (1987) for details.
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

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

void init_amplitudes_1h0p();

void const_terms_1h0p();

void calc_H1();

void calc_H2();

void calc_H3();

void folded_1h0p();

void t3_1h0p_contrib_to_folded();

void t3_1h0p_contrib_to_singles();

void t3_1h0p_contrib_to_doubles();


/*******************************************************************************
 * sector10
 *
 * Solves amplitude equations and find correlation energy for ground state
 * (Fock space sector 1h0p)
 ******************************************************************************/
int sector10(cc_options_t *opts)
{
    double diff1, diff2, diff3;
    int diffmax1_idx[CC_DIAGRAM_MAX_RANK];
    int diffmax2_idx[CC_DIAGRAM_MAX_RANK];
    int diffmax3_idx[CC_DIAGRAM_MAX_RANK];
    int max_t1_idx[2];
    int max_t2_idx[4];
    int max_t3_idx[4];
    double max_t1, max_t2, max_t3;
    int converged;
    int it;
    double t1, t2;
    int triples;

    printf("\n");
    printf("\t\t\t\t*****************\n");
    printf("\t\t\t\t** Sector 1h0p **\n");
    printf("\t\t\t\t*****************\n");
    printf("\n");

    opts->curr_sector_h = 1;
    opts->curr_sector_p = 0;

    triples = (opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;

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
    if (opts->cc_model >= CC_MODEL_CCSD_T3) {
        request_sorting("ppgh", "pphh", "0010", "1234");
    }
    perform_sorting();

    printf("\n Preparing T1 and T2 amplitudes ...\n");
    reorder("t1c", "t1r", "21");
    reorder("t2c", "t2r", "3412");
    reorder("pp", "ppr", "21");

    if (opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }

    printf(" Construction of S^(1,0)-independent contributions to the IPCCSD-equations ...\n");
    const_terms_1h0p();

    //calc_T1();
    //copy("t1nw", "t1c");
    //diveps("t1c");

    init_amplitudes_1h0p();

    printf(" Solution of amplitude equations (sector 1h0p)\t\t");
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

        if (opts->skip_sector[1][0]) {
            converged = 1;
            break;
        }

        reorder("h2c", "h2cr", "3412");
        reorder("h1c", "h1cr", "21");

        calc_H1();
        calc_H2();
#ifdef VERSION_DEVEL
        if (triples) {
            calc_H3();
        }
#endif

        closed("h1nw", "veff10");

        folded_1h0p();

        //vscc_update_spectator("h2nw", "t1c");

        diveps("h1nw");
        diveps("h2nw");
        if (triples) {
            diveps("h3nw");
        }

        if (cc_opts->cc_model == CC_MODEL_CCS) {
            clear("h2nw");
        }
        else if (cc_opts->cc_model == CC_MODEL_CCD) {
            clear("h1nw");
        }

        diffmax("h1c", "h1nw", &diff1, diffmax1_idx);
        diffmax("h2c", "h2nw", &diff2, diffmax2_idx);
        findmax("h1nw", &max_t1, max_t1_idx);
        findmax("h2nw", &max_t2, max_t2_idx);
        if (triples) {
            diffmax("h3c", "h3nw", &diff3, diffmax3_idx);
            findmax("h3nw", &max_t3, max_t3_idx);
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

        if (check_divergence("h1nw") == 1 ||
            check_divergence("h2nw") == 1) {
            printf("\n\t\t\tDIVERGED\n");
            return EXIT_FAILURE;
        }

        if (opts->diis_enabled) {
            diis_put(diis_queue, "h1nw", "h1c", "h2nw", "h2c", "h3nw", "h3c", it);
            if (it >= 2) {
                diis_truncate(diis_queue, opts->diis_dim);
                diis_extrapolate(diis_queue, "h1nw", "h2nw", "h3nw");
            }
        }

        // damping
        damping(1, 0, "h1c", "h1nw", it);
        damping(1, 0, "h2c", "h2nw", it);
        if (triples) {
            damping(1, 0, "h3c", "h3nw", it);
        }

        copy("h1nw", "h1c");
        copy("h2nw", "h2c");
        if (triples) {
            copy("h3nw", "h3c");
        }

        end_iter:

        // print time spent for iteration
        it_t2 = abs_time();
        double iter_time = it_t2 - it_t1;
        printf("%9.1f", iter_time);

        // flush amplitudes to disk if needed
        // (each opts->do_flush_iter iterations)
        if (opts->do_flush_iter && it % opts->do_flush_iter == 0) {
            diagram_write(diagram_stack_find("h1c"), "h1c.dg");
            diagram_write(diagram_stack_find("h2c"), "h2c.dg");
            diagram_write(diagram_stack_find("veff10"), "veff10.dg");
            if (triples) {
                diagram_write(diagram_stack_find("h3c"), "h3c.dg");
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
    print_max_ampl("Max S{10}_1 amplitude (s{10}_ia)", "h1c");
    print_max_ampl("Max S{10}_2 amplitude (s{10}_ijab)", "h2c");
    if (triples) {
        print_max_ampl("Max S{10}_3 amplitude (s{10}_ijkabc)", "h3c");
    }
    printf("\n");

    // norms of cluster operators
    print_ampl_norm("Norm |S{10}_1|", "h1c");
    print_ampl_norm("Norm |S{10}_2|", "h2c");
    if (triples) {
        print_ampl_norm("Norm |S{10}_3|", "h3c");
    }
    printf("\n");

    delete_diis_queue(diis_queue);

    // flush amplitudes and effective interaction to disk
    diagram_write(diagram_stack_find("h1c"), "h1c.dg");
    diagram_write(diagram_stack_find("h2c"), "h2c.dg");
    if (triples) {
        diagram_write(diagram_stack_find("h3c"), "h3c.dg");
    }
    diagram_write(diagram_stack_find("veff10"), "veff10.dg");

    // construct and diagonalize effective Hamiltonian
    // analyze its eigenvectors & eigenvalues
    diag_heff(1, 0, "veff10");

    // perturbative correction to the effective interaction
    // Heff will be re-constructed and diagonalized again with corrections added
    /*if (cc_opts->cc_model == CC_MODEL_CCSD_T3) {
        t3corr_0h1p();
    }*/

    return EXIT_SUCCESS;
}


/*******************************************************************************
 * init_amplitudes_1h0p
 *
 * Constructs initial guess to H1 and H2 cluster amplitudes (1h0p sector) or
 * read them from disk (if required)
 ******************************************************************************/
void init_amplitudes_1h0p()
{
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;
    int calc_s1 = 1, calc_s2 = 1, calc_s3 = 1, calc_veff = 1;

    printf("\n Initialization of S{10}_1 and S{10}_2 amplitudes ...\n");

    if (cc_opts->reuse_1h0p) {
        printf(" Trying to read amplitudes from disk ...\n");
        // Singles
        if (diagram_read("h1c.dg") != NULL) {
            printf(" S{10}_1 amplitudes successfully read from disk\n");
            calc_s1 = 0;
        }
        else {
            printf(" S{10}_1 amplitudes will be calculated\n");
        }
        // Doubles
        if (diagram_read("h2c.dg") != NULL) {
            printf(" S{10}_2 amplitudes successfully read from disk\n");
            calc_s2 = 0;
        }
        else {
            printf(" S{10}_2 amplitudes will be calculated\n");
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
                printf(" S{10}_3 amplitudes successfully read from disk\n");
                calc_s3 = 0;
            }
            else {
                printf(" S{10}_3 amplitudes will be calculated\n");
            }
        }
    }

    // calculate initial approximation to S{10} amplitudes (if needed)
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


/*******************************************************************************
 * const_terms_1h0p
 *
 * In sector 1h0p amplitude equations there are diagrams which are independent
 * on the S{10}_1 and S{10}_2 amplitudes and can be evaluated only once -- at
 * start of the calculation. They include only f, V, T{00}_1 and T{00}_2 oper-s.
 ******************************************************************************/
void const_terms_1h0p()
{
    dg_stack_pos_t pos;
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;

    timer_new_entry("const_s10", "Constant part of 1h0p amplitudes");
    timer_start("const_s10");

    tmplt("h1_0", "hh", "01", "12", NOT_PERM_UNIQUE);
    tmplt("h2_0", "hhhp", "0010", "1234", IS_PERM_UNIQUE);

    timer_stop("const_s10");
}


/*******************************************************************************
 * calc_H1
 *
 * Constructs the next approximation to S{10}_1 amplitudes ("Singles").
 ******************************************************************************/
void calc_H1()
{
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
#ifdef VERSION_DEVEL
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_1h0p_contrib_to_singles();
    }
#endif
}


/*******************************************************************************
 * calc_H2
 *
 * Constructs the next approximation to S{10}_2 amplitudes ("Doubles").
 ******************************************************************************/
void calc_H2()
{
    copy("h2_0", "h2nw");
    dg_stack_pos_t pos = get_stack_pos();

    // D1
    copy("hhgp", "h2nw");

    // diagram D2a-1
    reorder("pg", "pgr", "21");
    mult("t2c", "pgr", "r1", 1);
    reorder("r1", "r2", "2143");
    update("h2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // diagram D2a-2
    mult("h2c", "ppr", "r1", 1);
    update("h2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // D2b
    mult("h2cr", "hh", "r1", 1);
    reorder("r1", "r2", "3412");
    perm("r2", "(12)");
    update("h2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D2c
    mult("ppgp", "t2c", "r1", 2);
    reorder("r1", "r2", "3412");
    update("h2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // D2d
    mult("hhhh", "h2cr", "r1", 2);
    update("h2nw", 0.5, "r1");
    restore_stack_pos(pos);

    // diagram D2e-1
    reorder("t2c", "r1", "1324");
    reorder("phhg", "r2", "2431");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1324");
    reorder("r4", "r5", "2143");
    perm("r5", "(12)");
    update("h2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // diagram D2e-2
    reorder("h2c", "r1", "1324");
    mult("r1", "phhp", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(12)");
    update("h2nw", 1.0, "r4");
    restore_stack_pos(pos);

    // D3a
    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 2);
    mult("r1", "h2cr", "r2", 2);
    update("h2nw", 0.25, "r2");
    restore_stack_pos(pos);

    // D3b
    reorder("h2c", "r1", "1324");
    reorder("t2c", "r2", "2413");
    reorder("pphh", "v1", "4231");
    mult("r1", "v1", "i1", 2);
    mult("i1", "r2", "i2", 2);
    reorder("i2", "r3", "1324");
    perm("r3", "(12)");
    update("h2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // D3c
    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 3);
    reorder("h2c", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    perm("r3", "(12)");
    update("h2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // diagram D3d-1
    mult("t2r", "pphh", "r1", 3);
    mult("h2c", "r1", "r2", 1);
    update("h2nw", -0.5, "r2");
    restore_stack_pos(pos);

    // diagram D3d-2
    mult("h2cr", "pphh", "r1", 3);
    mult("t2c", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("h2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // D4a
    reorder("phgp", "r1", "2341");
    mult("t1c", "r1", "r2", 1);
    perm("r2", "(12)");
    update("h2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // diagram D4b-1
    reorder("hhhg", "r1", "1243");
    mult("r1", "t1r", "r2", 1);
    reorder("r2", "r3", "1243");
    reorder("r3", "r4", "2143");
    update("h2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // diagram D4b-2
    reorder("hhhp", "r1", "1243");
    mult("r1", "h1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    update("h2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D5a
    reorder("ph", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    mult("h2cr", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    perm("r4", "(12)");
    update("h2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // diagram D5b-1
    mult("t1r", "ph", "r1", 1);
    mult("h2c", "r1", "r2", 1);
    update("h2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // diagram D5b-2
    mult("h1cr", "ph", "r1", 1);
    mult("t2c", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("h2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // diagram D5c-1
    reorder("pphg", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    reorder("t2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    reorder("r5", "r6", "2143");
    perm("r6", "(12)");
    update("h2nw", 1.0, "r6");
    restore_stack_pos(pos);
    // diagram D5c-2
    reorder("pphp", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    reorder("h2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12)");
    update("h2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // diagram D5d-1
    reorder("phhh", "r1", "2314");
    mult("t1r", "r1", "i1", 1);
    reorder("h2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    perm("r3", "(12)");
    update("h2nw", -1.0, "r3");
    restore_stack_pos(pos);
    // diagram D5d-2
    reorder("phhh", "r1", "2314");
    mult("h1cr", "r1", "i1", 1);
    reorder("t2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    reorder("r3", "r4", "2143");
    perm("r4", "(12)");
    update("h2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // diagram D5e-1
    reorder("pphg", "r1", "4312");
    mult("t2c", "r1", "i1", 2);
    mult("i1", "t1r", "r2", 1);
    reorder("r2", "r3", "1243");
    reorder("r3", "r4", "2143");
    update("h2nw", -0.5, "r4");
    restore_stack_pos(pos);
    // diagram D5e-2
    reorder("pphp", "r1", "4312");
    mult("t2c", "r1", "i1", 2);
    mult("i1", "h1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    update("h2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // D5f
    reorder("phhh", "r1", "2341");
    mult("t1c", "r1", "i1", 1);
    mult("i1", "h2cr", "r2", 2);
    perm("r2", "(12)");
    update("h2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // diagram D5g-1
    reorder("pphg", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("h2nw", 1.0, "r3");
    restore_stack_pos(pos);
    // diagram D5g-2
    reorder("pphp", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("h2c", "i1", "r2", 1);
    update("h2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D5h
    reorder("phhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("h2c", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    perm("r3", "(12)");
    update("h2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D6a
    mult("ppgp", "t1c", "i1", 1);
    reorder("i1", "r1", "4123");
    mult("t1c", "r1", "r2", 1);
    update("h2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D6b
    mult("t1r", "hhhh", "i1", 1);
    mult("h1cr", "i1", "i2", 1);
    reorder("i2", "r1", "3412");
    update("h2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // diagram D6c-1
    reorder("phhg", "r1", "2431");
    mult("t1c", "r1", "i1", 1);
    mult("i1", "t1r", "i2", 1);
    reorder("i2", "r2", "1243");
    reorder("r2", "r3", "2143");
    perm("r3", "(12)");
    update("h2nw", -1.0, "r3");
    restore_stack_pos(pos);
    // diagram D6c-2
    mult("t1c", "phhp", "i1", 1);
    mult("i1", "h1cr", "i2", 1);
    reorder("i2", "r2", "1243");
    perm("r2", "(12)");
    update("h2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D7a
    reorder("pphh", "v1  ", "3412");
    mult("t1c", "v1  ", "r1  ", 1);
    mult("t1c", "r1  ", "r2  ", 1);
    mult("r2  ", "h2cr", "r3  ", 2);
    update("h2nw", 0.5, "r3  ");
    restore_stack_pos(pos);

    // D7b
    reorder("pphh", "v1  ", "3412");
    mult("t2c", "v1  ", "r1  ", 2);
    mult("t1r", "r1  ", "r2  ", 1);
    mult("h1cr", "r2  ", "r3  ", 1);
    reorder("r3  ", "r4  ", "3412");
    update("h2nw", 0.5, "r4  ");
    restore_stack_pos(pos);

    // diagram D7c-1
    reorder("pphh", "v   ", "3241");
    mult("t1c", "v   ", "r1  ", 1);
    reorder("t2c", "t2x ", "2431");
    mult("t2x ", "r1  ", "r2  ", 2);
    mult("r2  ", "h1cr", "r3  ", 1);
    reorder("r3  ", "r4  ", "3142");
    perm("r4  ", "(12)");
    update("h2nw", -1.0, "r4  ");
    restore_stack_pos(pos);
    // diagram D7c-2
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

    // D7d
    reorder("pphh", "v   ", "2431");
    mult("v   ", "t1c", "r1  ", 2);
    reorder("r1  ", "r2  ", "21");
    mult("t1c", "r2  ", "r3  ", 1);
    reorder("h2c", "h2r ", "2341");
    mult("r3  ", "h2r ", "r4  ", 1);
    perm("r4  ", "(12)");
    update("h2nw", -1.0, "r4  ");
    restore_stack_pos(pos);

    // diagram D7e-1
    reorder("pphh", "v   ", "2431");
    mult("v   ", "t1c", "r1  ", 2);
    mult("h1cr", "r1  ", "r2  ", 1);
    reorder("t2c", "t2x ", "1243");
    mult("t2x ", "r2  ", "r3  ", 1);
    reorder("r3  ", "r4  ", "1243");
    update("h2nw", -1.0, "r4  ");
    restore_stack_pos(pos);
    // diagram D7e-2
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

    // diagram D8a-1
    reorder("pphg", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    mult("r3", "t1r", "r4", 1);
    reorder("r4", "r5", "1243");
    reorder("r5", "r6", "2143");
    update("h2nw", -1.0, "r6");
    restore_stack_pos(pos);
    // diagram D8a-2
    reorder("pphp", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    mult("r3", "h1cr", "r4", 1);
    reorder("r4", "r5", "1243");
    update("h2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // D8b
    reorder("phhh", "r1", "2341");
    mult("t1c", "r1", "i1", 1);
    mult("t1r", "i1", "i2", 1);
    mult("h1cr", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update("h2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // D9
    reorder("pphh", "v1  ", "3412");
    mult("t1c", "v1  ", "r1  ", 1);
    mult("t1c", "r1  ", "r2  ", 1);
    mult("t1r", "r2  ", "r3  ", 1);
    mult("h1cr", "r3  ", "r4  ", 1);
    reorder("r4  ", "r5  ", "3412");
    update("h2nw", 1.0, "r5  ");
    restore_stack_pos(pos);

    // Triples contribution to Doubles
#ifdef VERSION_DEVEL
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_1h0p_contrib_to_doubles();
    }
#endif
}


/*******************************************************************************
 * folded
 *
 * Evaluates folded diagrams (1 for Singles and 1 for Doubles)
 ******************************************************************************/
void folded_1h0p()
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
#ifdef VERSION_DEVEL
    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        t3_1h0p_contrib_to_folded();
    }
#endif
}
