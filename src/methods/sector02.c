/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2020 The EXP-T developers.
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
 * sector02.c
 * ==========
 *
 * Multireference Fock-space coupled cluster method for model spaces
 * with two particles (Fock space sector 0h2p)
 *
 * Cluster Operator: T = T1 + T2 + S1 + S2 + X2
 * X2 =    |     |
 *         ^     ^
 *         |     |
 *         =======
 *         ^     ^
 *         ^     ^
 *         |     |
 * (excluding amplitudes with four valence indices)
 *
 * This code was obtained by modification of the CC(0,0) program:
 * all we need is to flip hole creation lines to turn them to valence particle
 * annihilation lines. See Kaldor, J. Comp. Chem. V. 8, P.448 (1987) for details.
 * Only two-particle diagrams are required.
 *
 * 2019 Alexander Oleynichenko
 ******************************************************************************/

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ccutils.h"
#include "datamodel.h"
#include "diis.h"
#include "engine.h"
#include "heff.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "utils.h"

void init_amplitudes_0h2p();

void const_terms_0h2p();

void t3corr_0h2p();

void calc_X2();

void calc_X3();

void folded_0h2p();

void t3_0h2p_const_contrib_to_doubles();

void t3_0h2p_contrib_to_doubles();

void t3_0h2p_contrib_to_folded();


/*******************************************************************************
 * sector02
 *
 * Solves amplitude equations in Fock space sector 0h2p
 ******************************************************************************/
int sector02(cc_options_t *opts)
{
    double diff2, diff3;
    int diffmax2_idx[CC_DIAGRAM_MAX_RANK];
    int diffmax3_idx[CC_DIAGRAM_MAX_RANK];
    int max_t2_idx[4];
    int max_t3_idx[6];
    double max_t2, max_t3;
    int converged;
    int it;
    double t1, t2;
    int triples;

    printf("\n");
    printf("\t\t\t\t*****************\n");
    printf("\t\t\t\t** Sector 0h2p **\n");
    printf("\t\t\t\t*****************\n");
    printf("\n");

    opts->curr_sector_h = 0;
    opts->curr_sector_p = 2;

    triples = (opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;

    // sorting -- interaction diagrams with valence lines
    request_sorting("vvpp", "pppp", "1100", "1234");
    request_sorting("vvhh", "pphh", "1100", "1234");
    request_sorting("vvhp", "pphp", "1100", "1234");
    if (opts->cc_model >= CC_MODEL_CCSD_T3) {
        request_sorting("vhhp", "phhp", "1000", "1234");
    }
    perform_sorting();

    printf("\n Preparing T1, T2, S{01}_1, S{01}_2 amplitudes ...\n");
    reorder("s1c", "s1r", "21");
    reorder("s2c", "s2r", "3412");

    if (opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }

    const_terms_0h2p();

    init_amplitudes_0h2p();

    printf(" Solution of amplitude equations (sector 0h2p)\n");
    print_asctime();
    if (triples == 0) {
        printf(" ---------------------------------------------------------\n");
        printf(" it.       diffmax(S2)     max(S2)    t,sec       mem,Gb\n");
        printf(" ---------------------------------------------------------\n");
    }
    else{
        printf(" ---------------------------------------------------------------------------------------\n");
        printf(" it.       diffmax(S2)       diffmax(S3)     max(S2)     max(S3)    t,sec       mem,Gb\n");
        printf(" ---------------------------------------------------------------------------------------\n");
    }
    diis_queue_t *diis_queue = new_diis_queue(0, 1, opts->diis_triples);
    converged = 0;
    t1 = abs_time();
    for (it = 1; it <= opts->maxiter; it++){

        double it_t1, it_t2;
        it_t1 = abs_time();

        reorder("x2c", "x2cr", "3412");

        calc_X2();
#ifdef VERSION_DEVEL
        if (triples) {
            calc_X3();
        }
#endif

        closed("x2nw", "veff02");

        folded_0h2p();

        diveps("x2nw");
        if (triples) {
            diveps("x3nw");
        }

#ifdef VERSION_DEVEL
        if (opts->do_relax) {
            remove_core_correlation("x2nw");
            if (triples) {
                remove_core_correlation("x3nw");
            }
        }
#endif

        if (cc_opts->cc_model == CC_MODEL_CCS) {
            clear("x2nw");
        }
#ifdef VERSION_DEVEL
        if (cc_opts->restrict_triples && triples) {
            restrict_triples("x3nw", cc_opts->restrict_triples_e1, cc_opts->restrict_triples_e2);
        }
#endif

        diffmax("x2c", "x2nw", &diff2, diffmax2_idx);
        findmax("x2nw", &max_t2, max_t2_idx);
        if (triples) {
            diffmax("x3c", "x3nw", &diff3, diffmax3_idx);
            findmax("x3nw", &max_t3, max_t3_idx);
        }

        // CCSD
        if (triples == 0) {
            printf(" %3d%18.12f%12.6f", it, diff2, max_t2);
            if (fabs(diff2) < opts->conv) {
                converged = 1;
                goto end_iter;
            }
        }
        // CCSDT
        else {
            printf(" %3d%18.12f%18.12f%12.6f%12.6f", it, diff2, diff3, max_t2, max_t3);
            if (fabs(diff2) < opts->conv && fabs(diff3) < opts->conv) {
                converged = 1;
                goto end_iter;
            }
        }

        if (check_divergence("x2nw") == 1) {
            printf("\n\t\t\tDIVERGED\n");
            return EXIT_FAILURE;
        }

        if (opts->diis_enabled) {
            diis_put(diis_queue, "", "", "x2nw", "x2c", "x3nw", "x3c", it);
            if (it >= 2) {
                diis_truncate(diis_queue, opts->diis_dim);
                diis_extrapolate(diis_queue, "", "x2nw", "x3nw");
            }
        }

        // damping
        damping(0, 2, "x2c", "x2nw", it);
        if (triples) {
            damping(0, 2, "x3c", "x3nw", it);
        }

        copy("x2nw", "x2c");
        if (triples) {
            copy("x3nw", "x3c");
        }

        end_iter:

        // print time spent for iteration
        it_t2 = abs_time();
        double iter_time = it_t2 - it_t1;
        printf("%9.1f", iter_time);

        // flush amplitudes to disk if needed
        // (each opts->do_flush_iter iterations)
        if (opts->do_flush_iter && it % opts->do_flush_iter == 0) {
            diagram_write(diagram_stack_find("x2c"), "x2c.dg");
            diagram_write(diagram_stack_find("veff02"), "veff02.dg");
            if (triples) {
                diagram_write(diagram_stack_find("x3c"), "x3c.dg");
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
    // end of iterations

    t2 = abs_time();
    if (triples == 0) {
        printf(" ---------------------------------------------------------\n");
    }
    else{
        printf(" ---------------------------------------------------------------------------------------\n");
    }
    if (converged == 0) {
        printf("\tnot converged!\n");
        return EXIT_FAILURE;
    }
    else{
        printf("\tconverged in %d iterations\n", it);
    }

    printf(" average time per iteration = %.3f sec\n", (t2 - t1) / it);

    // extract max X2 amplitudes
    printf(" (absolute values)\n");
    print_max_ampl("Max S{02}_2 amplitude (s{02}_ijab)", "x2c");
    if (triples) {
        print_max_ampl("Max S{02}_3 amplitude (s{02}_ijkabc)", "x3c");
    }
    printf("\n");

    // norms of cluster operators
    print_ampl_norm("Norm |S{02}_2|", "x2c");
    if (triples) {
        print_ampl_norm("Norm |S{02}_3|", "x3c");
    }
    printf("\n");

    delete_diis_queue(diis_queue);

    // flush amplitudes and effective interaction to disk
    diagram_write(diagram_stack_find("x2c"), "x2c.dg");
    if (triples) {
        diagram_write(diagram_stack_find("x3c"), "x3c.dg");
    }
    diagram_write(diagram_stack_find("veff02"), "veff02.dg");

    // construct and diagonalize effective Hamiltonian
    // analyze its eigenvectors & eigenvalues
    diag_heff(0, 2, "veff01", "veff02");

    // perturbative correction to the effective interaction
    // Heff will be re-constructed and diagonalized again with corrections added
#ifdef VERSION_TRIPLES
    if (cc_opts->cc_model == CC_MODEL_CCSD_T3) {
        t3corr_0h2p();
    }
#endif

    return EXIT_SUCCESS;
}


/*******************************************************************************
 * init_amplitudes_0h2p
 *
 * Constructs initial guess to X2 cluster amplitudes (0h2p sector) or
 * read them from disk (if required)
 ******************************************************************************/
void init_amplitudes_0h2p()
{
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;
    int calc_t2 = 1, calc_t3 = 1, calc_veff = 1;

    printf("\n Initialization of S{02}_2 amplitudes ...\n");

    if (cc_opts->reuse_0h2p) {
        printf(" Trying to read amplitudes (sector 0h2p) from disk ...\n");
        if (diagram_read("x2c.dg") != NULL) {
            printf(" S{02}_2 amplitudes successfully read from disk\n");
            calc_t2 = 0;
        }
        else{
            printf(" S{02}_2 amplitudes will be calculated\n");
        }
        if (triples) {
            if (diagram_read("x3c.dg") != NULL) {
                printf(" S{02}_3 amplitudes successfully read from disk\n");
                calc_t3 = 0;
            }
            else{
                printf(" S{02}_3 amplitudes will be calculated\n");
            }
        }
    }

    // init amplitudes if needed
    if (calc_t2 && calc_veff) {
        copy("vvpp", "x2c");
        closed("x2c", "veff02");
        diveps("x2c");
        printf(" Calculated: S{02}_2 and Heff{02}\n");
    }
    else if (calc_t2) {
        copy("vvpp", "x2c");
        closed("x2c", "veff02_reserve");
        diveps("x2c");
        printf(" Calculated: S{02}_2\n");
    }
    else if (calc_veff) {
        copy("vvpp", "veff02_open");
        closed("veff02_open", "veff02");
        printf(" Calculated: Heff{02}\n");
    }
    if (triples && calc_t3) {
        copy("x3_0", "x3c");
        diveps("x3c");
        printf(" Calculated: S{02}_3\n");
    }

    if (cc_opts->cc_model == CC_MODEL_CCS) {
        clear("x2c");
    }
#ifdef VERSION_DEVEL
    if (triples && cc_opts->restrict_triples) {
        restrict_triples("x3c", cc_opts->restrict_triples_e1, cc_opts->restrict_triples_e2);
    }
#endif

    printf(" done\n\n");
}


/*******************************************************************************
 * const_terms_0h2p
 *
 * In sector 0h2p amplitude equations there are diagrams which are independent
 * on the X{01}_2 amplitudes and can be evaluated only once -- at
 * start of the calculation. They include only f, V, T{00}_1, T{00}_2, S{01}_1,
 * S{01}_2 oper-s.
 ******************************************************************************/
void const_terms_0h2p()
{
    double t1, t2;

    printf(" Construction of S^(0,2)-independent contributions to the FSCC-equations ...\n");
    timer_new_entry("const_s02", "Constant part of 0h2p amplitudes");
    timer_start("const_s02");

    t1 = abs_time();

    copy("vvpp", "x2_0");
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        tmplt("x3_0", "pphppp", "110000", "123456", IS_PERM_UNIQUE);
    }
    dg_stack_pos_t pos = get_stack_pos();
    //check_unique("x2_0");
    //clear_non_unique("x2_0");

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

#ifdef VERSION_DEVEL
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_0h2p_const_contrib_to_doubles();
    }
#endif

    t2 = abs_time();
    printf(" done in %.2f sec\n", t2-t1);

    timer_stop("const_s02");
}


/*******************************************************************************
 * calc_X2
 *
 * Constructs the next approximation to S{02}_2 amplitudes ("Doubles").
 * Only S{02}_2/3 dependent diagrams will be calculated; for independent terms
 * see const_terms_0h2p().
 ******************************************************************************/
void calc_X2()
{
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
#ifdef VERSION_DEVEL
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_0h2p_contrib_to_doubles();
    }
#endif
}


/*******************************************************************************
 * folded_0h2p
 *
 * Evaluates folded diagrams in the (0,2) sector.
 ******************************************************************************/
void folded_0h2p()
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
    /*copy("pp", "vvx");
    closed("vvx", "vv");
    add(-1.0, "vv", 1.0, "veff01", "veff01");*/

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

#ifdef VERSION_DEVEL
    if (triples) {
        t3_0h2p_contrib_to_folded();
    }
#endif
}
