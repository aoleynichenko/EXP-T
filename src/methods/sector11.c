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
 * sector11.c
 * ==========
 *
 * Multireference Fock-space coupled cluster method for model spaces
 * with one hole and one particle (Fock space sector 1h1p)
 * Calculates energies of excitations from the reference state (EE-FS-MRCC)
 *
 * Cluster Operator: E = T1 + T2 + S1 + S2 + H1 + H2 + E1 + E2
 *     E2 =    |     |
 *             V     ^
 *             |     |
 *             =======
 *             V     ^
 *             V     ^
 *             |     |
 * (excluding amplitudes with four valence indices)
 * E1 operator is not used for calculation of excitation energies, however,
 * it is of extreme importance for properties, it is calculated separately.
 *
 * This code was obtained by modification of the CC(0,0) program.
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

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

void init_amplitudes_1h1p();

void const_terms_1h1p();

void t3corr_1h1p();

void calc_E1();

void calc_E2();

void calc_E3();

void folded_1h1p_singles();

void folded_1h1p_doubles();


/*******************************************************************************
 * sector11
 *
 * Solves amplitude equations in Fock space sector 1h1p
 ******************************************************************************/
int sector11(cc_options_t *opts)
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
    printf("\t\t\t\t** Sector 1h1p **\n");
    printf("\t\t\t\t*****************\n");
    printf("\n");

    opts->curr_sector_h = 1;
    opts->curr_sector_p = 1;

    triples = (opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;

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
    if (opts->cc_model >= CC_MODEL_CCSD_T3) {
        //request_sorting("vhhp", "phhp", "1000", "1234");
        //request_sorting("vphh", "pphh", "1000", "1234");
    }
    perform_sorting();

    printf("\n Preparing T1, T2, S{01}_1, S{01}_2, S{10}_1, S{10}_2 amplitudes ...\n");
    reorder("s1c", "s1r", "21");
    reorder("s2c", "s2r", "3412");
    reorder("h1c", "h1r", "21");
    reorder("h2c", "h2r", "3412");

    if (opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }

    printf(" Construction of S^(1,1)-independent contributions to the FSCC-equations ...\n");
    const_terms_1h1p();

    init_amplitudes_1h1p();

    printf(" Solution of amplitude equations (sector 1h1p)\t\t");
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

    diis_queue_t *diis_queue = new_diis_queue(1, 1, 0);
    converged = 0;
    t1 = abs_time();
    for (it = 1; it <= opts->maxiter; it++) {
        double it_t1, it_t2;
        it_t1 = abs_time();

        if (opts->skip_sector[1][1]) {
            converged = 1;
            break;
        }

        //reorder("x2c", "x2cr", "3412");

        calc_E1();
        calc_E2();
        if (triples) {
            calc_E3();
        }

        closed("e2nw", "veff11");

        folded_1h1p_singles();
        folded_1h1p_doubles();

        diveps("e1nw");
        diveps("e2nw");
        if (triples) {
            diveps("e3nw");
        }

        if (cc_opts->cc_model == CC_MODEL_CCS) {
            clear("e2nw");
        }
        if (cc_opts->cc_model == CC_MODEL_CCD) {
            clear("e1nw");
        }

        diffmax("e1c", "e1nw", &diff1, diffmax1_idx);
        diffmax("e2c", "e2nw", &diff2, diffmax2_idx);
        findmax("e1nw", &max_t1, max_t1_idx);
        findmax("e2nw", &max_t2, max_t2_idx);
        if (triples) {
            diffmax("e3c", "e3nw", &diff3, diffmax3_idx);
            findmax("e3nw", &max_t3, max_t3_idx);
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

        if (check_divergence("e2nw") == 1) {
            printf("\n\t\t\tDIVERGED\n");
            return EXIT_FAILURE;
        }

        if (opts->diis_enabled) {
            diis_put(diis_queue, "e1nw", "e1c", "e2nw", "e2c", "", "", it);
            if (it >= 2) {
                diis_truncate(diis_queue, opts->diis_dim);
                diis_extrapolate(diis_queue, "e1nw", "e2nw", "");
            }
        }

        // damping
        damping(1, 1, "e1c", "e1nw", it);
        damping(1, 1, "e2c", "e2nw", it);
        if (triples) {
            damping(1, 1, "e3c", "e3nw", it);
        }

        copy("e1nw", "e1c");
        copy("e2nw", "e2c");
        if (triples) {
            copy("e3nw", "e3c");
        }

        end_iter:

        // print time spent for iteration
        it_t2 = abs_time();
        double iter_time = it_t2 - it_t1;
        printf("%9.1f", iter_time);

        // flush amplitudes to disk if needed
        // (each opts->do_flush_iter iterations)
        if (opts->do_flush_iter && it % opts->do_flush_iter == 0) {
            diagram_write(diagram_stack_find("e1c"), "e1c.dg");
            diagram_write(diagram_stack_find("e2c"), "e2c.dg");
            if (triples) {
                diagram_write(diagram_stack_find("e3c"), "e3c.dg");
            }
            diagram_write(diagram_stack_find("veff11"), "veff11.dg");
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
    if (!triples) {
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

    // extract max X2 amplitudes
    printf(" (absolute values)\n");
    print_max_ampl("Max S{11}_1 amplitude (s{11}_ia)", "e1c");
    print_max_ampl("Max S{11}_2 amplitude (s{11}_ijab)", "e2c");
    if (triples) {
        print_max_ampl("Max S{11}_3 amplitude (s{11}_ijkabc)", "e3c");
    }
    printf("\n");

    delete_diis_queue(diis_queue);

    // flush amplitudes and effective interaction to disk
    diagram_write(diagram_stack_find("e1c"), "e1c.dg");
    diagram_write(diagram_stack_find("e2c"), "e2c.dg");
    if (triples) {
        diagram_write(diagram_stack_find("e3c"), "e3c.dg");
    }
    diagram_write(diagram_stack_find("veff11"), "veff11.dg");

    // construct and diagonalize effective Hamiltonian
    // analyze its eigenvectors & eigenvalues
    //clear("veff11");
    //prt("veff11");
    diag_heff(1, 1, "veff01", "veff10", "veff11");

    // perturbative correction to the effective interaction
    // Heff will be re-constructed and diagonalized again with corrections added
    if (cc_opts->cc_model == CC_MODEL_CCSD_T3) {
        t3corr_1h1p();
    }

    return EXIT_SUCCESS;
}


/*******************************************************************************
 * init_amplitudes_1h1p
 *
 * Constructs initial guess to E2 cluster amplitudes (1h1p sector) or
 * read them from disk (if required)
 ******************************************************************************/
void init_amplitudes_1h1p()
{
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;
    int calc_t1 = 1, calc_t2 = 1, calc_t3 = 1, calc_veff = 1;

    printf("\n Initialization of S{11}_2 amplitudes ...\n");

    if (cc_opts->reuse_1h1p) {
        printf(" Trying to read amplitudes and effective interaction (sector 1h1p) from disk ...\n");
        // Singles
        if (diagram_read("e1c.dg") != NULL) {
            printf(" S{11}_1 amplitudes successfully read from disk\n");
            calc_t1 = 0;
        }
        else {
            printf(" S{11}_1 amplitudes will be calculated\n");
        }
        // Doubles
        if (diagram_read("e2c.dg") != NULL) {
            printf(" S{11}_2 amplitudes successfully read from disk\n");
            calc_t2 = 0;
        }
        else {
            printf(" S{11}_2 amplitudes will be calculated\n");
        }
        // Effective interaction
        if (diagram_read("veff11.dg") != NULL) {
            printf(" Heff{11} diagram successfully read from disk\n");
            calc_veff = 0;
        }
        else {
            printf(" Heff{11} diagram will be calculated\n");
        }
        /*if (triples) {
            if (diagram_read("e3c.dg") != NULL) {
                printf(" S{11}_3 amplitudes successfully read from disk\n");
                calc_t3 = 0;
            }
            else {
                printf(" S{11}_3 amplitudes will be calculated\n");
            }
        }*/
    }

    // init amplitudes if needed
    if (calc_t1) {
        copy("vg", "e1c");
        diveps("e1c");
        printf(" Calculated: S{11}_1\n");
    }
    if (calc_t2 && calc_veff) {
        copy("vhpg", "e2c");
        closed("e2c", "veff11");
        diveps("e2c");
        printf(" Calculated: S{11}_2 and Heff{11}\n");
    }
    else if (calc_t2) {
        copy("vhpg", "e2c");
        closed("e2c", "veff11_reserve");
        diveps("e2c");
        printf(" Calculated: S{11}_2\n");
    }
    else if (calc_veff) {
        copy("vhpg", "veff11_open");
        closed("veff11_open", "veff11");
        printf(" Calculated: Heff{11}\n");
    }
    if (triples && calc_t3) {
        /*tmplt("x3c", "pphppp", "110000", "123456");
        printf(" Calculated: S{02}_3\n");*/
    }

    if (cc_opts->cc_model == CC_MODEL_CCS) {
        clear("e2c");
    }
    if (cc_opts->cc_model == CC_MODEL_CCD) {
        clear("e1c");
    }

    printf(" done\n\n");
}


/*******************************************************************************
 * const_terms_1h1p
 *
 * In sector 0h2p amplitude equations there are diagrams which are independent
 * on the X{01}_2 amplitudes and can be evaluated only once -- at
 * start of the calculation. They include only f, V, T{00}_1, T{00}_2, S{01}_1,
 * S{01}_2 oper-s.
 ******************************************************************************/
void const_terms_1h1p()
{
    timer_new_entry("const_s11", "Constant part of 1h1p amplitudes");
    timer_start("const_s11");

    copy("vhpg", "e2_0");

    timer_stop("const_s11");
}


/*******************************************************************************
 * calc_E1
 *
 * Singles equations for the 1h1p sector.
 ******************************************************************************/
void calc_E1()
{
    timer_new_entry("00-E1", "1h1p -- Singles equations (Exc1)");
    timer_start("00-E1");

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
    /*if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {

    }*/

    timer_stop("00-E1");
}


/*******************************************************************************
 * calc_X2
 *
 * Constructs the next approximation to S{02}_2 amplitudes ("Doubles").
 * Only S{02}_2/3 dependent diagrams will be calculated; for independent terms
 * see const_terms_0h2p().
 ******************************************************************************/
void calc_E2()
{
    timer_new_entry("11-S2", "1h1p -- Doubles equations (S{11}_2)");
    timer_start("11-S2");

    // D1
    copy("e2_0", "e2nw");

    dg_stack_pos_t pos = get_stack_pos();

    // D2a-1
    reorder("pg", "pgr", "21");
    mult("s2c", "pgr", "r1", 1);
    update("e2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // D2a-2
    reorder("pp", "ppr_", "21");
    reorder("e2c", "e2_21", "2143");
    mult("e2_21", "ppr_", "r1", 1);
    reorder("r1", "r2", "2143");
    update("e2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D2b-1
    reorder("e2c", "e2r_", "3412");
    mult("e2r_", "hh", "r1", 1);
    reorder("r1", "r2", "3412");
    update("e2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D2b-2
    reorder("h2c", "h2r_", "3412");
    mult("h2r_", "vh", "r1", 1);
    reorder("r1", "r2", "3412");
    reorder("r2", "r3", "2143");
    update("e2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D2c
    mult("pppgr", "s2c", "r1", 2);
    reorder("r1", "r2", "3412");
    update("e2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // D2d
    reorder("h2c", "r1", "3412");
    mult("hvhh", "r1", "r2", 2);
    reorder("r2", "r3", "2143");
    update("e2nw", 0.5, "r3");
    restore_stack_pos(pos);

    // D2e-1
    reorder("e2c", "r1", "1423");
    reorder("phph", "r2", "2341");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1342");
    update("e2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D2e-2
    reorder("t2c", "r1", "1423");
    reorder("pvgh", "r2", "2341");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1342");
    reorder("r4", "r5", "2143");
    update("e2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // D2e-3
    reorder("s2c", "r1", "1324");
    reorder("phhg", "r0", "2431");
    mult("r1", "r0", "r3", 2);
    reorder("r3", "r4", "1324");
    update("e2nw", 1.0, "r4");
    restore_stack_pos(pos);

    // D2e-4
    reorder("h2c", "r1", "1324");
    reorder("pvhp", "r0", "2431");
    mult("r1", "r0", "r3", 2);
    reorder("r3", "r4", "1324");
    reorder("r4", "r5", "2143");
    update("e2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // D3a
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "h2r_", "3412");
    reorder("pphh", "v1", "3412");
    mult("s2c", "v1", "r1", 2);
    mult("r1", "h2r_", "r2", 2);
    update("e2nw", 0.25, "r2");
    restore_stack_pos(pos);

    // D3b-1
    reorder("s2c", "s2_r1", "1324");
    reorder("h2c", "h2_r1", "1324");
    reorder("pphh", "r2", "4231");
    mult("s2_r1", "r2", "r3", 2);
    mult("r3", "h2_r1", "r4", 2);
    reorder("r4", "r5", "1324");
    update("e2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // D3b-2
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

    // D3c-1
    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 3);
    reorder("e2c", "e2_21", "2143");
    reorder("e2_21", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    reorder("r3", "r4", "2143");
    update("e2nw", -0.5, "r4");
    restore_stack_pos(pos);

    // D3c-2
    reorder("pphh", "v1", "3412");
    mult("s2c", "v1", "r1", 3);
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    update("e2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // D3d-1
    reorder("h2c", "h2r_", "3412");
    mult("h2r_", "pphh", "r1", 3);
    mult("s2c", "r1", "r2", 1);
    update("e2nw", -0.5, "r2");
    restore_stack_pos(pos);

    // D3d-2
    reorder("t2c", "t2r_", "3412");
    mult("t2r_", "pphh", "r1", 3);
    reorder("e2c", "e2_21", "2143");
    mult("e2_21", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("e2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // D4a-1
    reorder("phpg", "r1", "2341");
    mult("s1c", "r1", "r2", 1);
    update("e2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D4a-2
    reorder("pvgp", "r1", "2341");
    mult("t1c", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("e2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // D4b-1
    reorder("hvhp", "r1", "1243");
    reorder("h1c", "h1_21", "21");
    mult("r1", "h1_21", "r2", 1);
    reorder("r2", "r3", "2134");
    update("e2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D4b-2
    reorder("vhhg", "r1", "1243");
    reorder("t1c", "t1_21", "21");
    mult("r1", "t1_21", "r2", 1);
    reorder("r2", "r3", "1243");
    update("e2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D5a-1
    reorder("ph", "r1", "21");
    reorder("e2c", "e2_21", "3412");
    mult("t1c", "r1", "r2", 1);
    mult("e2_21", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("e2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D5a-2
    reorder("ph", "r1", "21");
    reorder("h2c", "h2_21", "3412");
    mult("s1c", "r1", "r2", 1);
    mult("h2_21", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    reorder("r4", "r5", "2143");
    update("e2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // D5b-1
    reorder("t1c", "t1_21", "21");
    mult("t1_21", "ph", "r1", 1);
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("e2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D5b-2
    reorder("h1c", "h1_21", "21");
    mult("h1_21", "ph", "r1", 1);
    mult("s2c", "r1", "r2", 1);
    update("e2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D5c-1
    reorder("pphg", "r1", "4312");
    mult("s1c", "r1", "r2", 1);
    reorder("t2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    reorder("r5", "r6", "2134");
    update("e2nw", -1.0, "r6");   // minus
    restore_stack_pos(pos);

    // D5c-2
    reorder("pphp", "r1", "4312");
    mult("s1c", "r1", "r2", 1);
    reorder("h2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    reorder("r5", "r6", "2143");
    update("e2nw", 1.0, "r6");
    restore_stack_pos(pos);

    // D5c-3
    reorder("pphg", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    reorder("s2c", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    update("e2nw", 1.0, "r5");
    restore_stack_pos(pos);

    // D5c-4
    reorder("e2c", "e2c_m", "1243");
    reorder("pphp", "r1", "4312");
    mult("t1c", "r1", "r2", 1);
    reorder("e2c_m", "r3", "1324");
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    reorder("r5", "r6", "1243");
    update("e2nw", 1.0, "r6");   // minus
    restore_stack_pos(pos);

    // D5d-1
    reorder("pvhh", "r1", "2314");
    reorder("h1c", "h1_21", "21");
    mult("h1_21", "r1", "i1", 1);
    reorder("t2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    reorder("r3", "r4", "2134");
    update("e2nw", 1.0, "r4");   // minus
    restore_stack_pos(pos);

    // D5d-2
    reorder("pvhh", "r1", "2314");
    reorder("t1c", "t1_21", "21");
    mult("t1_21", "r1", "i1", 1);
    reorder("h2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    reorder("r3", "r4", "2143");
    update("e2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D5d-3
    reorder("phhh", "r1", "2314");
    reorder("h1c", "h1_21", "21");
    mult("h1_21", "r1", "i1", 1);
    reorder("s2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    update("e2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D5d-4
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

    // D5e-1
    reorder("pphp", "r1", "4312");
    reorder("h1c", "h1_21", "21");
    mult("s2c", "r1", "i1", 2);
    mult("i1", "h1_21", "r2", 1);
    update("e2nw", 0.5, "r2");   // minus
    restore_stack_pos(pos);

    // D5e-2
    reorder("pphg", "r1", "4312");
    reorder("t1c", "t1_21", "21");
    mult("s2c", "r1", "i1", 2);
    mult("i1", "t1_21", "r2", 1);
    reorder("r2", "r3", "1243");
    update("e2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // D5f-1
    reorder("phhh", "r1", "2341");
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "h2r", "3412");
    mult("s1c", "r1", "i1", 1);
    mult("i1", "h2r", "r2", 2);
    update("e2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // D5f-2
    reorder("pvhh", "r1", "2341");
    reorder("h2c", "h2r", "3412");
    mult("t1c", "r1", "i1", 1);
    mult("i1", "h2r", "r2", 2);
    reorder("r2", "r3", "2143");
    update("e2nw", 0.5, "r3");
    restore_stack_pos(pos);

    // D5g-1
    reorder("pphp", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("e2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // D5g-2
    reorder("pphg", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("s2c", "i1", "r2", 1);
    update("e2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D5h-1
    reorder("phhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("e2c", "e2c_21", "2143");
    reorder("e2c_21", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    reorder("r3", "r4", "2143");
    update("e2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D5h-2
    reorder("pvhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    update("e2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D6a
    mult("pppgr", "t1c", "i1", 1);
    reorder("i1", "r1", "4123");
    mult("s1c", "r1", "r2", 1);
    update("e2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D6b
    reorder("t1c", "t1_21", "21");
    reorder("h1c", "h1_21", "21");
    mult("h1_21", "vhhh", "i1", 1);
    mult("t1_21", "i1", "i2", 1);
    reorder("i2", "r1", "3412");
    update("e2nw", 1.0, "r1");
    restore_stack_pos(pos);

    // D6c-1
    reorder("phph", "r1", "2341");
    reorder("h1c", "r2", "21");
    mult("s1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    update("e2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D6c-2
    reorder("pvgh", "r1", "2341");
    reorder("t1c", "r2", "21");
    mult("t1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    reorder("r4", "r5", "2143");
    update("e2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // D6c-3
    reorder("phhg", "r0", "2431");
    mult("s1c", "r0", "i1", 1);
    reorder("t1c", "t1_21", "21");
    mult("i1", "t1_21", "i2", 1);
    reorder("i2", "r2", "1243");
    update("e2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // D6c-4
    reorder("pvhp", "r0", "2431");
    mult("t1c", "r0", "i1", 1);
    reorder("h1c", "h1_21", "21");
    mult("i1", "h1_21", "i2", 1);
    reorder("i2", "r2", "1243");
    reorder("r2", "r3", "2143");
    update("e2nw", -1.0, "r3");
    restore_stack_pos(pos);

    // D7a
    reorder("pphh", "vr1", "3412");
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "h2_21_r", "3412");
    mult("t1c", "vr1", "r1", 1);
    mult("s1c", "r1", "r2", 1);
    mult("r2", "h2_21_r", "r3", 2);
    update("e2nw", 0.5, "r3");
    restore_stack_pos(pos);

    // D7b
    reorder("pphh", "vr1", "3412");
    reorder("t1c", "t1_r", "21");
    reorder("h1c", "h1_r", "21");
    mult("s2c", "vr1", "r1", 2);
    mult("h1_r", "r1", "r2", 1);
    mult("t1_r", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("e2nw", 0.5, "r4");
    restore_stack_pos(pos);

    // D7c-1
    reorder("pphh", "r1", "3142");
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "r2", "2413");
    reorder("t1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("s1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    reorder("r6", "r7", "1243");
    update("e2nw", -1.0, "r7");
    restore_stack_pos(pos);

    // D7c-2
    reorder("pphh", "r1", "3142");
    reorder("s2c", "s2_21", "2143");
    reorder("s2_21", "r2", "2413");
    reorder("h1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("t1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    reorder("r6", "r7", "1243");
    reorder("r7", "r8", "2143");
    update("e2nw", -1.0, "r8");
    restore_stack_pos(pos);

    // D7c-3
    reorder("pphh", "r1", "3142");
    reorder("t2c", "r2", "2413");
    reorder("h1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("s1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    update("e2nw", 1.0, "r6");   // minus
    restore_stack_pos(pos);

    // D7c-4
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

    // D7d-1
    reorder("pphh", "v_", "2431");
    mult("v_", "t1c", "r1", 2);
    reorder("r1", "r2", "21");
    mult("s1c", "r2", "r3", 1);
    reorder("h2c", "h2_21", "2143");
    reorder("h2_21", "h2_21_r", "2341");
    mult("r3", "h2_21_r", "r4", 1);
    update("e2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D7d-2
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

    // D7e-1
    reorder("pphh", "v_", "2431");
    reorder("h1c", "h1_21", "21");
    mult("v_", "t1c", "r1", 2);
    mult("h1_21", "r1", "r2", 1);
    reorder("s2c", "s2_21", "2143");
    reorder("s2_21", "s2_21_r", "1243");
    mult("s2_21_r", "r2", "r3", 1);
    reorder("r3", "r4", "2134");
    update("e2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D7e-2
    reorder("pphh", "v_", "2431");
    reorder("t1c", "t1_21", "21");
    mult("v_", "t1c", "r1", 2);
    mult("t1_21", "r1", "r2", 1);
    reorder("e2c", "e2r_", "1243");
    mult("e2r_", "r2", "r3", 1);
    reorder("r3", "r4", "1243");
    update("e2nw", -1.0, "r4");
    restore_stack_pos(pos);

    // D8a-1
    reorder("pphg", "r1", "4312");
    reorder("t1c", "t1_21", "21");
    mult("t1c", "r1", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    mult("r3", "t1_21", "r4", 1);
    reorder("r4", "r5", "1243");
    update("e2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // D8a-2
    reorder("pphp", "r1", "4312");
    reorder("h1c", "h1_21", "21");
    mult("s1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    mult("r3", "h1_21", "r4", 1);
    reorder("r4", "r5", "2134");
    update("e2nw", -1.0, "r5");
    restore_stack_pos(pos);

    // D8b-1
    reorder("phhh", "r1", "2341");
    reorder("t1c", "t1_21", "21");
    reorder("h1c", "h1_21", "21");
    mult("s1c", "r1", "i1", 1);
    mult("h1_21", "i1", "i2", 1);
    mult("t1_21", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    update("e2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // D8b-2
    reorder("pvhh", "r1", "2341");
    reorder("t1c", "t1_21", "21");
    reorder("h1c", "h1_21", "21");
    mult("t1c", "r1", "i1", 1);
    mult("t1_21", "i1", "i2", 1);
    mult("h1_21", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    reorder("r3", "r4", "2143");
    update("e2nw", 1.0, "r4");
    restore_stack_pos(pos);

    // D9
    reorder("t1c", "t1_21", "21");
    reorder("h1c", "h1_21", "21");
    reorder("pphh", "vr1", "3412");
    mult("t1c", "vr1", "r1", 1);
    mult("s1c", "r1", "r2", 1);
    mult("h1_21", "r2", "r3", 1);
    mult("t1_21", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    update("e2nw", 1.0, "r5");
    restore_stack_pos(pos);

    timer_stop("11-S2");
}


void calc_E3()
{

}


/*******************************************************************************
 * folded_1h1p
 *
 * Evaluates folded diagrams in the (1,1) sector.
 ******************************************************************************/
void folded_1h1p_singles()
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


void folded_1h1p_doubles()
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


/*******************************************************************************
 * t3corr_0h2p
 *
 * Evaluates perturbative correction +T(3) and adds it to the effective
 * interaction operator.
 ******************************************************************************/
void t3corr_1h1p()
{

}
