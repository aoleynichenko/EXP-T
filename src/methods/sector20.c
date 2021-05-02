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
 * sector20.c
 * ==========
 *
 * Multireference Fock-space coupled cluster method for model spaces
 * with two holes (Fock space sector 2h0p)
 *
 * Cluster Operator: T = T1 + T2 + H1 + H2 + G2
 * H10 =   |     H02 = |    \  /
 *         |           |     \/
 *        ===          =======
 *         V           V
 *         V           V
 *         |           |
 *
 * G20=    |     |
 *         V     V
 *         |     |
 *         =======
 *         V     V
 *         V     V
 *         |     |
 * (excluding amplitudes with four valence indices)
 *
 * This code was obtained by modification of the CC(0,0) program:
 * all we need is to flip particle creation lines to turn them to valence hole
 * annihilation lines. See Kaldor, J. Comp. Chem. V. 8, P.448 (1987) for details.
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

void init_amplitudes_2h0p();

void const_terms_2h0p();

void calc_G2();

void folded_2h0p();


/*******************************************************************************
 * sector20
 *
 * Solves amplitude equations in Fock space sector 2h0p
 ******************************************************************************/
int sector20(cc_options_t *opts)
{
    double diff2;
    int diffmax2_idx[CC_DIAGRAM_MAX_RANK];
    int max_t2_idx[4];
    double max_t2;
    int converged;
    int it;
    double t1, t2;

    printf("\n");
    printf("\t\t\t\t*****************\n");
    printf("\t\t\t\t** Sector 2h0p **\n");
    printf("\t\t\t\t*****************\n");
    printf("\n");

    opts->curr_sector_h = 2;
    opts->curr_sector_p = 0;

    // sorting -- interaction diagrams with valence lines
    request_sorting("hhgg", "hhhh", "0011", "1234");
    request_sorting("ppgg", "pphh", "0011", "1234");
    request_sorting("phgg", "phhh", "0011", "1234");
    perform_sorting();

    printf("\n Preparing T1, T2, S{10}_1, S{10}_2, S{20}_2 amplitudes ...\n");
    reorder("h1c", "h1r", "21");
    reorder("h2c", "h2r", "3412");

    if (opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }

    const_terms_2h0p();

    init_amplitudes_2h0p();

    printf(" Solution of amplitude equations (sector 2h0p)\n");
    print_asctime();
    printf(" ---------------------------------------------------------\n");
    printf(" it.       diffmax(S2)     max(S2)    t,sec       mem,Gb\n");
    printf(" ---------------------------------------------------------\n");

    diis_queue_t *diis_queue = new_diis_queue(0, 1, 0);
    converged = 0;
    t1 = abs_time();
    for (it = 1; it <= opts->maxiter; it++) {
        double it_t1, it_t2;
        it_t1 = abs_time();

        if (opts->skip_sector[2][0]) {
            converged = 1;
            break;
        }

        reorder("g2c", "g2cr", "3412");

        calc_G2();
        closed("g2nw", "veff20");
        folded_2h0p();
        diveps("g2nw");

        if (cc_opts->cc_model == CC_MODEL_CCS) {
            clear("g2nw");
        }

        diffmax("g2c", "g2nw", &diff2, diffmax2_idx);
        findmax("g2nw", &max_t2, max_t2_idx);
        printf(" %3d%18.12f%12.6f", it, diff2, max_t2);
        if (fabs(diff2) < opts->conv) {
            converged = 1;
            goto end_iter;
        }

        if (check_divergence("g2nw") == 1) {
            printf("\n\t\t\tDIVERGED\n");
            return EXIT_FAILURE;
        }

        if (opts->diis_enabled) {
            diis_put(diis_queue, "", "", "g2nw", "g2c", "", "", it);
            if (it >= 2) {
                diis_truncate(diis_queue, opts->diis_dim);
                diis_extrapolate(diis_queue, "", "g2nw", "");
            }
        }

        // damping
        damping(2, 0, "g2c", "g2nw", it);

        copy("g2nw", "g2c");

        end_iter:

        // print time spent for iteration
        it_t2 = abs_time();
        double iter_time = it_t2 - it_t1;
        printf("%9.1f", iter_time);

        // flush amplitudes to disk if needed
        // (each opts->do_flush_iter iterations)
        if (opts->do_flush_iter && it % opts->do_flush_iter == 0) {
            diagram_write(diagram_stack_find("g2c"), "g2c.dg");
            diagram_write(diagram_stack_find("veff20"), "veff20.dg");
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
    printf(" ---------------------------------------------------------\n");
    if (converged == 0) {
        printf("\tnot converged!\n");
        return EXIT_FAILURE;
    }
    else {
        printf("\tconverged in %d iterations\n", it);
    }

    printf(" average time per iteration = %.3f sec\n", (t2 - t1) / it);

    // extract max G2 amplitudes
    printf(" (absolute values)\n");
    print_max_ampl("Max S{20}_2 amplitude (s{20}_ijab)", "g2c");
    printf("\n");

    // norms of cluster operators
    print_ampl_norm("Norm |S{20}_2|", "g2c");
    printf("\n");

    delete_diis_queue(diis_queue);

    // flush amplitudes and effective interaction to disk
    diagram_write(diagram_stack_find("g2c"), "g2c.dg");
    diagram_write(diagram_stack_find("veff20"), "veff20.dg");

    // construct and diagonalize effective Hamiltonian
    // analyze its eigenvectors & eigenvalues
    heff_analysis(2, 0, "veff10", "veff20");
    model_space_properties_and_natural_orbitals(2, 0);

    return EXIT_SUCCESS;
}


/*******************************************************************************
 * init_amplitudes_2h0p
 *
 * Constructs initial guess to G2 cluster amplitudes (2h0p sector) or
 * read them from disk (if required)
 ******************************************************************************/
void init_amplitudes_2h0p()
{
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;
    int calc_t2 = 1, calc_t3 = 1, calc_veff = 1;

    printf("\n Initialization of S{20}_2 amplitudes ...\n");

    if (cc_opts->reuse_2h0p) {
        printf(" Trying to read amplitudes (sector 2h0p) from disk ...\n");
        if (diagram_read("g2c.dg") != NULL) {
            printf(" S{20}_2 amplitudes successfully read from disk\n");
            calc_t2 = 0;
        }
        else {
            printf(" S{20}_2 amplitudes will be calculated\n");
        }
        /*if (triples) {
            if (diagram_read("g3c.dg") != NULL) {
                printf(" S{20}_3 amplitudes successfully read from disk\n");
                calc_t3 = 0;
            }
            else {
                printf(" S{20}_3 amplitudes will be calculated\n");
            }
        }*/
    }

    // init amplitudes if needed
    if (calc_t2 && calc_veff) {
        copy("hhgg", "g2c");
        closed("g2c", "veff20");
        diveps("g2c");
        printf(" Calculated: S{20}_2 and Heff{20}\n");
    }
    else if (calc_t2) {
        copy("hhgg", "g2c");
        closed("g2c", "veff20_reserve");
        diveps("g2c");
        printf(" Calculated: S{20}_2\n");
    }
    else if (calc_veff) {
        copy("hhgg", "veff20_open");
        closed("veff20_open", "veff20");
        printf(" Calculated: Heff{20}\n");
    }
    /*if (triples && calc_t3) {
        tmplt("g3c", "???? pphppp", "???? 110000", "123456");
        printf(" Calculated: S{20}_3\n");
    }*/

    if (cc_opts->cc_model == CC_MODEL_CCS) {
        clear("g2c");
    }

    printf(" done\n\n");
}


/*******************************************************************************
 * const_terms_2h0p
 *
 ******************************************************************************/
void const_terms_2h0p()
{
}


/*******************************************************************************
 * calc_G2
 *
 * Constructs the next approximation to S{20}_2 amplitudes ("Doubles").
 ******************************************************************************/
void calc_G2()
{
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
}


/*******************************************************************************
 * folded_2h0p
 *
 * Evaluates folded diagrams in the (2,0) sector.
 ******************************************************************************/
void folded_2h0p()
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
}
