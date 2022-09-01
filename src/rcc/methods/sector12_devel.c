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
 * Fock space sector 1h2p.
 *
 * 2021-2022 Alexander Oleynichenko
 */

#include "methods.h"

#include <math.h>
#include <stdlib.h>

#include "ccutils.h"
#include "diis.h"
#include "engine.h"
#include "heff.h"
#include "datamodel.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "utils.h"

void heff12_ccsd();
void solve_amplitude_equations_1h2p();
void linsys_solver();
int iterative_solver();
void add_eps_to_diagonal(char *name);
void calc_1h2p_intermediates();

void ampl_equations_1h2p_D1();
void ampl_equations_1h2p_D2a_1();
void ampl_equations_1h2p_D2a_2();
void ampl_equations_1h2p_D2b();
void ampl_equations_1h2p_D2c();
void ampl_equations_1h2p_D2d();
void ampl_equations_1h2p_D2e_1();
void ampl_equations_1h2p_D2e_2();
void ampl_equations_1h2p_D3a();
void ampl_equations_1h2p_D3b();
void ampl_equations_1h2p_D3c();
void ampl_equations_1h2p_D3d_1();
void ampl_equations_1h2p_D3d_2();
void ampl_equations_1h2p_D4a();
void ampl_equations_1h2p_D4b_1();
void ampl_equations_1h2p_D4b_2();
void ampl_equations_1h2p_D5a();
void ampl_equations_1h2p_D5b_1();
void ampl_equations_1h2p_D5b_2();
void ampl_equations_1h2p_D5c_1();
void ampl_equations_1h2p_D5c_2();
void ampl_equations_1h2p_D5d_1();
void ampl_equations_1h2p_D5d_2();
void ampl_equations_1h2p_D5e_1();
void ampl_equations_1h2p_D5e_2();
void ampl_equations_1h2p_D5f();
void ampl_equations_1h2p_D5g_1();
void ampl_equations_1h2p_D5g_2();
void ampl_equations_1h2p_D5h();
void ampl_equations_1h2p_D6a();
void ampl_equations_1h2p_D6b();
void ampl_equations_1h2p_D6c_1();
void ampl_equations_1h2p_D6c_2();
void ampl_equations_1h2p_D7a();
void ampl_equations_1h2p_D7b();
void ampl_equations_1h2p_D7c_1();
void ampl_equations_1h2p_D7c_2();
void ampl_equations_1h2p_D7d();
void ampl_equations_1h2p_D7e_1();
void ampl_equations_1h2p_D7e_2();
void ampl_equations_1h2p_D8a_1();
void ampl_equations_1h2p_D8a_2();
void ampl_equations_1h2p_D8b();
void ampl_equations_1h2p_D9();

void ampl_equations_1h2p_folded_F1();
void ampl_equations_1h2p_folded_F2();
void ampl_equations_1h2p_folded_F3();
void ampl_equations_1h2p_folded_F4();
void ampl_equations_1h2p_folded_F5();
void ampl_equations_1h2p_folded_F6();
void ampl_equations_1h2p_folded_F7();
void ampl_equations_1h2p_folded_F8();
void ampl_equations_1h2p_folded_F9();

void diag_heff_1h2p_i_jk_c_ab__inv_h2_p13(int pt_order);
void diag_heff_1h2p_i_jk_c_ab__inv_h3_p12(int pt_order);
void diag_heff_1h2p_i_jk_c_ab__inv_h1_p23(int pt_order);
void diag_heff_1h2p_i_jk_c_ab__inv_h3_p23(int pt_order);

void diag_heff_1h2p_ijk_c_ab__inv_h2_p13(int pt_order);
void diag_heff_1h2p_ijk_c_ab__inv_h3_p12(int pt_order);
void diag_heff_1h2p_ijk_c_ab__inv_h1_p23(int pt_order);
void diag_heff_1h2p_ijk_c_ab__inv_h3_p13(int pt_order);
void diag_heff_1h2p_ijk_c_ab__inv_h2_p12(int pt_order);
void diag_heff_1h2p_ijk_c_ab__inv_h3_p23(int pt_order);

void diag_heff_1h2p_k_ij_a_bc__inv_h1_p23(int pt_order);
void diag_heff_1h2p_k_ij_a_bc__inv_h3_p12(int pt_order);
void diag_heff_1h2p_k_ij_a_bc__inv_h2_p13(int pt_order);
void diag_heff_1h2p_k_ij_a_bc__inv_h1_p12(int pt_order);

void diag_heff_1h2p_k_ij_abc__inv_h1_p23(int pt_order);
void diag_heff_1h2p_k_ij_abc__inv_h2_p13(int pt_order);
void diag_heff_1h2p_k_ij_abc__inv_h3_p12(int pt_order);
void diag_heff_1h2p_k_ij_abc__inv_h1_p12(int pt_order);
void diag_heff_1h2p_k_ij_abc__inv_h2_p12(int pt_order);
void diag_heff_1h2p_k_ij_abc__inv_h3_p13(int pt_order);


int sector12(cc_options_t *opts)
{
    printf("\n");
    printf("\t\t\t\t*****************\n");
    printf("\t\t\t\t** Sector 1h2p **\n");
    printf("\t\t\t\t*****************\n");
    printf("\n");

    opts->curr_sector_h = 1;
    opts->curr_sector_p = 2;

    request_sorting("vphh", "pphh", "1000", "1234");
    request_sorting("phpg", "phph", "0001", "1234");
    request_sorting("pppg", "ppph", "0001", "1234");
    request_sorting("ppph", "ppph", "0000", "1234");
    request_sorting("ppgh", "pphh", "0010", "1234");
    request_sorting("vvhg", "pphh", "1101", "1234");
    request_sorting("vvpg", "ppph", "1101", "1234");
    request_sorting("pvpg", "ppph", "0101", "1234");
    perform_sorting();

    tmplt("m3nw", "pphpph", "110001", "123456", IS_PERM_UNIQUE);

    solve_amplitude_equations_1h2p();

    heff12_ccsd();
    closed("m3nw", "veff12");

    heff_analysis(1, 2, "veff01", "veff10", "veff11", "veff02", "veff12");

    return EXIT_SUCCESS;
}


void solve_amplitude_equations_1h2p()
{
    //linsys_solver();
    iterative_solver();
}


void calc_1h2p_intermediates()
{
    tmplt("A", "pp", "00", "12", IS_PERM_UNIQUE);

    // I1
    add_eps_to_diagonal("A");
    update("A", 1.0, "pp");

    // I2
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r1", "2341");
    reorder("t2c", "r2", "4123");
    mult("r1", "r2", "r3", 3);
    update("A", -0.5, "r3");
    restore_stack_pos(pos);
}


void linsys_solver()
{
    // lists of:
    // active particles
    int nactp;
    moindex_t act_particles[CC_MAX_SPINORS];
    // active holes
    int nacth;
    moindex_t act_holes[CC_MAX_SPINORS];

    get_active_holes_particles(&nacth, &nactp, act_holes, act_particles);

    // all particles
    int particles[CC_MAX_SPINORS];
    int npart = 0;
    for (int i = 0; i < get_num_spinors(); i++) {
        if (is_particle(i)) {
            particles[npart++] = i;
        }
    }

    printf("act holes:\n");
    for (int i = 0; i < nacth; i++) {
        printf("%3d", act_holes[i]);
    }
    printf("\n");
    printf("act particles:\n");
    for (int i = 0; i < nactp; i++) {
        printf("%3d", act_particles[i]);
    }
    printf("\n");
    printf("particles:\n");
    for (int i = 0; i < npart; i++) {
        printf("%3d", particles[i]);
    }
    printf("\n");

    calc_1h2p_intermediates();

    prt("A");
    exit(0);

    for (int h1 = 0; h1 < nacth; h1++) {
        for (int p1 = 0; p1 < nactp; p1++) {
            for (int p2 = 0; p2 < nacth; p2++) {

                // linear system: A*t + b = 0
                double complex *t = z_zeros(1, npart);
                double complex *b = z_zeros(1, npart);
                double complex *A = z_zeros(npart, npart);

            }
        }
    }
}


void add_eps_to_diagonal(char *name)
{
    diagram_t *dg = diagram_stack_find(name);
    if (dg == NULL) {
        errquit("add_eps_to_diagonal(): diagram '%s' not found", name);
    }

    int idx2[2];
    for (size_t i = 0; i < get_num_spinors(); i++) {
        double eps = get_eps(i);
        idx2[0] = i;
        idx2[1] = i;
        diagram_set(dg, eps+0.0*I, idx2);
    }
}


int iterative_solver()
{
    double diff2;
    int diffmax2_idx[CC_DIAGRAM_MAX_RANK];
    int max_t2_idx[4];
    double max_t2;
    int it;

    tmplt("m2nw", "ppph", "1101", "1234", IS_PERM_UNIQUE);
    copy("m2nw", "m2c");
    predict_intruders_for_diagram("m2nw", 10);

    if (cc_opts->reuse_amplitudes[1][2]) {
        printf(" Trying to read amplitudes and effective interaction (sector 1h2p) from disk ...\n");
        // Singles
        if (diagram_read("m2c.dg") != NULL) {
            printf(" T{12}_1 amplitudes successfully read from disk\n");
        }
        else {
            printf(" T{12}_1 amplitudes will be calculated\n");
        }
    }
    printf("\n");

    printf(" Solution of amplitude equations (sector 1h2p)\n");
    print_asctime();
    printf(" ---------------------------------------------------------\n");
    printf(" it.       diffmax(S2)     max(S2)    t,sec       mem,Gb\n");
    printf(" ---------------------------------------------------------\n");

    diis_queue_t *diis_queue = new_diis_queue(0, 1, 0);
    int converged = 0;
    double t1 = abs_time();
    for (it = 1; it <= cc_opts->maxiter; it++) {

        double it_t1, it_t2;
        it_t1 = abs_time();

        // direct terms
        ampl_equations_1h2p_D1();
        ampl_equations_1h2p_D2a_1();
        ampl_equations_1h2p_D2a_2();
        ampl_equations_1h2p_D2b();
        ampl_equations_1h2p_D2c();
        ampl_equations_1h2p_D2d();
        ampl_equations_1h2p_D2e_1();
        ampl_equations_1h2p_D2e_2();
        ampl_equations_1h2p_D3a();
        ampl_equations_1h2p_D3b();
        ampl_equations_1h2p_D3c();
        ampl_equations_1h2p_D3d_1();
        ampl_equations_1h2p_D3d_2();
        ampl_equations_1h2p_D4a();
        ampl_equations_1h2p_D4b_1();
        ampl_equations_1h2p_D4b_2();
        ampl_equations_1h2p_D5a();
        ampl_equations_1h2p_D5b_1();
        ampl_equations_1h2p_D5b_2();
        ampl_equations_1h2p_D5c_1();
        ampl_equations_1h2p_D5c_2();
        ampl_equations_1h2p_D5d_1();
        ampl_equations_1h2p_D5d_2();
        ampl_equations_1h2p_D5e_1();
        ampl_equations_1h2p_D5e_2();
        ampl_equations_1h2p_D5f();
        ampl_equations_1h2p_D5g_1();
        ampl_equations_1h2p_D5g_2();
        ampl_equations_1h2p_D5h();
        ampl_equations_1h2p_D6a();
        ampl_equations_1h2p_D6b();
        ampl_equations_1h2p_D6c_1();
        ampl_equations_1h2p_D6c_2();
        ampl_equations_1h2p_D7a();
        ampl_equations_1h2p_D7b();
        ampl_equations_1h2p_D7c_1();
        ampl_equations_1h2p_D7c_2();
        ampl_equations_1h2p_D7d();
        ampl_equations_1h2p_D7e_1();
        ampl_equations_1h2p_D7e_2();
        ampl_equations_1h2p_D8a_1();
        ampl_equations_1h2p_D8a_2();
        ampl_equations_1h2p_D8b();
        ampl_equations_1h2p_D9();

        // folded terms
        ampl_equations_1h2p_folded_F1();
        ampl_equations_1h2p_folded_F2();
        ampl_equations_1h2p_folded_F3();
        ampl_equations_1h2p_folded_F4();
        ampl_equations_1h2p_folded_F5();
        ampl_equations_1h2p_folded_F6();
        ampl_equations_1h2p_folded_F7();
        ampl_equations_1h2p_folded_F8();
        ampl_equations_1h2p_folded_F9();

        diveps("m2nw");

        diffmax("m2c", "m2nw", &diff2, diffmax2_idx);
        findmax("m2nw", &max_t2, max_t2_idx);

        printf(" %3d%18.12f%12.6f", it, diff2, max_t2);
        if (fabs(diff2) < cc_opts->conv_thresh) {
            converged = 1;
            goto end_iter;
        }

        if (check_divergence("m2nw") == 1) {
            printf("\n\t\t\tDIVERGED\n");
            return EXIT_FAILURE;
        }

        if (cc_opts->diis_enabled) {
            diis_put(diis_queue, "", "", "m2nw", "m2c", "", "", it);
            if (it >= 2) {
                diis_truncate(diis_queue, cc_opts->diis_dim);
                diis_extrapolate(diis_queue, "", "m2nw", "");
            }
        }

        copy("m2nw", "m2c");

        end_iter:

        // print time spent for iteration
        it_t2 = abs_time();
        double iter_time = it_t2 - it_t1;
        printf("%9.1f", iter_time);

        // print memory statistics if the print level is high enough
        double curr_usage = cc_get_current_memory_usage() / (1024.0 * 1024.0 * 1024.0);
        double peak_usage = cc_get_peak_memory_usage() / (1024.0 * 1024.0 * 1024.0);
        printf("%8.2f/%.2f\n", curr_usage, peak_usage);

        if (converged) {
            break;
        }
    }
    // end of iterations

    double t2 = abs_time();

    printf(" ---------------------------------------------------------\n");

    if (converged == 0) {
        printf("\tnot converged!\n");
        return EXIT_FAILURE;
    }
    else {
        printf("\tconverged in %d iterations\n", it);
    }

    printf(" average time per iteration = %.3f sec\n", (t2 - t1) / it);

    // extract max amplitudes
    printf(" (absolute values)\n");
    print_max_ampl("Max S{12}_2 amplitude (s{12}_ijab)", "m2c");
    printf("\n");

    // norms of cluster operators
    print_ampl_norm("Norm |S{12}_2|", "m2c");
    printf("\n");

    delete_diis_queue(diis_queue);

    // flush amplitudes and effective interaction to disk
    diagram_write(diagram_stack_find("m2c"), "m2c.dg");

    //prt("m2c");
    exit(0);
}

void ampl_equations_1h2p_D1()
{
    copy("vvpg", "m2nw");
}

void ampl_equations_1h2p_D2a_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pg", "pgr", "21");
    mult("x2c", "pgr", "r1", 1);
    update("m2nw", 1.0, "r1");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D2a_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pp", "ppr_", "21");
    reorder("m2c", "m2c_21", "2143");
    mult("m2c_21", "ppr_", "r1", 1);
    reorder("r1", "r2", "2143");
    update("m2nw", 1.0, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D2b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("e2c", "e2cr_", "3412");
    mult("e2cr_", "vh", "r1", 1);
    reorder("r1", "r2", "3412");
    perm("r2", "(12)");
    update("m2nw", -1.0, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D2c()
{
    dg_stack_pos_t pos = get_stack_pos();
    mult("ppgp", "x2c", "r1", 2);
    reorder("r1", "r2", "4321");
    update("m2nw", 0.5, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D2d()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h2c", "h2c_21", "2143");
    reorder("h2c_21", "r1", "3412");
    mult("vvhh", "r1", "r2", 2);
    update("m2nw", 0.5, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D2e_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("s2c", "r1", "1324");
    mult("r1", "pvhg", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(12)");
    update("m2nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D2e_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("e2c", "e2c_1243", "1243");
    reorder("e2c_1243", "r1", "1324");
    reorder("pvhp", "r0", "2431");
    mult("r1", "r0", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(12)");
    reorder("r4", "r5", "2143"); // 1 <-> 2
    update("m2nw", -1.0, "r5");  // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D3a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h2c", "h2c_21", "2143");
    reorder("h2c_21", "h2cr_", "3412");
    reorder("pphh", "v1", "3412");
    mult("x2c", "v1", "r1", 2);
    mult("r1", "h2cr_", "r2", 2);
    update("m2nw", 0.25, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D3b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("e2c", "r1", "1423");
    reorder("pphh", "r2", "3142");
    reorder("s2c", "r3", "1324");
    mult("r1", "r2", "r4", 2);
    mult("r3", "r4", "r5", 2);
    reorder("r5", "r6", "1324");
    perm("r6", "(12)");
    update("m2nw", -1.0, "r6"); // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D3c()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "v1", "3412");
    mult("s2c", "v1", "r1", 3);
    reorder("e2c", "e2c_2134", "2134");
    reorder("e2c_2134", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    perm("r3", "(12)");
    update("m2nw", 0.5, "r3");  // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D3d_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h2c", "t2cr_", "3412");
    mult("t2cr_", "pphh", "r1", 3);
    mult("x2c", "r1", "r2", 1);
    update("m2nw", -0.5, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D3d_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t2c", "t2cr_", "3412");
    mult("t2cr_", "pphh", "r1", 3);
    reorder("m2c", "m2c_21", "2143");
    mult("m2c_21", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("m2nw", -0.5, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D4a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvpg", "r1", "2341");
    mult("s1c", "r1", "r2", 1);
    perm("r2", "(12)");
    update("m2nw", 1.0, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D4b_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("vvhg", "r1", "1243");
    reorder("t1c", "t1cr", "21");
    mult("r1", "t1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    update("m2nw", -1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D4b_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("vvhp", "r1", "1243");
    reorder("h1c", "t1cr", "21");
    mult("r1", "t1cr", "r2", 1);
    reorder("r2", "r3", "2134");   // 1243 + interchange electrons 1 <-> 2
    update("m2nw", -1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("ph", "r1", "21");
    reorder("e2c", "r2", "1342");
    mult("s1c", "r1", "r3", 1);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "1423");
    perm("r5", "(12)");
    update("m2nw", -1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5b_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("m2c", "r1", "1243");
    reorder("ph", "r2", "21");
    reorder("t1c", "r3", "21");
    mult("r1", "r2", "r4", 1);
    mult("r4", "r3", "r5", 1);
    reorder("r5", "r6", "1243");
    update("m2nw", -1.0, "r6");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5b_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("ph", "r1", "21");
    reorder("h1c", "r2", "21");
    mult("x2c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    update("m2nw", -1.0, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5c_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("ppph", "r1", "3421");
    reorder("e2c", "e2c_2134", "2134");  // sign will be changed
    reorder("e2c_2134", "r2", "2413");
    mult("s1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12)");
    update("m2nw", -1.0, "r5"); // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5c_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("s2c", "r0", "2143");
    reorder("r0", "r3", "2413");
    reorder("ppgh", "r1", "3421");
    mult("s1c", "r1", "r2", 1);
    mult("r2", "r3", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12)");
    reorder("r5", "r6", "2143");
    update("m2nw", 1.0, "r6");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5d_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2314");
    reorder("t1c", "t1cr", "21");
    mult("t1cr", "r1", "i1", 1);
    reorder("e2c", "e2c_1243", "1243"); // sign will be changed
    reorder("e2c_1243", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    perm("r3", "(12)");
    reorder("r3", "r4", "2143");
    update("m2nw", 1.0, "r4"); // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5d_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2314");
    reorder("h1c", "t1cr", "21");
    mult("t1cr", "r1", "i1", 1);
    reorder("s2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    perm("r3", "(12)");
    update("m2nw", -1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5e_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphg", "r1", "4312");
    reorder("t1c", "t1cr", "21");
    mult("x2c", "r1", "i1", 2);
    mult("i1", "t1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    update("m2nw", -0.5, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5e_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphp", "r1", "4312");
    reorder("h1c", "t1cr", "21");
    mult("x2c", "r1", "i1", 2);
    mult("i1", "t1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    reorder("r3", "r4", "2143");
    update("m2nw", -0.5, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5f()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2341");
    reorder("h2c", "h2c_2143", "2143");
    reorder("h2c_2143", "t2cr", "3412");
    mult("s1c", "r1", "i1", 1);
    mult("i1", "t2cr", "r2", 2);
    perm("r2", "(12)");
    update("m2nw", 0.5, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5g_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphp", "r1", "4231");
    reorder("m2c", "r2", "1243");
    mult("r1", "t1c", "r3", 2);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "1243");
    update("m2nw", 1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5g_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphg", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("x2c", "i1", "r2", 1);
    update("m2nw", 1.0, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5h()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("e2c", "e2c_2134", "2134");  // the sign will be changed
    reorder("e2c_2134", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    perm("r3", "(12)");
    update("m2nw", 1.0, "r3"); // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D6a()
{
    dg_stack_pos_t pos = get_stack_pos();
    mult("ppgp", "s1c", "i1", 1);
    reorder("i1", "r1", "4123");
    mult("s1c", "r1", "r2", 1);
    reorder("r2", "r3", "2143"); // interchange electrons 1 <-> 2
    update("m2nw", 1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D6b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c", "t1c_21", "21");
    reorder("h1c", "h1c_21", "21");
    mult("h1c_21", "vvhh", "i1", 1);
    mult("t1c_21", "i1", "i2", 1);
    reorder("i2", "r1", "3412");
    update("m2nw", 1.0, "r1");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D6c_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c", "r1", "21");
    mult("s1c", "pvhg", "r2", 1);
    mult("r2", "r1", "r3", 1);
    reorder("r3", "r4", "1243");
    perm("r4", "(12)");
    update("m2nw", -1.0, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D6c_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhp", "r1", "2431");
    reorder("h1c", "r2", "21");
    mult("s1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    reorder("r4", "r5", "1243");
    perm("r5", "(12)");
    reorder("r5", "r6", "2143");  // interchange electrons 1 <-> 2
    update("m2nw", -1.0, "r6");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "vr1", "3412");
    reorder("h2c", "h2c_2143", "2143");
    reorder("h2c_2143", "t2cr", "3412");
    mult("s1c", "vr1", "r1", 1);
    mult("s1c", "r1", "r2", 1);
    mult("r2", "t2cr", "r3", 2);
    update("m2nw", 0.5, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "vr1", "3412");
    reorder("t1c", "t1cr", "21");
    reorder("h1c", "h1cr", "21");
    mult("x2c", "vr1", "r1", 2);
    mult("h1cr", "r1", "r2", 1);
    mult("t1cr", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("m2nw", 0.5, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7c_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r1", "3142");
    reorder("e2c", "e2c_2134", "2134");  // the sign will be changed
    reorder("e2c_2134", "r2", "2413");
    reorder("t1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("s1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    reorder("r6", "r7", "1243");
    perm("r7", "(12)");
    update("m2nw", 1.0, "r7");  // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7c_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r1", "3142");
    reorder("s2c", "s2c_2143", "2143");
    reorder("s2c_2143", "r2", "2413");
    reorder("h1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("s1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    reorder("r6", "r7", "1243");
    perm("r7", "(12)");
    reorder("r7", "r8", "2143");
    update("m2nw", -1.0, "r8");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7d()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "v_", "2431");
    mult("v_", "t1c", "r1", 2);
    reorder("r1", "r2", "21");
    mult("s1c", "r2", "r3", 1);
    reorder("e2c", "e2c_2134", "2134"); // sign will be changed
    reorder("e2c_2134", "t2r", "2341");
    mult("r3", "t2r", "r4", 1);
    perm("r4", "(12)");
    update("m2nw", 1.0, "r4");  // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7e_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r1", "2431");
    reorder("t1c", "r2", "21");
    reorder("m2c", "r3", "1243");
    mult("r1", "t1c", "r4", 2);
    mult("r2", "r4", "r5", 1);
    mult("r3", "r5", "r6", 1);
    reorder("r6", "r7", "1243");
    update("m2nw", -1.0, "r7");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7e_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "v_", "2431");
    reorder("h1c", "h1cr", "21");
    mult("v_", "t1c", "r1", 2);
    mult("h1cr", "r1", "r2", 1);
    reorder("x2c", "x2r", "1243");
    mult("x2r", "r2", "r3", 1);
    reorder("r3", "r4", "1243");
    reorder("r4", "r5", "2143");  // interchange electrons 1 <-> 2
    update("m2nw", -1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D8a_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphp", "r1", "4312");
    reorder("h1c", "h1cr", "21");
    mult("s1c", "r1", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    mult("r3", "h1cr", "r4", 1);
    reorder("r4", "r5", "2134");  // 1243 + interchange electrons 1<->2
    update("m2nw", -1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D8a_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphg", "r1", "4312");
    reorder("t1c", "t1cr", "21");
    mult("s1c", "r1", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    mult("r3", "t1cr", "r4", 1);
    reorder("r4", "r5", "1243");
    update("m2nw", -1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D8b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2341");
    reorder("t1c", "t1cr", "21");
    reorder("h1c", "h1cr", "21");
    mult("s1c", "r1", "i1", 1);
    mult("h1cr", "i1", "i2", 1);
    mult("t1cr", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update("m2nw", 1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D9()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c", "t1cr", "21");
    reorder("h1c", "h1cr", "21");
    reorder("pphh", "vr1", "3412");
    mult("s1c", "vr1", "r1", 1);
    mult("s1c", "r1", "r2", 1);
    mult("h1cr", "r2", "r3", 1);
    mult("t1cr", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    update("m2nw", 1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("m2c", "r1", "2341");
    mult("veff01", "r1", "r2", 1);
    perm("r2", "(12)");
    update("m2nw", -1.0, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff10", "r1", "21");
    mult("m2c", "r1", "r2", 1);
    update("m2nw", 1.0, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F3()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("m2c", "r1", "3412");
    mult("veff02", "r1", "r2", 2);
    update("m2nw", -0.5, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F4()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("e1c", "e1c_21", "21");
    mult("veff02", "e1c_21", "r1", 1);
    tmplt("r1_ext", "ppph", "1101", "1234", IS_PERM_UNIQUE);
    expand_diagram("r1", "r1_ext");
    update("m2nw", -1.0, "r1_ext");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F5()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("s1c", "s1c_21", "21");
    reorder("e1c", "e1c_21", "21");
    mult("e1c_21", "veff02", "r1", 1);
    mult("s1c_21", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    update("m2nw", -1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F6()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("m2c", "r1", "2341");
    reorder("veff11", "r2", "1423");
    mult("r2", "r1", "r3", 2);
    reorder("r3", "r4", "1342");
    perm("r4", "(12)");
    update("m2nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F7()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff11", "r1", "3412");
    mult("r1", "e1c", "r2", 1);
    reorder("r2", "r3", "3412");
    tmplt("r3_ext", "ppph", "1101", "1234", NOT_PERM_UNIQUE);
    expand_diagram("r3", "r3_ext");
    perm("r3_ext", "(12)");
    update("m2nw", 1.0, "r3_ext");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F8()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff11", "r1", "1423");
    reorder("s1c", "r2", "21");
    mult("r2", "r1", "r3", 1);
    mult("r3", "e1c", "r4", 1);
    reorder("r4", "r5", "2413");
    perm("r5", "(12)");
    update("m2nw", 1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F9()
{
    dg_stack_pos_t pos = get_stack_pos();

    restore_stack_pos(pos);
}


void heff12_ccsd()
{
    timer_new_entry("12-CCSD-Heff", "1h2p -- CCSD part of eff Hamiltonian");
    timer_start("12-CCSD-Heff");

    printf("\n Construction of the CCSD eff Hamiltonian in the 1h2p sector\n");

    diag_heff_1h2p_i_jk_c_ab__inv_h2_p13(PT_INF);
    diag_heff_1h2p_i_jk_c_ab__inv_h3_p12(PT_INF);
    diag_heff_1h2p_i_jk_c_ab__inv_h1_p23(PT_INF);
    diag_heff_1h2p_i_jk_c_ab__inv_h3_p23(PT_INF);

    diag_heff_1h2p_ijk_c_ab__inv_h2_p13(PT_INF);
    diag_heff_1h2p_ijk_c_ab__inv_h3_p12(PT_INF);
    diag_heff_1h2p_ijk_c_ab__inv_h1_p23(PT_INF);
    diag_heff_1h2p_ijk_c_ab__inv_h3_p13(PT_INF);
    diag_heff_1h2p_ijk_c_ab__inv_h2_p12(PT_INF);
    diag_heff_1h2p_ijk_c_ab__inv_h3_p23(PT_INF);

    diag_heff_1h2p_k_ij_a_bc__inv_h1_p23(PT_INF);
    diag_heff_1h2p_k_ij_a_bc__inv_h3_p12(PT_INF);
    diag_heff_1h2p_k_ij_a_bc__inv_h2_p13(PT_INF);
    //diag_heff_1h2p_k_ij_a_bc__inv_h1_p12(int pt_order);   // S2_1h2p is required

    diag_heff_1h2p_k_ij_abc__inv_h1_p23(PT_INF);
    diag_heff_1h2p_k_ij_abc__inv_h2_p13(PT_INF);
    diag_heff_1h2p_k_ij_abc__inv_h3_p12(PT_INF);
    //diag_heff_1h2p_k_ij_abc__inv_h1_p12(PT_INF);
    diag_heff_1h2p_k_ij_abc__inv_h2_p12(PT_INF);
    diag_heff_1h2p_k_ij_abc__inv_h3_p13(PT_INF);

    timer_stop("12-CCSD-Heff");
}


/*
 * Final T3 permutation: (1/23|6/45) aka (i/jk|c/ab)
 * Diagrams included:
 * (i/jk|c/ab): T1b T3a T3d T4b T5d T7c T8a T8e T10a
 */
void diag_heff_1h2p_i_jk_c_ab__inv_h2_p13(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hpph", "0100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("hvhp", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3a
        reorder("ph", "phr_", "21");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("s2c", "s2c_21", "2143");
        reorder("pphp", "r2", "4312");
        mult("s2c_21", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T4b
        reorder("t1c", "r2", "21");
        mult("hvhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T7c
        reorder("pphp", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        mult("t1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("t1c", "r3", "21");
        reorder("s2c", "s2c_21", "2143");
        mult("s2c_21", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        // T10a
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("s1c", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    perm("r4", "(13|46)");
    reorder("r4", "r5", "132465"); // interchange electrons 2 <-> 3
    update("m3nw", 1.0, "r5");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_i_jk_c_ab__inv_h3_p12(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phhh", "1010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("vhhg", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3a
        reorder("ph", "phr_", "21");
        reorder("e2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("pphg", "r2", "4312");
        mult("s2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T4b
        reorder("h1c", "r2", "21");
        mult("vhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T7c
        reorder("pphg", "r1", "4312");
        mult("t1c", "r1", "r2", 1);
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("e2c", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("h1c", "r3", "21");
        mult("s2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        // T10a
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("t1c", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    perm("r4", "(12)");
    update("m3nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_i_jk_c_ab__inv_h1_p23(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "ppph", "1100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("vvhp", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3a
        reorder("ph", "phr_", "21");
        reorder("x2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("pphp", "r2", "4312");
        mult("x2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T4b
        reorder("t1c", "r2", "21");
        mult("vvhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T7c
        reorder("pphp", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("x2c", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("t1c", "r3", "21");
        mult("x2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        // T10a
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("s1c", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    perm("r4", "(56)");
    reorder("r4", "r5", "321654");
    update("m3nw", 1.0, "r5");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_i_jk_c_ab__inv_h3_p23(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pphh", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("vvhg", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3a
        /*reorder("ph", "phr_", "21");
        reorder("t2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);*/

        // T3d
        reorder("pphg", "r2", "4312");
        mult("x2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T4b
        reorder("h1c", "r2", "21");
        mult("vvhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T7c
        reorder("pphg", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        /*reorder("pphh", "r1", "4231");
        reorder("t2c", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);*/

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("h1c", "r3", "21");
        mult("x2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        // T10a
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("s1c", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    reorder("r4", "r5", "231546");
    update("m3nw", -1.0, "r5");  // the sign is changed from +1 to -1
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123|6/45) aka (ijk|c/ab)
 * Diagrams included:
 * (ijk|c/ab):  T3c T4c T7b T8b
 */
void diag_heff_1h2p_ijk_c_ab__inv_h2_p13(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hpph", "0100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        reorder("pvhp", "r0", "2431");
        mult("t1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("t1c", "r1", "21");
        reorder("hphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(13|46)");
    reorder("r4", "r5", "132465");
    update("m3nw", 1.0, "r5");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_ijk_c_ab__inv_h3_p12(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phhh", "1010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        reorder("phhg", "r0", "2431");
        mult("s1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("h1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(12)");
    update("m3nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_ijk_c_ab__inv_h1_p23(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "ppph", "1100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        reorder("pvhp", "r0", "2431");
        mult("s1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("t1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(23|56)");
    reorder("r4", "r5", "321654");
    update("m3nw", 1.0, "r5");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_ijk_c_ab__inv_h3_p13(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphh", "0110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4"); // change sign: x -1
        restore_stack_pos(pos2);

        // T4c
        mult("t1c", "pvhg", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("h1c", "r1", "21");
        reorder("hphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("e2c", "e2c_2134", "2134");
        reorder("pphh", "r1", "3142");
        reorder("e2c_2134", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", 1.0, "r4"); // change sign: x -1
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(13)");
    reorder("r4", "r5", "132456");
    update("m3nw", -1.0, "r5"); // change sign: x -1
    restore_stack_pos(pos);
}

void diag_heff_1h2p_ijk_c_ab__inv_h2_p12(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phph", "1000", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("t2c", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        mult("s1c", "phhp", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("t1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("t2c", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(12|46)");
    reorder("r4", "r5", "123465");
    update("m3nw", -1.0, "r5"); // change sign: x -1
    restore_stack_pos(pos);
}

void diag_heff_1h2p_ijk_c_ab__inv_h3_p23(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pphh", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");  // change sign: x -1
        restore_stack_pos(pos2);

        // T4c
        mult("s1c", "pvhg", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("h1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("e2c", "e2c_2134", "2134");
        reorder("pphh", "r1", "3142");
        reorder("e2c_2134", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");  // change sign: x -1
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(23)");
    reorder("r4", "r5", "231546");
    update("m3nw", -1.0, "r5"); // change sign: x -1
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|4/56) aka (k/ij|a/bc)
 * Diagrams included:
 * (k/ij|a/bc): T1a T3e T4a T5e T7d T8d T10b
 */
void diag_heff_1h2p_k_ij_a_bc__inv_h1_p23(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pppp", "1000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("pvpp", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("t2c", "r1", "3412");
        mult("r1", "pvhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        mult("ppppr", "s1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("t1c", "r1", "21");
        mult("r1", "pvhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("t2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("s1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("t1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    perm("r3", "(23)");
    reorder("r3", "r4", "321654");
    update("m3nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_k_ij_a_bc__inv_h3_p12(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphp", "0010", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("phpg", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "3412");
        mult("r1", "phhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        reorder("pppg", "pppgr", "3412");
        mult("pppgr", "t1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("h1c", "r0", "21");
        reorder("t1c", "r1", "21");
        mult("r0", "phhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("h1c", "r0", "21");
        reorder("t1c", "r1", "21");
        mult("r0", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("x2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    perm("r3", "(45)");
    update("m3nw", 1.0, "r3");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_k_ij_a_bc__inv_h2_p13(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phpp", "1100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("pvgp", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("h2c", "r1", "3412");
        mult("r1", "pvhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        mult("ppgp", "s1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("t1c", "t1c_21", "21");
        reorder("h1c", "h1c_21", "21");
        mult("t1c_21", "pvhh", "r2", 1);
        mult("h1c_21", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("h2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("s1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("h1c", "h1c_21", "21");
        reorder("t1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("h1c_21", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("s2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    perm("r3", "(13|46)");
    reorder("r3", "r4", "132465");
    update("m3nw", 1.0, "r4");
    restore_stack_pos(pos);
}

// S_2^{1h2p} is required
/*void diag_heff_1h2p_k_ij_a_bc__inv_h1_p12(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hppp", "0000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("phpp", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("t2c", "r1", "3412");
        mult("r1", "phhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        mult("ppppr", "t1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("t1c", "r1", "21");
        mult("r1", "phhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("t2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("t1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    perm("r3", "(3/12|4/56)");
    update("t3nw", 1.0, "r3");
    restore_stack_pos(pos);
}*/


/*
 * Final T3 permutation: (3/12|456) aka (k/ij|abc)
 * Diagrams included:
 * (k/ij|abc):  T3b T4d T7a T8c
 */
void diag_heff_1h2p_k_ij_abc__inv_h1_p23(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pppp", "1000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppph", "r1", "1342");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("pvhp", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("t1c", "r1", "21");
        reorder("ppph", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c", "r1", "21");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(23|56)");
    reorder("r3", "r4", "321654");  // interchange electrons 1 <-> 3
    update("m3nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_k_ij_abc__inv_h2_p13(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phpp", "1100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("pvhp", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("t1c", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    mult("s2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(13|46)");
    reorder("r3", "r4", "132465"); // interchange electrons 2 <-> 3
    update("m3nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_k_ij_abc__inv_h3_p12(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphp", "0010", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppph", "r1", "1342");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("phhg", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("h1c", "r1", "21");
        reorder("ppph", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c", "r1", "21");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    mult("x2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(45)");
    update("m3nw", 1.0, "r3");
    restore_stack_pos(pos);
}

// S_2^{1h2p} is required
void diag_heff_1h2p_k_ij_abc__inv_h1_p12(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hppp", "0000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b
        reorder("ppph", "r1", "1342");
        reorder("t2c", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4d
        reorder("phhp", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7a
        reorder("t1c", "r1", "21");
        reorder("ppph", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c", "r1", "21");
        reorder("t2c", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|456)");
    update("t3nw", 1.0, "r3");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_k_ij_abc__inv_h2_p12(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhpp", "0100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("t2c", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("phhp", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("t1c", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
        reorder("t2c", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    mult("x2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(46)");
    reorder("r3", "r4", "123465");
    update("m3nw", -1.0, "r4"); // change sign
    restore_stack_pos(pos);
}

void diag_heff_1h2p_k_ij_abc__inv_h3_p13(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pphp", "1010", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppph", "r1", "1342");
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", -1.0, "r5");  // change sign: x -1
        restore_stack_pos(pos2);

        // T4d
        reorder("pvhg", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("h1c", "r1", "21");
        reorder("ppph", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c", "r1", "21");
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", 1.0, "r7");  // change sign: x -1
        restore_stack_pos(pos2);
    }

    finish:
    mult("s2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(13|45)");
    reorder("r3", "r4", "132456");
    update("m3nw", -1.0, "r4");  // change sign: x -1
    restore_stack_pos(pos);
}
