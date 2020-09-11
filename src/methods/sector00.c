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
 * sector00.c
 *
 * Single-reference coupled cluster method (Fock space sector 0h0p)
 *
 * 2018-2020 Alexander Oleynichenko
 ******************************************************************************/

#include "methods.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "ccutils.h"
#include "engine.h"
#include "datamodel.h"
#include "diis.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "utils.h"

void calc_T1();

void calc_T2();

void calc_T3();

void t3corr();

double mp2_energy();

double cc_energy();

void initial_guess();

void t3_0h0p_contrib_to_doubles();
void diag_T3_0h0p_c_ab();
void diag_T3_0h0p_a_bc();
void diag_T3_0h0p_abc();
void diag_T3_0h0p_k_ij_c_ab();
void diag_T3_0h0p_i_jk();
void diag_T3_0h0p_ijk();
void diag_T3_0h0p_k_ij();
void diag_T3_0h0p_k_ij_a_bc();
void diag_T3_0h0p_k_ij_abc();
void diag_T3_0h0p_i_jk_c_ab();
void diag_T3_0h0p_ijk_c_ab();



/*******************************************************************************
 * sector00
 *
 * Solves amplitude equations and find correlation energy for ground state
 * (Fock space sector 0h0p)
 ******************************************************************************/
int sector00(cc_options_t *opts)
{
    double diff1, diff2, diff3;
    int diffmax1_idx[CC_DIAGRAM_MAX_RANK];
    int diffmax2_idx[CC_DIAGRAM_MAX_RANK];
    int diffmax3_idx[CC_DIAGRAM_MAX_RANK];
    char total_energy_str[16];
    int converged;
    int it;
    double ecorr;
    double t1, t2;
    int triples;
    char cc_model_str[32];

    printf("\n");
    printf("\t\t\t\t*****************\n");
    printf("\t\t\t\t** Sector 0h0p **\n");
    printf("\t\t\t\t*****************\n");
    printf("\n");

    opts->curr_sector_h = 0;
    opts->curr_sector_p = 0;

    triples = (opts->cc_model < CC_MODEL_CCSDT_1A) ? 0 : 1;
    if (opts->cc_model == CC_MODEL_CCSD_T3) {  // for CCSD(T)
        strcpy(cc_model_str, "CCSD");
    }
    else {
        strcpy(cc_model_str, opts->cc_model_str);
    }

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
    if (triples || opts->cc_model == CC_MODEL_CCSD_T3) {
        request_sorting("ppph", "ppph", "0000", "1234");
    }

    perform_sorting();

    // hphp -> phph
    reorder("hphp", "phph", "2143");
    set_order("phph", "1234");

    reorder("pphh", "v", "3412");

    if (opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }

    // initialize T1 and T2 amplitudes (diagrams t1c and t2c)
    initial_guess();

    //diagram_perm_struct(diagram_stack_find("t3c"));

    printf(" Solution of amplitude equations (sector 0h0p)\t\t\t");
    print_asctime();
    if (!triples) {
        printf(" --------------------------------------------------------------------------------------------\n");
        printf(" it.            E(corr)              diffmax(T1)         diffmax(T2)    t,sec       mem,Gb\n");
        printf(" --------------------------------------------------------------------------------------------\n");
    }
    else {  // triples enabled
        printf(" ----------------------------------------------------------------------------------------------------------------\n");
        printf(" it.            E(corr)              diffmax(T1)         diffmax(T2)         diffmax(T3)    t,sec       mem,Gb\n");
        printf(" ----------------------------------------------------------------------------------------------------------------\n");
    }

    diis_queue_t *diis_queue = new_diis_queue(1, 1, opts->diis_triples);
    converged = 0;
    t1 = abs_time();

    for (it = 1; it <= opts->maxiter; it++) {
        double it_t1, it_t2;
        it_t1 = abs_time();

        calc_T1();
        calc_T2();
#ifdef VERSION_DEVEL
        if (triples) {
            calc_T3();
        }
#endif

        diveps("t1nw");
        diveps("t2nw");
        if (triples) {
            diveps("t3nw");
        }

#ifdef VERSION_DEVEL
        if (opts->do_remove_inner_core_corr) {
            remove_inner_core_correlation("t2nw");
        }

        if (opts->do_relax) {
            remove_core_correlation("t2nw");
            if (triples) {
                remove_core_correlation("t3nw");
            }
        }

        if (opts->restrict_triples) {
            restrict_triples("t3nw", opts->restrict_triples_e1, opts->restrict_triples_e2);
        }
#endif

        if (opts->cc_model == CC_MODEL_CCS) {
            clear("t2nw");
        }
        else if (opts->cc_model == CC_MODEL_CCD) {
            clear("t1nw");
        }

        // for vscc
        //clear("t1nw");

        ecorr = cc_energy();
        diffmax("t1c", "t1nw", &diff1, diffmax1_idx);
        diffmax("t2c", "t2nw", &diff2, diffmax2_idx);
        if (triples) {
            diffmax("t3c", "t3nw", &diff3, diffmax3_idx);
        }

        // CCSD
        if (!triples) {
            printf(" %3d%24.15f%20.14f%20.14f", it, ecorr, diff1, diff2);
            if (fabs(diff1) < opts->conv && fabs(diff2) < opts->conv) {
                converged = 1;
                goto end_iter;
            }
        }
            // CCSDT
        else {
            printf(" %3d%24.15f%20.14f%20.14f%20.14f", it, ecorr, diff1, diff2, diff3);
            if (fabs(diff1) < opts->conv && fabs(diff2) < opts->conv &&
                fabs(diff3) < opts->conv) {
                converged = 1;
                goto end_iter;
            }
        }

        if (check_divergence("t1nw") == 1 ||
            check_divergence("t2nw") == 1) {
            printf("\n\t\t\tDIVERGED\n");
            return EXIT_FAILURE;
        }

        if (opts->diis_enabled) {
            diis_put(diis_queue, "t1nw", "t1c", "t2nw", "t2c", "t3nw", "t3c", it);
            if (it >= 2) {
                diis_truncate(diis_queue, opts->diis_dim);
                diis_extrapolate(diis_queue, "t1nw", "t2nw", "t3nw");
            }
        }

        damping(0, 0, "t1c", "t1nw", it);
        damping(0, 0, "t2c", "t2nw", it);
        if (triples) {
            damping(0, 0, "t3c", "t3nw", it);
        }

        copy("t1nw", "t1c");
        //clear("t1c");
        copy("t2nw", "t2c");
        if (triples) {
            copy("t3nw", "t3c");
        }

end_iter:

        // print time spent for iteration
        it_t2 = abs_time();
        double iter_time = it_t2 - it_t1;
        printf("%9.1f", iter_time);

        // flush amplitudes to disk if needed
        // (each opts->do_flush_iter iterations)
        if (opts->do_flush_iter && it % opts->do_flush_iter == 0) {
            diagram_write(diagram_stack_find("t1c"), "t1c.dg");
            diagram_write(diagram_stack_find("t2c"), "t2c.dg");
            if (triples) {
                diagram_write(diagram_stack_find("t3c"), "t3c.dg");
            }
        }

        // print memory statistics if the print level is high enough
        double curr_usage = cc_get_current_memory_usage() / (1024.0 * 1024.0 * 1024.0);
        double peak_usage = cc_get_peak_memory_usage() / (1024.0 * 1024.0 * 1024.0);
        printf("%8.2f/%.2f\n", curr_usage, peak_usage);
        //printf("%.3f %.3f\n", curr_usage*1e3, peak_usage*1e3);

        if (converged) {
            break;
        }
    }
    t2 = abs_time();
    if (!triples) {
        printf(" --------------------------------------------------------------------------------------------\n");
    }
    else {
        printf(" ----------------------------------------------------------------------------------------------------------------\n");
    }
    if (converged == 0) {
        printf("\tnot converged!\n");
        return EXIT_FAILURE;
    }
    else {
        printf("\tconverged in %d iterations\n", it);
    }

    delete_diis_queue(diis_queue);

    diagram_write(diagram_stack_find("t1c"), "t1c.dg");
    diagram_write(diagram_stack_find("t2c"), "t2c.dg");
    if (triples) {
        diagram_write(diagram_stack_find("t3c"), "t3c.dg");
    }

    strcpy(total_energy_str, " Total ");
    strcat(total_energy_str, cc_model_str);
    printf("\n");
    printf("          SCF reference energy = %20.12f\n", opts->escf);
    printf(" %10s correlation energy = %20.12f\n", cc_model_str, ecorr);
    printf("  %21s energy = %20.12f\n", total_energy_str, opts->escf + ecorr);
    printf("\n");

    // write total energy to the file with formatted Heff
    // file will be truncated
    FILE *hefff = fopen("HEFF", "w");
    if (carith) {
        fprintf(hefff, "complex      # arithmetic\n");
    }
    else {
        fprintf(hefff, "real         # arithmetic\n");
    }
    fprintf(hefff, "0h0p         # sector\n");
    fprintf(hefff, "   1     1   # rep No & heff size\n");
    if (carith) {
        fprintf(hefff, "%21.12E%21.12E\n", opts->escf + ecorr, 0.0);
    }
    else {
        fprintf(hefff, "%21.12E\n", opts->escf + ecorr);
    }
    fclose(hefff);

    printf(" average time per iteration = %.3f sec\n", (t2 - t1) / it);

    // extract max T1 and T2 amplitudes
    printf(" (absolute values)\n");
    print_max_ampl("Max T1 amplitude (t_ia)", "t1c");
    print_max_ampl("Max T2 amplitude (t_ijab)", "t2c");
    if (triples) {
        print_max_ampl("Max T3 amplitude (t_ijkabc)", "t3c");
    }
    printf("\n");

    // norms of cluster operators
    print_ampl_norm("Norm |T1|", "t1c");
    print_ampl_norm("Norm |T2|", "t2c");
    if (triples) {
        print_ampl_norm("Norm |T3|", "t3c");
    }
    printf("\n");

    // print T1 diagnostic
    {
        int nelec = get_num_electrons();
        double t1d = t1_diagnostic("t1c", nelec);
        printf(" T1 diagnostic = %.8f\n", t1d);
        printf(" More on T1 diagnostic: [T.J.Lee,P.R.Taylor,IJQC 23,199(1989)]\n\n");
    }

    // [T] and (T) perturbative corrections
#ifdef VERSION_DEVEL
    if (opts->cc_model == CC_MODEL_CCSD_T3) {
        t3corr();
    }
#endif

    if (triples) {
        diagram_stack_erase("t3nw");
    }

    return EXIT_SUCCESS;
}


/*******************************************************************************
 * initial_guess
 *
 * Constructs initial guess (MP2) to T1 and T2 cluster amplitudes
 * (diagrams 't1c' and 't2c', respectively) or reads these diagrams from disk.
 ******************************************************************************/
void initial_guess()
{
    int triples = (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) ? 1 : 0;
    int calc_t1 = 1, calc_t2 = 1, calc_t3 = 1;

    printf("\n Initial guess\n");
    printf(" -------------\n");

    if (cc_opts->reuse_0h0p) {
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

#ifdef VERSION_DEVEL
    if (cc_opts->restrict_triples) {
        restrict_triples("t3c", cc_opts->restrict_triples_e1, cc_opts->restrict_triples_e2);
    }
#endif
}


/*******************************************************************************
 * cc_mp2
 *
 * Calculates MP2 energy.
 ******************************************************************************/
double mp2_energy()
{
    double complex emp2 = 0.0 + 0.0 * I;
    double complex et1 = 0.0 + 0.0 * I;
    double complex et2 = 0.25 * scalar_product("N", "N", "v", "t2c");

    dg_stack_pos_t pos = get_stack_pos();
    reorder("ph", "phr", "21");
    et1 = scalar_product("N", "N", "phr", "t1c");
    restore_stack_pos(pos);

    if (cc_opts->print_level >= CC_PRINT_HIGH) {
        printf(" Contributions to the MP2 energy:\n");
        printf("   one-particle  %15.10f   two-particle  %15.10f\n", creal(et1), creal(et2));
    }

    emp2 = et1 + et2;

    return creal(emp2);
}


/*******************************************************************************
 * cc_energy
 *
 * Calculates CCSD energy on the current iteration.
 ******************************************************************************/
double cc_energy()
{
    double complex ecc = 0.0 + 0.0 * I;
    double complex et1 = 0.0 + 0.0 * I;
    double complex et1_2 = 0.0 + 0.0 * I;
    double complex et2 = 0.25 * scalar_product("N", "N", "v", "t2nw");

    dg_stack_pos_t pos = get_stack_pos();
    reorder("ph", "phr", "21");
    et1 = scalar_product("N", "N", "phr", "t1nw");
    reorder("pphh", "r1", "3142");
    mult("r1", "t1nw", "r2", 2);
    et1_2 = 0.5 * scalar_product("N", "N", "r2", "t1nw");
    restore_stack_pos(pos);

    /*if (cc_opts->print_level >= CC_PRINT_HIGH) {
        printf(" Contributions to the CC energy:\n");
        printf("   T1    %15.10f   T1^2  %15.10f   T2    %15.10f\n", creal(et1), creal(et1_2), creal(et2));
    }*/

    ecc = et1 + et1_2 + et2;
    cc_opts->eref = cc_opts->escf + creal(ecc);  // save for higher sectors

    return creal(ecc);
}


/*******************************************************************************
 * calc_T1
 *
 * Singles equations.
 ******************************************************************************/
void calc_T1()
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
#ifdef VERSION_DEVEL
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_0h0p_contrib_to_singles();
    }
#endif

    timer_stop("00-T1");
}


/*******************************************************************************
 * calc_T2
 *
 * Doubles amplitude equations.
 ******************************************************************************/
void calc_T2()
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

    // D2d
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

    // D3a
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
    /*mult("t1c", "phhp", "i1", 1);
    reorder("t1c", "t1cr", "21");
    mult("i1", "t1cr", "i2", 1);
    reorder("i2", "r2", "1243");
    perm("r2", "(12|34)");
    add(-1.0, "r2", 1.0, "t2nw", "t2nw");
    restore_stack_pos(pos);*/

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

    // TODO: too many reorderings
    /*reorder("pphh", "v_", "3241");
    reorder("t1c", "t1cr", "21");
    mult("t1c", "v_", "r1", 1);
    reorder("t2c", "t2r", "2431");
    mult("t2r", "r1", "r2", 2);
    mult("r2", "t1cr", "r3", 1);
    reorder("r3", "r4", "3142");
    reorder("r4", "r5", "2143");
    perm("r5", "(12|34)");
    add(-1.0, "r5", 1.0, "t2nw", "t2nw");
    restore_stack_pos(pos);*/

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
#ifdef VERSION_DEVEL
    if (cc_opts->cc_model >= CC_MODEL_CCSDT_1A) {
        t3_0h0p_contrib_to_doubles();
    }
#endif

    timer_stop("00-T2");
}
