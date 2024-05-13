/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2024 The EXP-T developers.
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
 * Utility functions common for different FSCC models.
 */

#include "ccutils.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "codata.h"
#include "crop.h"
#include "diis.h"
#include "engine.h"
#include "methods.h"
#include "options.h"
#include "symmetry.h"
#include "utils.h"


/**
 * Finds solution of coupled cluster equations by the Jacobi procedure
 */
int solve_amplitude_equations(
        int sector_h, int sector_p,
        char *singles, char *singles_buf,
        char *doubles, char *doubles_buf,
        char *triples, char *triples_buf,
        char *veff,
        void (*construct_singles)(),
        void (*construct_doubles)(),
        void (*construct_triples)(),
        void (*construct_folded)()
        )
{
    double diff1, diff2, diff3;
    int diffmax1_idx[CC_DIAGRAM_MAX_RANK];
    int diffmax2_idx[CC_DIAGRAM_MAX_RANK];
    int diffmax3_idx[CC_DIAGRAM_MAX_RANK];
    int max_t1_idx[2];
    int max_t2_idx[4];
    int max_t3_idx[6];
    double max_t1, max_t2, max_t3;
    int iter;

    int do_singles = (singles != NULL) ? 1 : 0;
    int do_doubles = (doubles != NULL) ? 1 : 0;
    int do_triples = (triples != NULL) ? 1 : 0;

    /*
     * setup DIIS/CROP
     */
    crop_queue_t *crop_queue = NULL;
    diis_queue_t *diis_queue = NULL;
    if (cc_opts->diis_enabled) {
        diis_queue = new_diis_queue(do_singles, do_doubles, do_triples && cc_opts->diis_triples);
    }
    else if (cc_opts->crop_enabled) {
        crop_queue = new_crop_queue(do_singles, do_doubles, do_triples && cc_opts->crop_triples);
    }

    /*
     * beautiful header before iterations
     */
    print_iterations_table_header(sector_h, sector_p, do_singles, do_doubles, do_triples);

    /*
     * main CC loop
     */
    int converged = 0;
    int diverged = 0;
    double time_start = abs_time();

    for (iter = 1; iter <= cc_opts->maxiter; iter++) {
        double it_t1, it_t2;
        it_t1 = abs_time();

        /*
         * skip any calculations in the current sector
         */
        if (cc_opts->skip_sector[sector_h][sector_p]) {
            printf(" skipping solution in sector %dh%dp\n", sector_h, sector_p);
            converged = 1;
            break;
        }

        /*
         * construct next approximation to amplitudes
         */
        if (do_singles) {
            construct_singles();
        }
        if (do_doubles) {
            construct_doubles();
        }
        if (do_triples) {
            construct_triples(PT_INF);
        }

        /*
         * get effective interaction operator
         */
        if (veff) {
            if (sector_h + sector_p == 1) {
                closed(singles_buf, veff);
            }
            else if (sector_h + sector_p == 2) {
                closed(doubles_buf, veff);
            }
            else if (sector_h + sector_p == 3 && veff) {
                closed(triples_buf, veff);
            }
        }

        /*
         * add folded terms to cluster amplitudes
         */
        if (construct_folded != NULL) {
            construct_folded();
        }

        /*
         * divide amplitudes by energy denominators.
         * Intermediate Hamiltonian is also applied here.
         */
        if (do_singles) {
            diveps(singles_buf);
        }
        if (do_doubles) {
            diveps(doubles_buf);
        }
        if (do_triples) {
            diveps(triples_buf);
        }

        /*
         * some cluster amplitudes should be set to zero "by hands"
         */
        if (do_singles) {
            apply_selections(sector_h, sector_p, singles_buf);
        }
        if (do_doubles) {
            apply_selections(sector_h, sector_p, doubles_buf);
        }
        if (do_triples) {
            apply_selections(sector_h, sector_p, triples_buf);
        }

        if (cc_opts->cc_model == CC_MODEL_CCS) {
            if (do_doubles) {
                clear(doubles_buf);
            }
            if (do_triples) {
                clear(triples_buf);
            }
        }
        else if (cc_opts->cc_model == CC_MODEL_CCD) {
            if (do_singles) {
                clear(singles_buf);
            }
            if (do_triples) {
                clear(triples_buf);
            }
        }

        /*
         * calculate difference with the previous step
         */
        if (do_singles) {
            diffmax(singles, singles_buf, &diff1, diffmax1_idx);
            findmax(singles_buf, &max_t1, max_t1_idx);
        }
        if (do_doubles) {
            diffmax(doubles, doubles_buf, &diff2, diffmax2_idx);
            findmax(doubles_buf, &max_t2, max_t2_idx);
        }
        if (do_triples) {
            diffmax(triples, triples_buf, &diff3, diffmax3_idx);
            findmax(triples_buf, &max_t3, max_t3_idx);
        }

        /*
         * print results of the current iteration
         */
        printf(" %3d", iter);
        if (do_singles) {
            printf("%18.12f", diff1);
        }
        if (do_doubles) {
            printf("%18.12f", diff2);
        }
        if (do_triples) {
            printf("%18.12f", diff3);
        }
        if (do_singles) {
            printf("%12.6f", max_t1);
        }
        if (do_doubles) {
            printf("%12.6f", max_t2);
        }
        if (do_triples) {
            printf("%12.6f", max_t3);
        }

        /*
         * check for convergence
         */
        converged = 1;
        if (do_singles && fabs(diff1) > cc_opts->conv_thresh) {
            converged = 0;
        }
        if (do_doubles && fabs(diff2) > cc_opts->conv_thresh) {
            converged = 0;
        }
        if (do_triples && fabs(diff3) > cc_opts->conv_thresh) {
            converged = 0;
        }

        /*
         * check for divergence
         */
        if (do_singles && fabs(max_t1) > cc_opts->div_thresh) {
            diverged = 1;
        }
        if (do_doubles && fabs(max_t2) > cc_opts->div_thresh) {
            diverged = 1;
        }
        if (do_triples && fabs(max_t3) > cc_opts->div_thresh) {
            diverged = 1;
        }

        if (converged || diverged) {
            goto end_of_iter;
        }

        /*
         * DIIS/CROP extrapolation step
         */
        if (cc_opts->diis_enabled) {
            diis_put(diis_queue, singles_buf, singles, doubles_buf, doubles, triples_buf, triples, iter);
            if (iter >= 2) {
                diis_truncate(diis_queue, cc_opts->diis_dim);
                diis_extrapolate(diis_queue, singles_buf, doubles_buf, triples_buf);
            }
        }
        else if (cc_opts->crop_enabled) {
            crop_put(crop_queue, singles_buf, singles, doubles_buf, doubles, triples_buf, triples, iter);
            if (iter >= 2) {
                crop_truncate(crop_queue, cc_opts->crop_dim);
                crop_extrapolate(crop_queue, singles_buf, doubles_buf, triples_buf);
            }
        }

        /*
         * damping: combine new amplitudes with the previous ones
         */
        if (do_singles) {
            damping(sector_h, sector_p, singles, singles_buf, iter);
        }
        if (do_doubles) {
            damping(sector_h, sector_p, doubles, doubles_buf, iter);
        }
        if (do_triples) {
            damping(sector_h, sector_p, triples, triples_buf, iter);
        }

        /*
         * update amplitudes
         */
        if (do_singles) {
            copy(singles_buf, singles);
        }
        if (do_doubles) {
            copy(doubles_buf, doubles);
        }
        if (do_triples) {
            copy(triples_buf, triples);
        }

        /*
         * flush amplitudes to disk if needed
         */
        if (cc_opts->do_flush_iter && (iter % cc_opts->do_flush_iter == 0)) {
            save_cluster_amplitudes(sector_h, sector_p, singles, doubles, triples, veff);
        }

        /*
         * print time and memory used for the iteration
         */
        end_of_iter:
        it_t2 = abs_time();
        double iter_time = it_t2 - it_t1;
        double curr_usage = (double) cc_get_current_memory_usage() / (1024.0 * 1024.0 * 1024.0);
        double peak_usage = (double) cc_get_peak_memory_usage() / (1024.0 * 1024.0 * 1024.0);
        printf("%9.1f%8.2f/%.2f\n", iter_time, curr_usage, peak_usage);

        if (converged || diverged) {
            break;
        }
    }

    /*
     * finalize DIIS/CROP
     */
    if (diis_queue) {
        delete_diis_queue(diis_queue);
    }
    if (crop_queue) {
        delete_crop_queue(crop_queue);
    }

    /*
     * print footer
     */
    print_iterations_table_footer(do_singles, do_doubles, do_triples, iter, converged);
    if (diverged) {
        printf(" diverged!\n");
    }
    else if (!converged) {
        printf(" no convergence\n");
    }
    else {
        printf(" converged in %d iterations\n\n", iter);
    }

    printf(" average time per iteration = %.3f sec\n\n", (abs_time() - time_start) / iter);

    if (converged) {
        return EXIT_SUCCESS;
    }
    else {
        return EXIT_FAILURE;
    }
}


/**
 * Prints header banner for the sector
 */
void print_sector_banner(int sect_h, int sect_p)
{
    printf("\n\n");
    printf("\t\t\t\t*****************\n");
    printf("\t\t\t\t** Sector %dh%dp **\n", sect_h, sect_p);
    printf("\t\t\t\t*****************\n");
    printf("\n\n");
}


/**
 * Damping of cluster amplitudes
 */
void damping(int h, int p, char *old_ampl, char *new_ampl, int iter)
{
    assert(h >= 0 && h < MAX_SECTOR_RANK);
    assert(p >= 0 && p < MAX_SECTOR_RANK);

    if (cc_opts->damping[h][p].enabled &&
        iter <= cc_opts->damping[h][p].stop) {
        double alpha = cc_opts->damping[h][p].factor;
        if (cc_opts->print_level >= CC_PRINT_DEBUG) {
            printf("\n damping: %.3f %s[iter=%d] + %.3f %s[iter=%d]\n",
                   1.0 - alpha, new_ampl, iter, alpha, old_ampl, iter - 1);
        }
        add(1.0 - alpha, new_ampl, alpha, old_ampl, new_ampl);
    }
}


/**
 * Prints norm of the operator represented by the diagram 'dg_name'.
 * Norm = sqrt ( sum |a_i|^2 for a_i in diagram 'a' )
 */
void print_ampl_norm(char *label, char *dg_name)
{
    double complex norm2z = scalar_product("C", "N", dg_name, dg_name);
    double norm2 = creal(norm2z);
    double norm = sqrt(norm2);
    printf(" %s = ", label);
    if (norm > 1e-4) {
        printf("%.8f\n", norm);
    }
    else {
        printf("%.4E\n", norm);
    }
}


/**
 * Prints max (by abs value) element of diagram 'dg_name'
 */
void print_max_ampl(char *label, char *dg_name)
{
    double max_t = 0.0;
    int max_idx[CC_DIAGRAM_MAX_RANK];

    int rk = rank(dg_name);
    findmax(dg_name, &max_t, max_idx);
    printf(" %s [ ", label);
    for (int i = 0; i < rk / 2; i++) {
        printf("%d ", max_idx[i] + 1);
    }
    printf("-> ");
    for (int i = rk / 2; i < rk; i++) {
        printf("%d ", max_idx[i] + 1);
    }
    printf("] = ");
    if (max_t > 1e-4) {
        printf("%.8f\n", max_t);
    }
    else {
        printf("%.4E\n", max_t);
    }

    // print more information about these spinors in order to simplify analysis
    for (int i = 0; i < rk; i++) {
        int occ, active, repno;
        double eps;

        get_spinor_info(max_idx[i], &repno, &occ, &active, &eps);
        printf("   [%4d] eps=%16.8f rep=%s %s %s\n",
               max_idx[i] + 1, eps, get_irrep_name(repno),
               (active == 1) ? "active" : "inactive",
               (occ == 1) ? "occ" : "virt");
    }
}


/**
 * T1 diagnostic.
 * For details, see T. J. Lee, P. R. Taylor, Int. J. Quantum Chem. 23, 199 (1989)
 */
double t1_diagnostic(char *dg_name, int nelec)
{
    double t1diagn = 0.0;

    diagram_t *dg = diagram_stack_find(dg_name);
    if (dg == NULL) {
        errquit("t1_diagnostic(): diagram '%s' not found", dg_name);
    }

    assert(dg->rank == 2);

    for (size_t isb = 0; isb < dg->n_blocks; isb++) {
        block_t *sb = dg->blocks[isb];
        block_load(sb);
        for (size_t i = 0; i < sb->size; i++) {
            double abs_val = (arith == CC_ARITH_COMPLEX) ? cabs(sb->buf[i]) : fabs(((double *) sb->buf)[i]);
            t1diagn += abs_val * abs_val;
        }
        block_unload(sb);
    }

    t1diagn = sqrt(t1diagn / nelec);
    return t1diagn;
}


/**
 * Checks if one of cluster amplitudes |t| > 1
 * If |t| > 1 returns 1, if all is OK returns 0
 */
int check_divergence(char *dg_name)
{
    double max_t;
    int max_idx[CC_DIAGRAM_MAX_RANK];

    findmax(dg_name, &max_t, max_idx);

    if (max_t > cc_opts->div_thresh) {
        return 1;
    }
    else {
        return 0;
    }
}


/**
 * Returns string with the name of the RCC model
 */
void get_cc_model_name(int sect_h, int sect_p, cc_model_t cc_model, char *cc_model_name)
{
    cc_model_name[0] = '\0';

    if (sect_h > 0 || sect_p > 0) {
        strcat(cc_model_name, "FS-");
    }

    switch (cc_model) {
        case CC_MODEL_CCS:
            strcpy(cc_model_name, "CCS");
            break;
        case CC_MODEL_CCD:
            strcpy(cc_model_name, "CCD");
            break;
        case CC_MODEL_CCSD:
            strcpy(cc_model_name, "CCSD");
            break;
        case CC_MODEL_CCSD_T3:
            strcpy(cc_model_name, "CCSD+T(3)");
            break;
        case CC_MODEL_CCSD_T4:
            strcpy(cc_model_name, "CCSD+T(4)");
            break;
        case CC_MODEL_CCSDT_1A:
            strcpy(cc_model_name, "CCSDT-1a");
            break;
        case CC_MODEL_CCSDT_1B:
            strcpy(cc_model_name, "CCSDT-1b");
            break;
        case CC_MODEL_CCSDT_1B_PRIME:
            strcpy(cc_model_name, "CCSDT-1b'");
            break;
        case CC_MODEL_CCSDT_2:
            strcpy(cc_model_name, "CCSDT-2");
            break;
        case CC_MODEL_CCSDT_3:
            strcpy(cc_model_name, "CCSDT-3");
            break;
        case CC_MODEL_CCSDT:
            strcpy(cc_model_name, "CCSDT");
            break;
        default:
            strcpy(cc_model_name, "UNKNOWN");
            break;
    }

    /*
     * special case of the CCSDT-1b' model
     */
    if (cc_model == CC_MODEL_CCSDT_1B_PRIME) {
        strcpy(cc_model_name, "CCSDT-1b");
    }

    /*
     * for CCSD(T) in the 0h0p sector
     */
    if (cc_model == CC_MODEL_CCSD_T3 || cc_model == CC_MODEL_CCSD_T4) {  // for CCSD(T)
        strcpy(cc_model_name, "CCSD");
    }
    if (cc_model == CC_MODEL_CCSDT_1B_PRIME) {
        strcpy(cc_model_name, "CCSDT-1b");
    }

    if (sect_h > 0 || sect_p > 0) {
        char sector_label[10];

        sprintf(sector_label, "(%dh%dp)", sect_h, sect_p);
        strcat(cc_model_name, sector_label);
    }
}


/**
 * Print beautiful string with ionization potential in a.u., cm-1, eV to stdout.
 */
void print_ionization_potential(char *label, double ion_pot)
{
    printf(" %s = %18.12f a.u. = %8.4f eV = %10.2f cm^-1\n",
           label, ion_pot, ion_pot * CODATA_AU_TO_EV, ion_pot * CODATA_AU_TO_CM);
}


/**
 * Summarizes information about calculated cluster amplitudes:
 * - max amplitudes (by abs value)
 * - norm of cluster operators
 * - distribution analysis
 */
void print_cluster_operator_analysis(int sect_h, int sect_p, char *oper_t1, char *oper_t2, char *oper_t3)
{
    /*
     * max amplitudes
     */
    char max_label_singles[256];
    char max_label_doubles[256];
    char max_label_triples[256];

    sprintf(max_label_singles, "Max T{%dh%dp}_1 amplitude (t{%d%d}_ia)", sect_h, sect_p, sect_h, sect_p);
    sprintf(max_label_doubles, "Max T{%dh%dp}_2 amplitude (t{%d%d}_ijab)", sect_h, sect_p, sect_h, sect_p);
    sprintf(max_label_triples, "Max T{%dh%dp}_3 amplitude (t{%d%d}_ijkabc)", sect_h, sect_p, sect_h, sect_p);

    printf("\n");
    printf(" (absolute values)\n");
    if (oper_t1 != NULL) {
        print_max_ampl(max_label_singles, oper_t1);
    }
    if (oper_t2 != NULL) {
        print_max_ampl(max_label_doubles, oper_t2);
    }
    if (oper_t3 != NULL) {
        print_max_ampl(max_label_triples, oper_t3);
    }
    printf("\n");

    /*
     * norm of cluster operators
     */
    char norm_label_singles[256];
    char norm_label_doubles[256];
    char norm_label_triples[256];

    sprintf(norm_label_singles, "Norm |T{%dh%dp}_1|", sect_h, sect_p);
    sprintf(norm_label_doubles, "Norm |T{%dh%dp}_2|", sect_h, sect_p);
    sprintf(norm_label_triples, "Norm |T{%dh%dp}_3|", sect_h, sect_p);

    if (oper_t1 != NULL) {
        print_ampl_norm(norm_label_singles, oper_t1);
    }
    if (oper_t2 != NULL) {
        print_ampl_norm(norm_label_doubles, oper_t2);
    }
    if (oper_t3 != NULL) {
        print_ampl_norm(norm_label_triples, oper_t3);
    }

    /*
     * print distribution of amplitudes
     */
    if (oper_t3 != NULL) {
        print_amplitude_distribution_analysis(oper_t3);
    }

    printf("\n");
}


/**
 * for the given sector: flushes all cluster amplitudes and effective interaction to disk
 */
void save_cluster_amplitudes(int sect_h, int sect_p, char *oper_t1, char *oper_t2, char *oper_t3, char *oper_veff)
{
    char file_name[256];

    if (oper_t1) {
        sprintf(file_name, "%s.dg", oper_t1);
        diagram_write_binary(diagram_stack_find(oper_t1), file_name);
    }
    if (oper_t2) {
        sprintf(file_name, "%s.dg", oper_t2);
        diagram_write_binary(diagram_stack_find(oper_t2), file_name);
    }
    if (oper_t3) {
        sprintf(file_name, "%s.dg", oper_t3);
        diagram_write_binary(diagram_stack_find(oper_t3), file_name);
    }
    if (oper_veff) {
        sprintf(file_name, "%s.dg", oper_veff);
        diagram_write_binary(diagram_stack_find(oper_veff), file_name);
    }
}


/*
 * for the given sector: flushes all cluster amplitudes to formatted txt files
 */
void save_cluster_amplitudes_formatted(int sect_h, int sect_p, char *oper_t1, char *oper_t2, char *oper_t3)
{
    char txt_file_name_singles[256];
    char txt_file_name_doubles[256];
    char txt_file_name_triples[256];

    sprintf(txt_file_name_singles, "T1_%dh%dp.txt", sect_h, sect_p);
    sprintf(txt_file_name_doubles, "T2_%dh%dp.txt", sect_h, sect_p);
    sprintf(txt_file_name_triples, "T3_%dh%dp.txt", sect_h, sect_p);

    if (oper_t1) {
        diagram_write_formatted(diagram_stack_find(oper_t1), txt_file_name_singles);
    }
    if (oper_t2) {
        diagram_write_formatted(diagram_stack_find(oper_t2), txt_file_name_doubles);
    }
    if (oper_t3) {
        diagram_write_formatted(diagram_stack_find(oper_t3), txt_file_name_triples);
    }
}


void intruder_state_analysis(char *oper_t1, char *oper_t2, char *oper_t3)
{
    if (oper_t1) {
        predict_intruders_for_diagram(oper_t1, 5);
    }
    if (oper_t2) {
        predict_intruders_for_diagram(oper_t2, 5);
    }
    if (oper_t3) {
        predict_intruders_for_diagram(oper_t3, 5);
    }
}


void print_iterations_table_header(int sect_h, int sect_p, int singles_column, int doubles_column, int triples_column)
{
    const int COL_WIDTH_ITER = 3;
    const int COL_WIDTH_DIFFMAX = 18;
    const int COL_WIDTH_MAX = 12;
    const int COL_WIDTH_TIME = 9;
    const int COL_WIDTH_MEM = 15;

    int num_diffmax_col = singles_column + doubles_column + triples_column;
    int num_max_col = num_diffmax_col;
    int num_hyphens = COL_WIDTH_ITER +
                      num_diffmax_col * COL_WIDTH_DIFFMAX +
                      num_max_col * COL_WIDTH_MAX +
                      COL_WIDTH_TIME +
                      COL_WIDTH_MEM;

    printf(" solution of amplitude equations (sector %dh%dp)\t\t", sect_h, sect_p);
    print_asctime();

    printf(" ");
    for (int i = 0; i < num_hyphens; i++) {
        printf("-");
    }
    printf("\n");

    printf(" it.");
    printf("%s%s%s",
           singles_column ? "       diffmax(T1)" : "",
           doubles_column ? "       diffmax(T2)" : "",
           triples_column ? "       diffmax(T3)" : "");
    printf("%s%s%s",
           singles_column ? "     max(T1)" : "",
           doubles_column ? "     max(T2)" : "",
           triples_column ? "     max(T3)" : "");
    printf("    t,sec       mem,Gb\n");

    printf(" ");
    for (int i = 0; i < num_hyphens; i++) {
        printf("-");
    }
    printf("\n");
}


void print_iterations_table_footer(int singles_column, int doubles_column, int triples_column, int num_iterations,
                                   int converged)
{
    const int COL_WIDTH_ITER = 3;
    const int COL_WIDTH_DIFFMAX = 18;
    const int COL_WIDTH_MAX = 12;
    const int COL_WIDTH_TIME = 9;
    const int COL_WIDTH_MEM = 15;

    int num_diffmax_col = singles_column + doubles_column + triples_column;
    int num_max_col = num_diffmax_col;
    int num_hyphens = COL_WIDTH_ITER +
                      num_diffmax_col * COL_WIDTH_DIFFMAX +
                      num_max_col * COL_WIDTH_MAX +
                      COL_WIDTH_TIME +
                      COL_WIDTH_MEM;

    printf(" ");
    for (int i = 0; i < num_hyphens; i++) {
        printf("-");
    }
    printf("\n");
}
