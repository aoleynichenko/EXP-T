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
 * Perturbative corrections to account for triple excitations: +T(CCSD) and (T)
 * For a detailed description of the theoretical part, see:
 *   (1) I. Shavitt, R. J. Bartlett, "Many-Body Methods in Chemistry and
 *       Physics: MBPT and Coupled-Cluster Theory", 1st edition
 *   (2) M. J. O. Deegan, P. J. Knowles, Chem. Phys. Lett. V. 227, P. 321 (1994)
 */

#include "methods.h"

#include <stdio.h>
#include <string.h>

#include "engine.h"
#include "options.h"
#include "spinors.h"

void construct_triples_0h0p(int pt_order);

void t3_contrib_to_singles_0h0p(int pt_order);

void t3_contrib_to_doubles_0h0p(int pt_order);

void sector_0h0p_perturbative_triples_connected();

void sector_0h0p_perturbative_triples_disconnected();

void sector_0h0p_perturbative_triples_energy_corrections(double *de_t, double *de_st_dt);

void construct_disconnected_rank4_rank2(char *src1_name, char *src2_name, char *dst_name);

void rank2_diagram_as_matrix(diagram_t *diag, double complex *matrix);


/**
 * Calculates the (T) and [T] perturbative corrections for the ground-state CCSD
 * (FS sector 0h0p) and prints them out.
 * NOTE: this is an experimental code; permutational symmetry of T3 amplitudes
 * is not employed and T3 amplitudes by default are stored on disk;
 * so the performance of this code is still very low
 */
void sector_0h0p_perturbative_triples_correction()
{
    timer_new_entry("00-(T)", "0h0p -- CCSD(T) perturbative triples");
    timer_start("00-(T)");

    printf("\n");
    printf(" Perturbative corrections to account for triple excitations\n");
    printf(" ----------------------------------------------------------\n\n");

    printf(" disclaimer: expressions for the CCSD(T) energy correction and gradients\n");
    printf(" were derived in assumption that the orbitals are canonical/semi-canonical.\n");
    printf(" this is not so in open-shell calculations with Kramers-restricted reference.\n");
    printf(" more on the details: J. D. Watts et al, J. Chem. Phys. 88(11), 8718 (1993).\n\n");

    /*
     * construction of "connected" and "disconnected" T3 amplitudes
     */
    sector_0h0p_perturbative_triples_connected();    // => t3c
    sector_0h0p_perturbative_triples_disconnected(); // => t3(d)

    /*
     * calculate energy corrections [T] and (T)
     */
    double de_t = 0.0;
    double de_st_dt = 0.0;
    sector_0h0p_perturbative_triples_energy_corrections(&de_t, &de_st_dt);

    printf("\n");
    printf("          SCF reference energy = %20.12f\n", cc_opts->escf);
    printf("       CCSD correlation energy = %20.12f\n", cc_opts->eref - cc_opts->escf);
    printf("             Total CCSD energy = %20.12f\n", cc_opts->eref);
    printf("  4th order triples [T] corr-n = %20.12f\n", de_t);
    printf("  5th order triples (T) corr-n = %20.12f\n", de_st_dt);
    printf("     Total CCSD+T(CCSD) energy = %20.12f\n", cc_opts->eref + de_t);
    printf("          Total CCSD(T) energy = %20.12f\n", cc_opts->eref + de_t + de_st_dt);
    printf("\n");

    timer_stop("00-(T)");
}


/**
 * construction of "connected" T3 amplitudes, t3(c)
 */
void sector_0h0p_perturbative_triples_connected()
{
    double time_start = abs_time();

    // contribution T1a
    reorder("phpp", "r1", "2341");
    mult("t2c", "r1", "r2", 1);

    // contribution T1b
    reorder("t2c", "r3", "2341");
    mult("hhph", "r3", "r4", 1);

    // combine diagrams: T1a - T1b
    update("r2", -1.0, "r4");
    diagram_stack_erase("r1");
    diagram_stack_erase("r3");
    diagram_stack_erase("r4");
    reorder("r2", "t3c", "124356");
    diagram_stack_erase("r2");

    // permutation operator P(k/ij)P(a/bc)
    perm("t3c", "(3/12|4/56)");

    // divide by energy denominators
    diveps("t3c");

    double curr_mem_usage = (double) cc_get_current_memory_usage() / (1024.0 * 1024.0 * 1024.0);
    double peak_mem_usage = (double) cc_get_peak_memory_usage() / (1024.0 * 1024.0 * 1024.0);
    printf(" %-25s%8.3f sec\tmemory usage %.1f/%.1f gb\n", "time for construction of T3(c)",
           abs_time() - time_start, curr_mem_usage, peak_mem_usage);
}


/**
 * construction of "disconnected" T3 amplitudes, t3(d)
 */
void sector_0h0p_perturbative_triples_disconnected()
{
    double time_start = abs_time();

    construct_disconnected_rank4_rank2("t2c", "hp", "t3(d)");   // => E_DT
    construct_disconnected_rank4_rank2("hhpp", "t1c", "t3(d)"); // => E_ST
    perm("t3(d)", "(3/12|6/45)");
    diveps("t3(d)");

    double curr_mem_usage = (double) cc_get_current_memory_usage() / (1024.0 * 1024.0 * 1024.0);
    double peak_mem_usage = (double) cc_get_peak_memory_usage() / (1024.0 * 1024.0 * 1024.0);
    printf(" %-25s%8.3f sec\tmemory usage %.1f/%.1f gb\n", "time for construction of T3(d)",
           abs_time() - time_start, curr_mem_usage, peak_mem_usage);
}


/**
 * calculate perturbative energy corrections [T] and (T)
 * NOTE: permutational symmetry is not employed
 */
void sector_0h0p_perturbative_triples_energy_corrections(double *de_t, double *de_st_dt)
{
    double time_start = abs_time();

    // extract one-electron energies
    int nspinors = get_num_spinors();
    double *eps = (double *) cc_malloc(nspinors * sizeof(double));
    for (size_t i = 0; i < nspinors; i++) {
        eps[i] = spinor_info[i].eps;
    }

    *de_t = 0.0;
    *de_st_dt = 0.0;

    diagram_t *dg_t3c = diagram_stack_find("t3c");
    diagram_t *dg_t3d = diagram_stack_find("t3(d)");

    for (size_t isb = 0; isb < dg_t3c->n_blocks; isb++) {

        // "connected" part of triples
        block_t *block_t3c = dg_t3c->blocks[isb];
        block_load(block_t3c);
        double complex *buf_t3c = block_t3c->buf;

        // "disconnected" part of triples
        block_t *block_t3d = dg_t3d->blocks[isb];
        block_load(block_t3d);
        double complex *buf_t3d = block_t3d->buf;

        int dim_i = block_t3c->shape[0];
        int dim_j = block_t3c->shape[1];
        int dim_k = block_t3c->shape[2];
        int dim_a = block_t3c->shape[3];
        int dim_b = block_t3c->shape[4];
        int dim_c = block_t3c->shape[5];
        int *block_indices_i = block_t3c->indices[0];
        int *block_indices_j = block_t3c->indices[1];
        int *block_indices_k = block_t3c->indices[2];
        int *block_indices_a = block_t3c->indices[3];
        int *block_indices_b = block_t3c->indices[4];
        int *block_indices_c = block_t3c->indices[5];

        // energy denominators are constructed stepwise
        // to reduce the number of memory access operations
        // denom = ei + ej + ek - ea - eb - ec
        size_t index = 0;
        for (int i0 = 0; i0 < dim_i; i0++) {
            int idx_i = block_indices_i[i0];
            double denom_i = eps[idx_i];
            for (int i1 = 0; i1 < dim_j; i1++) {
                int idx_j = block_indices_j[i1];
                double denom_ij = denom_i + eps[idx_j];
                for (int i2 = 0; i2 < dim_k; i2++) {
                    int idx_k = block_indices_k[i2];
                    double denom_ijk = denom_ij + eps[idx_k];
                    for (int i3 = 0; i3 < dim_a; i3++) {
                        int idx_a = block_indices_a[i3];
                        double denom_ijka = denom_ijk - eps[idx_a];
                        for (int i4 = 0; i4 < dim_b; i4++) {
                            int idx_b = block_indices_b[i4];
                            double denom_ijkab = denom_ijka - eps[idx_b];
                            for (int i5 = 0; i5 < dim_c; i5++) {
                                int idx_c = block_indices_c[i5];
                                double denom_ijkabc = denom_ijkab - eps[idx_c];

                                double complex t3c_ijkabc = 0.0 + 0.0 * I;
                                double complex t3d_ijkabc = 0.0 + 0.0 * I;

                                if (arith == CC_ARITH_COMPLEX) {
                                    t3c_ijkabc = buf_t3c[index];
                                    t3d_ijkabc = buf_t3d[index];
                                }
                                else {
                                    t3c_ijkabc = ((double *) buf_t3c)[index];
                                    t3d_ijkabc = ((double *) buf_t3d)[index];
                                }

                                /*
                                 * dE(T) += [t(c)* + t(d)*] x t(c) x D_ijk^abc
                                 */

                                *de_t += creal(conj(t3c_ijkabc) * t3c_ijkabc) * denom_ijkabc;
                                *de_st_dt += creal(conj(t3d_ijkabc) * t3c_ijkabc) * denom_ijkabc;

                                index++;
                            }
                        }
                    }
                }
            }
        }

        block_unload(block_t3c);
        block_unload(block_t3d);
    }

    *de_t /= 36.0;
    *de_st_dt /= 36.0;

    cc_free(eps);

    printf(" %-25s%8.3f sec\n", "time for (T) energy correction", abs_time() - time_start);
}



/**
 * Constructor of disconnected diagrams:
 * <ijk|abc> = <ij|ab> * <k|c>
 */
void construct_disconnected_rank4_rank2(char *src1_name, char *src2_name, char *dst_name)
{
    /*
     * get source diagrams
     */
    diagram_t *src_diagram_1 = diagram_stack_find(src1_name);
    diagram_t *src_diagram_2 = diagram_stack_find(src2_name);

    if (src_diagram_1 == NULL) {
        errquit("construct_disconnected_rank4_rank2(): source diagram '%s' not found", src1_name);
    }
    if (src_diagram_2 == NULL) {
        errquit("construct_disconnected_rank4_rank2(): source diagram '%s' not found", src2_name);
    }

    /*
     * construct template of the destination diagram (if it doesn't exist)
     */
    diagram_t *dst_diagram = diagram_stack_find(dst_name);

    if (dst_diagram == NULL) {
        char dst_qparts[] = "xxxxxx\0";
        dst_qparts[0] = src_diagram_1->qparts[0];
        dst_qparts[1] = src_diagram_1->qparts[1];
        dst_qparts[2] = src_diagram_2->qparts[0];
        dst_qparts[3] = src_diagram_1->qparts[2];
        dst_qparts[4] = src_diagram_1->qparts[3];
        dst_qparts[5] = src_diagram_2->qparts[1];

        char dst_valence[] = "xxxxxx\0";
        dst_valence[0] = src_diagram_1->valence[0] + '0';
        dst_valence[1] = src_diagram_1->valence[1] + '0';
        dst_valence[2] = src_diagram_2->valence[0] + '0';
        dst_valence[3] = src_diagram_1->valence[2] + '0';
        dst_valence[4] = src_diagram_1->valence[3] + '0';
        dst_valence[5] = src_diagram_2->valence[1] + '0';

        tmplt(dst_name, dst_qparts, dst_valence, "123456", NOT_PERM_UNIQUE);
        dst_diagram = diagram_stack_find(dst_name);
    }


    /*
     * matrix elements of the second (one-electron) operator are saved
     * into the square matrix to make faster the access to them
     */
    int nspinors = get_num_spinors();
    double complex *t1_matrix = z_zeros(nspinors, nspinors);
    rank2_diagram_as_matrix(src_diagram_2, t1_matrix);

    /*
     * for each block:
     * for each matrix element:
     * <ijk|abc> = <ij|ab> * <k|c>
     */
    for (size_t iblock = 0; iblock < dst_diagram->n_blocks; iblock++) {
        block_t *block = dst_diagram->blocks[iblock];
        block_load(block);
        if (block->is_unique == 0) {
            continue;
        }

        int dim_i = block->shape[0];
        int dim_j = block->shape[1];
        int dim_k = block->shape[2];
        int dim_a = block->shape[3];
        int dim_b = block->shape[4];
        int dim_c = block->shape[5];

        int *block_indices_i = block->indices[0];
        int *block_indices_j = block->indices[1];
        int *block_indices_k = block->indices[2];
        int *block_indices_a = block->indices[3];
        int *block_indices_b = block->indices[4];
        int *block_indices_c = block->indices[5];

        size_t index = 0;

        for (int i1 = 0; i1 < dim_i; i1++) {
            int idx_i = block_indices_i[i1];
            for (int i2 = 0; i2 < dim_j; i2++) {
                int idx_j = block_indices_j[i2];

                // intermediate with fixed i,j indices - to make access faster
                double complex *t2_ij_matrix = z_zeros(nspinors, nspinors);
                for (int i4 = 0; i4 < dim_a; i4++) {
                    for (int i5 = 0; i5 < dim_b; i5++) {
                        int idx_a = block_indices_a[i4];
                        int idx_b = block_indices_b[i5];
                        t2_ij_matrix[idx_a * nspinors + idx_b] = diagram_get_4(src_diagram_1, idx_i, idx_j, idx_a, idx_b);
                    }
                }

                for (int i3 = 0; i3 < dim_k; i3++) {
                    int idx_k = block_indices_k[i3];
                    for (int i4 = 0; i4 < dim_a; i4++) {
                        int idx_a = block_indices_a[i4];
                        for (int i5 = 0; i5 < dim_b; i5++) {
                            int idx_b = block_indices_b[i5];
                            for (int i6 = 0; i6 < dim_c; i6++) {
                                int idx_c = block_indices_c[i6];

                                double complex t_ijab = t2_ij_matrix[idx_a * nspinors + idx_b];
                                double complex t_kc = t1_matrix[nspinors * idx_k + idx_c];
                                double complex t_ijkabc = t_ijab * t_kc;

                                if (arith == CC_ARITH_COMPLEX) {
                                    block->buf[index] += t_ijkabc;
                                }
                                else {
                                    ((double *) block->buf)[index] += creal(t_ijkabc);
                                }

                                index++;
                            }
                        }
                    }
                }

                cc_free(t2_ij_matrix);
            }
        }

        block_store(block);
    }

    /*
     * clean up
     */
    cc_free(t1_matrix);
}


/**
 * copy matrix elements of a one-electron operator into the square matrix
 */
void rank2_diagram_as_matrix(diagram_t *diag, double complex *matrix)
{
    int nspinors = get_num_spinors();

    for (int p = 0; p < nspinors; p++) {
        for (int q = 0; q < nspinors; q++) {
            matrix[p * nspinors + q] = diagram_get_2(diag, p, q);
        }
    }
}


/**
 * Constructs triples amplitudes which are correct up to 2nd PT order
 * (energy is correct up to 4th PT order)
 */
void sector_0h0p_ccsd_t3()
{
    timer_new_entry("00-T(3)", "0h0p -- CCSD+T(3) perturbative triples correction");
    timer_start("00-T(3)");

    printf(" 0h0p: 2nd order triples\n");

    tmplt("t3nw", "hhhppp", "000000", "123456", IS_PERM_UNIQUE);
    construct_triples_0h0p(PT_2);
    rename_diagram("t3nw", "t3c");
    diveps("t3c");

    timer_stop("00-T(3)");
}


/**
 * Constructs triples amplitudes which are correct up to 3rd PT order
 * (energy is correct up to 5th PT order)
 */
void sector_0h0p_ccsd_t4()
{
    dg_stack_pos_t pos;
    timer_new_entry("00-T(4)", "0h0p -- CCSD+T(4) perturbative triples correction");
    timer_start("00-T(4)");

    copy("t1c", "t1c_ccsd");
    copy("t2c", "t2c_ccsd");

    // calculate PT2 perturbative triples
    sector_0h0p_ccsd_t3();

    // calculate PT3 corrections
    tmplt("t1nw", "hp", "00", "12", NOT_PERM_UNIQUE);
    tmplt("t2nw", "hhpp", "0000", "1234", IS_PERM_UNIQUE);
    printf(" 0h0p: 3rd order correction to singles\n");
    t3_contrib_to_singles_0h0p(PT_3);
    printf(" 0h0p: 3rd order correction to doubles\n");
    t3_contrib_to_doubles_0h0p(PT_3);
    rename_diagram("t1nw", "t1_pt3_corr");
    rename_diagram("t2nw", "t2_pt3_corr");
    diveps("t1_pt3_corr");
    diveps("t2_pt3_corr");

    // update T1 and T2 with PT3 corrections
    update("t1c", 1.0, "t1_pt3_corr");
    update("t2c", 1.0, "t2_pt3_corr");
    // at this point 't1c' and 't2c' contain amplitudes which are correct up to PT3

    // calculate PT3 perturbative triples
    printf(" 0h0p: 3rd order triples\n");
    tmplt("t3nw", "hhhppp", "000000", "123456", IS_PERM_UNIQUE);
    construct_triples_0h0p(PT_3);
    rename_diagram("t3nw", "t3c");
    diveps("t3c");
    // at this point 't3c' contain amplitudes which are correct up to PT3

    timer_stop("00-T(4)");
}
