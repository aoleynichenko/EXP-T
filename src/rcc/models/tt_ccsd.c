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


/*
 * Interface to the tensor-train-based CCSD codes of Artem S. Rumyantsev
 */

#include "tt_ccsd.h"

#ifdef TENSOR_TRAIN

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "cc_properies.h"
#include "ccutils.h"
#include "engine.h"
#include "diis.h"
#include "heff.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"
#include "timer.h"
#include "io.h"

#include "../engine/tt.h"

void print_cc_energy(double escf, double ecorr);


double estimated_ram_gb(int nspinors)
{
    double n = nspinors;
    return n * n * n * n * 16 / (1024.0 * 1024.0 * 1024.0);
}


double complex kramers_unrestricted_mp2_energy(int nspinors, int nocc, double complex *eri, double complex *fock, double *eps)
{
    double complex e_corr_mp2 = 0.0 + 0.0 * I;
    size_t nspinors64 = nspinors;
    size_t nsp3 = nspinors64 * nspinors64 * nspinors64;
    size_t nsp2 = nspinors64 * nspinors64;
    size_t nsp1 = nspinors64;

    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int a = nocc; a < nspinors; a++) {
                for (int b = nocc; b < nspinors; b++) {
                    double complex v_ijab = eri[nsp3 * i + nsp2 * j + nsp1 * a + b];
                    double d_ijab = eps[i] + eps[j] - eps[a] - eps[b];
                    e_corr_mp2 += 0.25 * conj(v_ijab) * v_ijab / d_ijab;
                }
            }
        }
    }

    return e_corr_mp2;
}


int artem_to_dirac(int nspinors, int i)
{
    if (i % 2 == 0) { // alpha
        return i / 2;
    }
    else { // beta
        return nspinors / 2 + (i - 1) / 2;
    }
}


int dirac_to_artem(int nspinors, int i)
{
    if (i < nspinors / 2) {
        return 2 * i;
    }
    else {
        return 2 * (i - nspinors / 2) + 1;
    }
}


void extract_two_electron_integrals(char *diagram_name, double *eri_ex)
{
    int nspinors = get_num_spinors();

    assert_diagram_exists(diagram_name);
    diagram_t *dg_eri = diagram_stack_find(diagram_name);

    for (int iblock = 0; iblock < dg_eri->n_blocks; iblock++) {
        block_t *block = dg_eri->blocks[iblock];
        if (!block->is_unique) {
            continue;
        }
        block_load(block);

        int dims_0 = block->shape[0];
        int dims_1 = block->shape[1];
        int dims_2 = block->shape[2];
        int dims_3 = block->shape[3];

        int coef_0 = dims_1 * dims_2 * dims_3;
        int coef_1 = dims_2 * dims_3;
        int coef_2 = dims_3;
        //int coef_3 = 1;

        int *block_indices_0 = block->indices[0];
        int *block_indices_1 = block->indices[1];
        int *block_indices_2 = block->indices[2];
        int *block_indices_3 = block->indices[3];

        // loop over matrix elements
        for (size_t i = 0; i < block->size; i++) {
            // if matrix element is zero we can skip it
            double v_ijkl = ((double *) block->buf)[i];
            if (fabs(v_ijkl) < 1e-16) {
                continue;
            }

            // calculate compound index
            int64_t index = i;
            int idx_0 = (int) (index / (int64_t) coef_0);
            index = index % ((int64_t) coef_0);
            int idx_1 = (int) (index / (int64_t) coef_1);
            index = index % ((int64_t) coef_1);
            int idx_2 = (int) (index / (int64_t) coef_2);
            index = index % ((int64_t) coef_2);
            int idx_3 = (int) index;    // since coef_3 == 1

            // convert relative indices into spinor indices
            idx_0 = block_indices_0[idx_0];
            idx_1 = block_indices_1[idx_1];
            idx_2 = block_indices_2[idx_2];
            idx_3 = block_indices_3[idx_3];

            int art_0 = dirac_to_artem(nspinors, idx_0);
            int art_1 = dirac_to_artem(nspinors, idx_1);
            int art_2 = dirac_to_artem(nspinors, idx_2);
            int art_3 = dirac_to_artem(nspinors, idx_3);

            int shape[] = {nspinors, nspinors, nspinors, nspinors};
            int art_idx[4];

            art_idx[0] = art_0;
            art_idx[1] = art_1;
            art_idx[2] = art_2;
            art_idx[3] = art_3;
            size_t idx_ijkl = tensor_index_to_linear(4, shape, art_idx);
            eri_ex[idx_ijkl] = v_ijkl;

            art_idx[0] = art_1;
            art_idx[1] = art_0;
            art_idx[2] = art_2;
            art_idx[3] = art_3;
            idx_ijkl = tensor_index_to_linear(4, shape, art_idx);
            eri_ex[idx_ijkl] = - v_ijkl;

            art_idx[0] = art_0;
            art_idx[1] = art_1;
            art_idx[2] = art_3;
            art_idx[3] = art_2;
            idx_ijkl = tensor_index_to_linear(4, shape, art_idx);
            eri_ex[idx_ijkl] = - v_ijkl;

            art_idx[0] = art_1;
            art_idx[1] = art_0;
            art_idx[2] = art_3;
            art_idx[3] = art_2;
            idx_ijkl = tensor_index_to_linear(4, shape, art_idx);
            eri_ex[idx_ijkl] = v_ijkl;
        }

        block_unload(block);
    }
}


void extract_two_electron_integrals_complex(char *diagram_name, double complex *eri_ex, int do_conj)
{
    int nspinors = get_num_spinors();

    assert_diagram_exists(diagram_name);
    diagram_t *dg_eri = diagram_stack_find(diagram_name);

    for (int iblock = 0; iblock < dg_eri->n_blocks; iblock++) {
        block_t *block = dg_eri->blocks[iblock];
        if (!block->is_unique) {
            continue;
        }
        block_load(block);

        int dims_0 = block->shape[0];
        int dims_1 = block->shape[1];
        int dims_2 = block->shape[2];
        int dims_3 = block->shape[3];

        int coef_0 = dims_1 * dims_2 * dims_3;
        int coef_1 = dims_2 * dims_3;
        int coef_2 = dims_3;
        //int coef_3 = 1;

        int *block_indices_0 = block->indices[0];
        int *block_indices_1 = block->indices[1];
        int *block_indices_2 = block->indices[2];
        int *block_indices_3 = block->indices[3];

        // loop over matrix elements
        for (size_t i = 0; i < block->size; i++) {
            // if matrix element is zero we can skip it
            double complex v_ijkl = block->buf[i];
            if (cabs(v_ijkl) < 1e-16) {
                continue;
            }

            // calculate compound index
            int64_t index = i;
            int idx_0 = (int) (index / (int64_t) coef_0);
            index = index % ((int64_t) coef_0);
            int idx_1 = (int) (index / (int64_t) coef_1);
            index = index % ((int64_t) coef_1);
            int idx_2 = (int) (index / (int64_t) coef_2);
            index = index % ((int64_t) coef_2);
            int idx_3 = (int) index;    // since coef_3 == 1

            // convert relative indices into spinor indices
            idx_0 = block_indices_0[idx_0];
            idx_1 = block_indices_1[idx_1];
            idx_2 = block_indices_2[idx_2];
            idx_3 = block_indices_3[idx_3];

            int art_0 = dirac_to_artem(nspinors, idx_0);
            int art_1 = dirac_to_artem(nspinors, idx_1);
            int art_2 = dirac_to_artem(nspinors, idx_2);
            int art_3 = dirac_to_artem(nspinors, idx_3);

            int shape[] = {nspinors, nspinors, nspinors, nspinors};
            int art_idx[4];

            art_idx[0] = art_0;
            art_idx[1] = art_1;
            art_idx[2] = art_2;
            art_idx[3] = art_3;
            size_t idx_ijkl = tensor_index_to_linear(4, shape, art_idx);
            eri_ex[idx_ijkl] = do_conj ? conj(v_ijkl) : v_ijkl;

            art_idx[0] = art_1;
            art_idx[1] = art_0;
            art_idx[2] = art_2;
            art_idx[3] = art_3;
            idx_ijkl = tensor_index_to_linear(4, shape, art_idx);
            eri_ex[idx_ijkl] = do_conj ? conj(- v_ijkl) : (- v_ijkl);

            art_idx[0] = art_0;
            art_idx[1] = art_1;
            art_idx[2] = art_3;
            art_idx[3] = art_2;
            idx_ijkl = tensor_index_to_linear(4, shape, art_idx);
            eri_ex[idx_ijkl] = do_conj ? conj(- v_ijkl) : (- v_ijkl);

            art_idx[0] = art_1;
            art_idx[1] = art_0;
            art_idx[2] = art_3;
            art_idx[3] = art_2;
            idx_ijkl = tensor_index_to_linear(4, shape, art_idx);
            eri_ex[idx_ijkl] = do_conj ? conj(v_ijkl) : v_ijkl;
        }

        block_unload(block);
    }
}


#define CHUNK_SIZE_BYTES (4096 * 4)
#define CHUNK_SIZE_COMPLEX ((CHUNK_SIZE_BYTES) / (sizeof(double complex)))


void tt_ccsd()
{
    int nspinors = get_num_spinors();
    int nocc = get_num_electrons();
    double tol = cc_opts->tt_options.svd_thresh;
    double tol_eri = cc_opts->tt_options.cholesky_thresh;

    printf("\n");
    printf("\t\t\t**************************************\n");
    printf("\t\t\t          tensor-train ccsd           \n");
    printf("\t\t\t a. s. rumyantsev, a. v. oleynichenko \n");
    printf("\t\t\t**************************************\n");
    printf("\n");
    printf("\t\t\t tensor trains             %10s\n", cc_opts->tt_options.use_tensor_trains ? "enabled" : "disabled");
    printf("\t\t\t maxiter                   %10d\n", cc_opts->maxiter);
    printf("\t\t\t conv_thresh               %10.1e\n", cc_opts->conv_thresh);
    printf("\t\t\t ttsvd tolerance           %10.1e\n", cc_opts->tt_options.svd_thresh);
    printf("\t\t\t cholesky tolerance        %10.1e\n", cc_opts->tt_options.cholesky_thresh);
    printf("\t\t\t arithmetic                   %7s\n", (arith == CC_ARITH_REAL) ? "real" : "complex");
    printf("\t\t\t use pyscf integrals          %7s\n", cc_opts->tt_options.use_pyscf_integrals ? "yes" : "no");
    printf("\t\t\t path to pyscf integrals      %s\n", cc_opts->tt_options.pyscf_integrals_path);
    printf("\n");


    /*if (cc_opts->tt_options.use_pyscf_integrals) {

        printf("\ntesting routine for integral sorting\n\n");

        int nspinors = 0;
        int nocc = 0;
        double enuc = 0.0;
        double escf = 0.0;
        double *eps = NULL;
        int *occ = NULL;
        double complex *eri = NULL;
        size_t *eri_index = NULL;

        read_pyscf_integrals(cc_opts->tt_options.pyscf_integrals_path, &nspinors, &nocc, &enuc, &escf, &eps, &occ,
                                  &eri, &eri_index);

        double complex ecorr_mp2 = 0.0 + 0.0 * I;
        for (int i = 0; i < nocc; i++) {
            for (int j = 0; j < nocc; j++) {
                for (int a = nocc; a < nspinors; a++) {
                    for (int b = nocc; b < nspinors; b++) {
                        double complex v_ijab = get_eri(nspinors, eri, eri_index, i, j, a, b);
                        double d_ijab = eps[i] + eps[j] - eps[a] - eps[b];
                        ecorr_mp2 += 0.25 * conj(v_ijab) * v_ijab / d_ijab;
                    }
                }
            }
        }

        printf("corr energy mp2 = %25.16f%25.16f\n", creal(ecorr_mp2), cimag(ecorr_mp2));

        printf("\nend of testing routine for integral sorting\n\n");
    }*/

    /*
     * tensor of two-electron integrals
     */
    void *eri_ex = NULL;
    void *fm_ex = NULL;
    printf("                                                  sec\n");

    /*
     * read integrals from pyscf output - binary file
     */
    if (cc_opts->tt_options.use_pyscf_integrals) {
        printf(" reading transformed integrals from pyscf  ");
        double t_begin_eris = abs_time();
        arith = CC_ARITH_COMPLEX;

        char pyscf_integrals_path[256];
        strcpy(pyscf_integrals_path, "../");
        strcat(pyscf_integrals_path, cc_opts->tt_options.pyscf_integrals_path);
        int f = io_open(pyscf_integrals_path, "r");
        if (f == -1) {
            errquit("cannot open file with pyscf integrals '%s'", cc_opts->tt_options.pyscf_integrals_path);
        }

        nspinors = 0;
        io_read(f, &nspinors, sizeof(int));
        int64_t nspinors64 = nspinors;

        eri_ex = calloc(nspinors64 * nspinors64 * nspinors64 * nspinors64, sizeof(double complex));
        if (eri_ex == NULL) {
            errquit("not enough memory (estimated ram usage: %.3f gb)", estimated_ram_gb(nspinors));
        }

        double *pyscf_eps = calloc(nspinors64, sizeof(double));
        int *occ = calloc(nspinors64, sizeof(int));
        double enuc = 0.0;
        double escf = 0.0;
        nocc = 0;

        io_read(f, &enuc, sizeof(double));
        io_read(f, &escf, sizeof(double));
        cc_opts->escf = escf;

        for (int i = 0; i < nspinors64; i++) {
            io_read(f, occ + i, sizeof(int));
            io_read(f, pyscf_eps + i, sizeof(double));
            nocc += occ[i];
        }

        double complex *buf = calloc(CHUNK_SIZE_COMPLEX, sizeof(double complex));
        size_t ibuf = CHUNK_SIZE_COMPLEX;

        size_t nsp3 = nspinors64 * nspinors64 * nspinors64;
        size_t nsp2 = nspinors64 * nspinors64;
        size_t nsp1 = nspinors64;

        for (int p = 0; p < nspinors64; p++) {
            for (int q = p + 1; q < nspinors64; q++) {
                for (int r = p; r < nspinors64; r++) {
                    for (int s = r + 1; s < nspinors64; s++) {

                        if (ibuf == CHUNK_SIZE_COMPLEX) {
                            // load a new chunk of raw binary data
                            io_read_unsafe(f, buf, CHUNK_SIZE_BYTES);
                            ibuf = 0;
                        }

                        double complex val = buf[ibuf];
                        double complex conj_val = conj(val);
                        ibuf++;

                        double complex *dst = (double complex *) eri_ex;

                        dst[nsp3 * p + nsp2 * q + nsp1 * r + s] = val;
                        dst[nsp3 * q + nsp2 * p + nsp1 * r + s] = - val;
                        dst[nsp3 * p + nsp2 * q + nsp1 * s + r] = - val;
                        dst[nsp3 * q + nsp2 * p + nsp1 * s + r] = val;

                        dst[nsp3 * r + nsp2 * s + nsp1 * p + q] = conj_val;
                        dst[nsp3 * r + nsp2 * s + nsp1 * q + p] = - conj_val;
                        dst[nsp3 * s + nsp2 * r + nsp1 * p + q] = - conj_val;
                        dst[nsp3 * s + nsp2 * r + nsp1 * q + p] = conj_val;
                    }
                }
            }
        }

        io_close(f);

        /*
         * reading integrals to the 1d array, employing 8-fold permutational symmetry
         */

        printf("%10.2f\n", abs_time() - t_begin_eris);

        /*
         * Fock matrix
         * for canonical orbitals only
         */
        printf(" fock matrix construction                  ");
        double t_begin_fock = abs_time();

        fm_ex = calloc(nspinors64 * nspinors64, sizeof(double complex));
        for (int i = 0; i < nspinors64; i++) {
            ((double complex *) fm_ex)[nspinors64 * i + i] = pyscf_eps[i] + 0.0 * I;
        }

        printf("%10.2f\n", abs_time() - t_begin_fock);

        /*
         * calculate mp2 energy
         */
        double complex e_corr_mp2 = kramers_unrestricted_mp2_energy(nspinors, nocc, eri_ex, fm_ex, pyscf_eps);

        /*
         * print data about scf reference and spinors
         */
        printf("\n");
        printf(" nuclear repulsion energy     %25.16f a.u.\n", enuc);
        printf(" scf energy                   %25.16f a.u.\n", escf);
        printf(" mp2 correlation energy (re)  %25.16f a.u.\n", creal(e_corr_mp2));
        printf(" mp2 correlation energy (im)  %25.16f a.u.\n", cimag(e_corr_mp2));
        printf(" mp2 energy                   %25.16f a.u.\n", escf + creal(e_corr_mp2));
        printf("\n");
        printf(" molecular spinors:\n\n");
        printf("   no  occ      one-el energy, a.u.\n");
        for (int i = 0; i < nspinors; i++) {
            printf("%5d%5d%25.16f\n", i + 1, occ[i], pyscf_eps[i]);
        }
        printf("\n");
    }
    else {
        printf(" conversion of sorted eris to 4d tensor    ");
        double t_begin_eris = abs_time();

        if (arith == CC_ARITH_REAL) {
            int64_t nspinors64 = nspinors;
            eri_ex = calloc(nspinors64 * nspinors64 * nspinors64 * nspinors64, sizeof(double));

            extract_two_electron_integrals("hhhh", eri_ex);

            extract_two_electron_integrals("hphh", eri_ex);
            extract_two_electron_integrals("hhhp", eri_ex);

            extract_two_electron_integrals("hphp", eri_ex);
            extract_two_electron_integrals("hhpp", eri_ex);
            extract_two_electron_integrals("pphh", eri_ex);

            extract_two_electron_integrals("phpp", eri_ex);
            extract_two_electron_integrals("pphp", eri_ex);

            extract_two_electron_integrals("ppppr", eri_ex);
        }
        else {
            int64_t nspinors64 = nspinors;
            eri_ex = calloc(nspinors64 * nspinors64 * nspinors64 * nspinors64, sizeof(double complex));

            extract_two_electron_integrals_complex("hhhh", eri_ex, 0);

            extract_two_electron_integrals_complex("hphh", eri_ex, 0);
            extract_two_electron_integrals_complex("hhhp", eri_ex, 0);

            extract_two_electron_integrals_complex("hphp", eri_ex, 0);
            extract_two_electron_integrals_complex("hhpp", eri_ex, 0);
            extract_two_electron_integrals_complex("pphh", eri_ex, 0);

            extract_two_electron_integrals_complex("phpp", eri_ex, 0);
            extract_two_electron_integrals_complex("pphp", eri_ex, 0);

            extract_two_electron_integrals_complex("ppppr", eri_ex, 1);
        }

        printf("%10.2f\n", abs_time() - t_begin_eris);

        /*
         * Fock matrix
         * for canonical orbitals only
         */
        printf(" fock matrix construction                  ");
        double t_begin_fock = abs_time();

        if (arith == CC_ARITH_REAL) {
            fm_ex = calloc(nspinors * nspinors, sizeof(double));
            for (int i = 0; i < nspinors; i++) {
                ((double *) fm_ex)[nspinors * i + i] = get_eps(artem_to_dirac(nspinors, i));
            }
        }
        else {
            fm_ex = calloc(nspinors * nspinors, sizeof(double complex));
            for (int i = 0; i < nspinors; i++) {
                ((double complex *) fm_ex)[nspinors * i + i] = get_eps(artem_to_dirac(nspinors, i)) + 0.0 * I;
            }
        }

        printf("%10.2f\n", abs_time() - t_begin_fock);
    }

    /*
     * from tensors to trains: fock
     */
    /*printf(" conversion fock -> train                  ");
    double t_begin_conversion_fock = abs_time();
    void *fm = NULL;
    if (arith == CC_ARITH_REAL) {
        fm = d_fm_create(fm_ex, nspinors, nocc);
    }
    else {
        fm = z_fm_create(fm_ex, nspinors, nocc);
    }
    printf("%10.2f\n", abs_time() - t_begin_conversion_fock);*/

    /*
     * from tensors to trains: 2-electron integrals
     */
    /*printf(" conversion eri -> train                   ");
    double t_begin_conversion_eri = abs_time();
    void *eri = NULL;
    if (arith == CC_ARITH_REAL) {
        eri = d_eri_create(eri_ex, nspinors, nocc, tol, tol_eri);
    }
    else {
        eri = z_eri_create(eri_ex, nspinors, nocc, tol, tol_eri);
    }
    printf("%10.2f\n", abs_time() - t_begin_conversion_eri);*/

    /*
     * CCSD iterations
     */
    printf(" initialization of amplitudes              ");
    double t_begin_ccsd_init = abs_time();
    void *ccsd = NULL;

    if (arith == CC_ARITH_REAL) {
        ccsd = d_ccsd_create(eri_ex, fm_ex, nspinors, nocc, tol, cc_opts->tt_options.use_tensor_trains);
    }
    else {
        ccsd = z_ccsd_create(eri_ex, fm_ex, nspinors, nocc, tol, cc_opts->tt_options.use_tensor_trains);
    }
    printf("%10.2f\n", abs_time() - t_begin_ccsd_init);

    /*
     * MP2 energy
     */
    double mp2_energy = 0.0;
    if (arith == CC_ARITH_REAL) {
        mp2_energy = d_mp2_energy(ccsd, eri_ex, fm_ex);
    }
    else {
        mp2_energy = z_mp2_energy(ccsd, eri_ex, fm_ex);
    }
    printf(" MP2 correlation energy = %10.10f\n", mp2_energy);

    //print_all(ccsd, eri, fm);

    printf("\n");
    printf(" solution of amplitude equations (sector 0h0p)\n");
    printf(" ---------------------------------------------------------------------------------------\n");
    printf(" it.                 e(corr)                 e(ccsd)                   delta     t,sec  \n");
    printf(" ---------------------------------------------------------------------------------------\n");

    double ecorr_prev = 0.0;
    double ecorr = 0.0;
    int converged = 0;
    int iter = 0;
    double t_begin_iterations = abs_time();

    for (iter = 1; iter <= cc_opts->maxiter; iter++) {

        double t_begin_iter = abs_time();
        if (arith == CC_ARITH_REAL) {
            d_ccsd_update(ccsd, eri_ex, fm_ex, tol);
            ecorr = d_ccsd_energy(ccsd, eri_ex, fm_ex);
        }
        else {
            z_ccsd_update(ccsd, eri_ex, fm_ex, tol);
            ecorr = z_ccsd_energy(ccsd, eri_ex, fm_ex);
        }
        double t_end_iter = abs_time();

        printf(" %3d%24.16f%24.16f%24.16f%10.2f\n", iter, ecorr, cc_opts->escf + ecorr, ecorr - ecorr_prev,
               t_end_iter - t_begin_iter);

        if (fabs(ecorr - ecorr_prev) < cc_opts->conv_thresh) {
            converged = 1;
            break;
        }

        ecorr_prev = ecorr;
    }

    printf(" ---------------------------------------------------------------------------------------\n");
    if (converged) {
        printf(" converged in %d iterations\n\n", iter);
    }
    else {
        printf(" no convergence after %d interations\n\n", iter);
    }

    printf(" average time per iteration = %.3f sec\n\n", (abs_time() - t_begin_iterations) / iter);

    /*
     * Total CC energy
     */
    print_cc_energy(cc_opts->escf, ecorr);
    write_formatted_heff_file_0h0p(cc_opts->escf + ecorr);
}

#endif // TENSOR_TRAIN

