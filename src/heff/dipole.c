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
 * dipole.c
 *
 * Dipole-length transition moments calculations.
 * Model-space estimation of properties.
 *
 * 2019 Alexander Oleynichenko
 ******************************************************************************/

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "error.h"
#include "io.h"
#include "linalg.h"
#include "mvcoef.h"
#include "slater.h"
#include "spinors.h"
#include "options.h"

double complex get_dm(double complex *mat, int *indices);

void tdms_ground_to_excited_1h1p(double complex **dip_mat,
                                 char *rep_name, size_t ms_size, size_t nroots, slater_det_t *dets,
                                 double complex *eigval, double *energy_cm,
                                 double complex *vl, double complex *vr);


/**
 * Reads property matrix elements (in the basis of molecular spinors)
 */
void read_property_formatted(char *file_name, int nspinors, double complex *matrix, int swap_re_im)
{
    FILE *f;
    int idx1, idx2;
    double re, im;

    f = fopen(file_name, "r");
    if (f == NULL) {
        errquit("read_property_formatted(): file '%s' with property integrals not found", file_name);
    }
    while (fscanf(f, "%d%d%lf%lf", &idx1, &idx2, &re, &im) == 4) {
        if (!swap_re_im) {
            matrix[(idx1 - 1) * nspinors + (idx2 - 1)] = re + im*I;
        }
        else {
            matrix[(idx1 - 1) * nspinors + (idx2 - 1)] = im + re*I; // re <-> im
        }
    }
    fclose(f);
}


/**
 * Construct property matrix in the basis of Slater determinants.
 *
 * @param sect_h FS sector: number of holes
 * @param sect_p FS sector: number of particles
 * @param nspinors
 * @param prp_spinors property matrix in the basis of molecular spinors
 *                    dim = nspinors x nspinors
 * @param prp_slater  property matrix in the basis of Slater determinants
 *                    dim = n_det_bra x n_det_ket
 * @param n_det_bra   number of Slater determinants in the bra vector
 * @param bra_dets    Slater determinants - basis for the bra vector
 * @param n_det_ket   the same, but for ket
 * @param ket_dets    the same, but for ket
 */
void matrix_slater_basis(int sect_h, int sect_p, int nspinors,
        double complex *prp_spinor, double complex *prp_slater,
        slater_det_t *bra_dets, size_t n_det_bra, slater_det_t *ket_dets, size_t n_det_ket)
{
    setup_slater(prp_spinor, (matrix_getter_fun) get_dm, sect_h, sect_p, sect_h, sect_p, 1);

    for (size_t i = 0; i < n_det_bra; i++) {
        for (size_t j = 0; j < n_det_ket; j++) {
            slater_det_t *bra = bra_dets + i;
            slater_det_t *ket = ket_dets + j;
            int bra_vac = is_vacuum_det(bra);
            int ket_vac = is_vacuum_det(ket);
            double complex prp_ij;

            if (bra_vac && ket_vac) {
                prp_ij = 0.0 + 0.0 * I;
            }
            else if (bra_vac) {
                size_t i1 = ket->indices[0];
                size_t i2 = ket->indices[1];
                prp_ij = prp_spinor[i1 * nspinors + i2];
            }
            else if (ket_vac) {
                size_t i1 = bra->indices[0];
                size_t i2 = bra->indices[1];
                prp_ij = prp_spinor[i2 * nspinors + i1];
            }
            else{
                prp_ij = slater(bra, ket);
            }
            prp_slater[i * n_det_ket + j] = prp_ij;
        }
    }
}


/**
 * Model-space estimation of properties
 */
void model_space_property(char *prop_name, int swap_re_im)
{
    int sect_h = cc_opts->sector_h;
    int sect_p = cc_opts->sector_p;

    printf("\n");
    printf("\t*********************************************************\n");
    printf("\t*     MODEL-SPACE ESTIMATION OF PROPERTY: %-8s      *\n", prop_name);
    printf("\t*********************************************************\n");
    printf("\n");

    if (!io_file_exists(prop_name)) {
        printf(" Error: file with property integrals not found\n");
        printf(" Calculation will be skipped\n");
        return;
    }

    int nrep = 0;
    struct mv_block mv_blocks[CC_MAX_NREP];

    read_model_vectors_unformatted(sect_h, sect_p, &nrep, mv_blocks);

    double complex *prp_spinor = xzeros(CC_COMPLEX, nspinors, nspinors);

    // read property matrix elements (in the basis of molecular spinors)
    read_property_formatted(prop_name, nspinors, prp_spinor, swap_re_im);

    for (int irep1 = 0; irep1 < nrep; irep1++){
        for (int irep2 = 0; irep2 < nrep; irep2++){
            struct mv_block *block1 = &mv_blocks[irep1];
            struct mv_block *block2 = &mv_blocks[irep2];
            printf(" < %-4s | prop | %-4s >", block1->rep_name, block2->rep_name);

            double complex *prp_slater = xzeros(CC_COMPLEX, block1->ms_size, block2->ms_size);
            double complex *prop = xzeros(CC_COMPLEX, block1->nroots, block2->nroots);

            // construct property matrix in the basis of Slater determinants
            matrix_slater_basis(sect_h, sect_p, nspinors, prp_spinor, prp_slater,
                    block1->dets, block1->ms_size, block2->dets, block2->ms_size);

            // construct property matrix in the basis of model vectors
            size_t nroots_i = block1->nroots;
            size_t nroots_f = block2->nroots;
            size_t ms_size_1 = block1->ms_size;
            size_t ms_size_2 = block2->ms_size;
            double complex *c_bra = block1->vl;
            double complex *c_ket = block2->vr;

            for (size_t i = 0; i < nroots_i; i++) {
                for (size_t j = 0; j < nroots_f; j++) {
                    double complex prop_if = 0.0 + 0.0 * I;
                    for (size_t a = 0; a < ms_size_2; a++) {
                        for (size_t b = 0; b < ms_size_1; b++) {
                            prop_if += conj(c_bra[ms_size_1 * i + b]) * c_ket[ms_size_2 * j + a] *
                                       prp_slater[ms_size_2 * b + a];
                        }
                    }
                    prop[i * nroots_f + j] = prop_if;
                }
            }

            // print property matrix
            int non_zero = 0;
            for (size_t i = 0; i < nroots_i; i++) {
                for (size_t f = 0; f < nroots_f; f++) {
                    double e_i = block1->energy_cm[i];
                    double e_f = block2->energy_cm[f];

                    double complex prop_if = prop[i * nroots_f + f];
                    if (cabs(prop_if) < 1e-6) {
                        continue;
                    }
                    non_zero++;
                    if (non_zero == 1) {
                        printf("  E1(cm-1)  E2(cm-1)         Re            Im          |prop|\n");
                    }
                    printf("%2d (%-4s) -> %2d (%-4s) %10.2f%10.2f%14.6f%14.6f%14.6f\n",
                           i + 1, block1->rep_name, f + 1, block2->rep_name, e_i, e_f,
                           creal(prop_if), cimag(prop_if), cabs(prop_if));
                }
            }

            if (non_zero == 0) {
                printf("     no allowed transitions\n");
            }
            printf("\n");

            // cleanup
            cc_free(prp_slater);
            cc_free(prop);
        }
    }

    // cleanup
    cc_free(prp_spinor);
    for (size_t irep = 0; irep < nrep; irep++) {
        struct mv_block *b = &mv_blocks[irep];
        cc_free(b->dets);
        cc_free(b->eigval);
        cc_free(b->energy_cm);
        cc_free(b->vl);
        cc_free(b->vr);
    }
}


/*******************************************************************************
 * dipole_length_tdms
 *
 * calculates dipole-length TDMs
 * uses disk files:
 * MVCOEF**   coefficients of model vectors
 * *DIPLEN    matrix elements of the dipole moment operator (electronic part only)
 ******************************************************************************/
void dipole_length_tdms(int sect_h, int sect_p)
{
    char mvcoef_file_name[64];
    double eigval_0;
    double const AU2CM = 219474.6313702;

    printf("\n *** DL-TDM SUBROUTINE ***\n");

    if (!io_file_exists("XDIPLEN") ||
        !io_file_exists("YDIPLEN") ||
        !io_file_exists("ZDIPLEN")) {
        printf(" error: no files with dipole moment integrals (*DIPLEN)\n");
        printf(" DL-TDM calculation will be skipped\n");
        return;
    }

    sprintf(mvcoef_file_name, "MVCOEF%d%d", sect_h, sect_p);
    printf(" MVCOEF file = %s\n", mvcoef_file_name);

    int nrep = 0;
    struct mv_block mv_blocks[64];

    // read model vectors from unformatted file
    int f_mvcoef = io_open(mvcoef_file_name, "r");
    while (1) {
        struct mv_block *curr_block = &mv_blocks[nrep];
        size_t rep_name_len;
        size_t ms_size, nroots;

        io_read(f_mvcoef, &rep_name_len, sizeof(rep_name_len));
        io_read(f_mvcoef, curr_block->rep_name, rep_name_len);
        if (strcmp(curr_block->rep_name, "EOF") == 0) {
            break;
        }

        io_read(f_mvcoef, &ms_size, sizeof(ms_size));
        io_read(f_mvcoef, &nroots, sizeof(nroots));
        curr_block->ms_size = ms_size;
        curr_block->nroots = nroots;

        curr_block->dets = (slater_det_t *) cc_malloc(sizeof(slater_det_t) * ms_size);
        curr_block->eigval = (double complex *) cc_malloc(sizeof(double complex) * nroots);
        curr_block->energy_cm = (double *) cc_malloc(sizeof(double) * nroots);
        curr_block->vl = (double complex *) cc_malloc(sizeof(double complex) * ms_size * nroots);
        curr_block->vr = (double complex *) cc_malloc(sizeof(double complex) * ms_size * nroots);

        io_read(f_mvcoef, curr_block->dets, sizeof(slater_det_t) * ms_size);
        io_read(f_mvcoef, curr_block->eigval, sizeof(double complex) * nroots);
        size_t nb = sizeof(double complex) * nroots * ms_size;
        io_read(f_mvcoef, curr_block->vr, nb);
        io_read(f_mvcoef, curr_block->vl, nb);

        nrep++;
    }
    io_read(f_mvcoef, &eigval_0, sizeof(eigval_0));
    io_close(f_mvcoef);

    // calculate energy of each state wrt ground state (eigval_0)
    for (size_t irep = 0; irep < nrep; irep++){
        struct mv_block *b = &mv_blocks[irep];
        for (size_t i = 0; i < b->nroots; i++){
            b->energy_cm[i] = (creal(b->eigval[i]) - eigval_0) * AU2CM;
        }
    }

    printf(" reading dipole moment integrals in spinor basis:\n");
    double complex *d_spinor[3];
    for (int i = 0; i < 3; i++){
        d_spinor[i] = (double complex *) cc_malloc(sizeof(double complex) * nspinors * nspinors);
        memset(d_spinor[i], 0, sizeof(double complex) * nspinors * nspinors);
    }

    char *file_names[] = {"XDIPLEN", "YDIPLEN", "ZDIPLEN"};
    for (int i = 0; i < 3; i++){
        printf(" %s ... ", file_names[i]);
        FILE *f = fopen(file_names[i], "r");
        int idx1, idx2;
        double re, im;
        while (fscanf(f, "%d%d%lf%lf", &idx1, &idx2, &re, &im) == 4) {
            d_spinor[i][(idx2 - 1) * nspinors + (idx1 - 1)] = re + I * im;  // together with transposition
        }
        fclose(f);
        printf("done\n");
    }

    printf(" NOTE: only model vectors are used to estimate TDMs\n");
    printf(" The DL-TDM technique will be documented in A. Hehn, L. Visscher (to be published)\n");
    printf(" Transition dipole moments:\n");

    for (int irep1 = 0; irep1 < nrep; irep1++){
        for (int irep2 = 0; irep2 < nrep; irep2++){
            struct mv_block *block1 = &mv_blocks[irep1];
            struct mv_block *block2 = &mv_blocks[irep2];
            printf(" -----------------------------------------------------------------------------------------------------------\n");
            printf(" < %4s | d | %4s >", block1->rep_name, block2->rep_name);

            double complex *d_slater[3];
            for (int i = 0; i < 3; i++){
                size_t nbytes = sizeof(double complex) * block1->ms_size * block2->ms_size;
                d_slater[i] = (double complex *) cc_malloc(nbytes);
                memset(d_slater[i], 0, nbytes);   // d_slater ?
            }

            // construct dipole moment components (dx,dy,dz) matrices
            // in the basis of Slater determinants
            for (int icoord = 0; icoord < 3; icoord++){
                setup_slater(d_spinor[icoord], (matrix_getter_fun) get_dm, sect_h, sect_p, sect_h, sect_p, 1);
                size_t ms_size_bra = block1->ms_size;
                size_t ms_size_ket = block2->ms_size;

                for (int i = 0; i < ms_size_bra; i++){
                    for (int j = 0; j < ms_size_ket; j++){
                        slater_det_t *bra = &block1->dets[i];
                        slater_det_t *ket = &block2->dets[j];
                        int bra_vac = is_vacuum_det(bra);
                        int ket_vac = is_vacuum_det(ket);
                        double complex d_ij;

                        if (bra_vac && ket_vac) {
                            d_ij = 0.0 + 0.0 * I;
                        }
                        else if (bra_vac) {
                            size_t i1 = ket->indices[0];
                            size_t i2 = ket->indices[1];
                            d_ij = d_spinor[icoord][i1 * nspinors + i2];
                        }
                        else if (ket_vac) {
                            size_t i1 = bra->indices[0];
                            size_t i2 = bra->indices[1];
                            d_ij = d_spinor[icoord][i2 * nspinors + i1];
                        }
                        else{
                            d_ij = slater(bra, ket);
                        }
                        d_slater[icoord][i * ms_size_ket + j] = d_ij;
                    }
                }
            }

            double complex *tdm[3];
            for (int i = 0; i < 3; i++){
                size_t nbytes = sizeof(double complex) * block1->nroots * block2->nroots;
                tdm[i] = (double complex *) cc_malloc(nbytes);
                memset(tdm[i], 0, nbytes);
            }

            // Slater det-s basis -> basis of model vectors
            for (int icoord = 0; icoord < 3; icoord++){
                size_t nroots_i = block1->nroots;
                size_t nroots_f = block2->nroots;
                size_t ms_size_1 = block1->ms_size;
                size_t ms_size_2 = block2->ms_size;
                double complex *c_bra = block1->vl;
                double complex *c_ket = block2->vr;

                for (size_t i = 0; i < nroots_i; i++){
                    for (size_t j = 0; j < nroots_f; j++){
                        double complex tdm_if = 0.0 + 0.0 * I;
                        for (size_t a = 0; a < ms_size_2; a++){
                            for (size_t b = 0; b < ms_size_1; b++){
                                tdm_if += conj(c_bra[ms_size_1 * i + b]) * c_ket[ms_size_2 * j + a] *
                                          d_slater[icoord][ms_size_2 * b + a];
                            }
                        }
                        tdm[icoord][i * nroots_f + j] = tdm_if;
                    }
                }
            }

            size_t nroots_i = block1->nroots;
            size_t nroots_f = block2->nroots;
            int non_zero = 0;

            // print TDM matrix
            for (size_t i = 0; i < nroots_i; i++){
                for (size_t f = 0; f < nroots_f; f++){
                    double e_i = block1->energy_cm[i];
                    double e_f = block2->energy_cm[f];
                    double d_e = e_f - e_i;
                    if (fabs(d_e) < 1e-3) { // skip transitions between degenerate states
                        continue;
                    }
                    if (irep1 == irep2 && i == f) { // skip diagonal elements
                        continue;
                    }
                    // if all dx,dy,dz == 0 => skip this element
                    double dx = cabs(tdm[0][i * nroots_f + f]);
                    double dy = cabs(tdm[1][i * nroots_f + f]);
                    double dz = cabs(tdm[2][i * nroots_f + f]);
                    if (dx < 1e-6 && dy < 1e-6 && dz < 1e-6) {
                        continue;
                    }
                    non_zero++;
                    if (non_zero == 1) {
                        printf("     E1(cm-1)  E2(cm-1)       |d|       osc str          |dx|        |dy|        |dz|\n");
                    }
                    printf("%2d (%4s) -> %2d (%4s) %10.2f%10.2f",
                           i + 1, block1->rep_name, f + 1, block2->rep_name, e_i, e_f);
                    double d2 = 0.0;
                    for (int icoord = 0; icoord < 3; icoord++){
                        double complex tdm_if = tdm[icoord][i * nroots_f + f];
                        d2 += pow(cabs(tdm_if), 2);
                    }
                    // |d|
                    printf("%12.6f", sqrt(d2));
                    // intensity ? osc str I = 2/3 * |d|^2 * excitation energy
                    // in atomic units
                    double osc_str = 2.0 / 3.0 * d2 * (fabs(d_e) / AU2CM);
                    printf("%12.6f    ", osc_str);
                    // |dx|, |dy|, |dz|
                    for (int icoord = 0; icoord < 3; icoord++){
                        double complex tdm_if = tdm[icoord][i * nroots_f + f];
                        printf("%12.6f", cabs(tdm_if));
                    }
                    printf("\n");
                }
            }

            if (non_zero == 0) {
                printf("     no allowed transitions\n");
            }

            // cleanup
            for (int i = 0; i < 3; i++){
                cc_free(d_slater[i]);
                cc_free(tdm[i]);
            }
        }
    }

    printf(" -----------------------------------------------------------------------------------------------------------\n");

    if (sect_h == 1 && sect_p == 1 && cc_opts->mixed == 0) {
        printf("\n");
        printf("  transition moments: ground state -> excited states\n");
        printf(" ----------------------------------------------------\n");
        printf("\n ene level       |d|       osc str          |dx|        |dy|        |dz|\n");
        for (size_t irep = 0; irep < nrep; irep++){
            struct mv_block *b = &mv_blocks[irep];
            tdms_ground_to_excited_1h1p(d_spinor, b->rep_name, b->ms_size, b->nroots,
                                        b->dets, b->eigval, b->energy_cm, b->vl, b->vr);
        }
        printf("\n");
    }

    // cleanup
    cc_free(d_spinor[0]);
    cc_free(d_spinor[1]);
    cc_free(d_spinor[2]);
    for (size_t irep = 0; irep < nrep; irep++){
        struct mv_block *b = &mv_blocks[irep];
        cc_free(b->dets);
        cc_free(b->eigval);
        cc_free(b->energy_cm);
        cc_free(b->vl);
        cc_free(b->vr);
    }

    printf(" *** END OF DL-TDM SUBROUTINE ***\n\n");
}


/*******************************************************************************
 * tdms_ground_to_excited_1h1p
 *
 * Calculates transition dipole moments 0h0p (ground state) -> 1h1p (excited)
 * via the dipole-length transition dipole moments technique (only model vectors
 * are used for these estimates)
 *
 * |exc> = \sum_i C_i |det_i>
 * <0|d|exc> = \sum_i C_i <0|d|det_i>
 ******************************************************************************/
void tdms_ground_to_excited_1h1p(double complex **dip_mat,
                                 char *rep_name, size_t ms_size, size_t nroots, slater_det_t *dets,
                                 double complex *eigval, double *energy_cm,
                                 double complex *vl, double complex *vr)
{
    int icoord, i, j, a, iroot, idet;
    int n_nonzero = 0;
    double const AU2CM = 219474.6313702;

    printf(" < ground | d | %4s >", rep_name);

    // FSCC effective is not Hermitian => <0|d|model vect> != <model vect|d|0>
    // "bra": <left model vect|d|0> (coef-s must be complex conjugated)
    // "ket": <0|d|right model vect>
    for (iroot = 0; iroot < nroots; iroot++){
        double complex *coef_r = &vr[ms_size * iroot];
        double complex *coef_l = &vl[ms_size * iroot];
        double complex d_bra[3] = {0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I};
        double complex d_ket[3] = {0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I};
        double d2_bra = 0.0, d2_ket = 0.0;
        double d_e = fabs(energy_cm[iroot]) / AU2CM;
        for (icoord = 0; icoord < 3; icoord++){
            for (idet = 0; idet < ms_size; idet++){
                slater_det_t *det = &dets[idet];
                i = det->indices[0];
                a = det->indices[1];
                d_bra[icoord] += conj(coef_l[idet]) * dip_mat[icoord][nspinors * a + i];
                d_ket[icoord] += coef_r[idet] * dip_mat[icoord][nspinors * i + a];
            }
            d2_bra += pow(cabs(d_bra[icoord]), 2);
            d2_ket += pow(cabs(d_ket[icoord]), 2);
        }
        double osc_str_bra = 2.0 / 3.0 * d2_bra * d_e;
        double osc_str_ket = 2.0 / 3.0 * d2_ket * d_e;
        double abs_d_ket = sqrt(d2_ket);
        double abs_d_bra = sqrt(d2_bra);
        if (abs_d_ket >= 1e-6 || abs_d_bra >= 1e-6) {
            n_nonzero++;
            // 0h0p -> 1h1p
            printf("\n%10.2f%12.6f%12.6f%16.6f%12.6f%12.6f\n", energy_cm[iroot],
                   sqrt(d2_ket), osc_str_ket, cabs(d_ket[0]), cabs(d_ket[1]), cabs(d_ket[2]));
            // 1h1p -> 0h0p
            printf("%22.6f%12.6f%16.6f%12.6f%12.6f",
                   sqrt(d2_bra), osc_str_bra, cabs(d_bra[0]), cabs(d_bra[1]), cabs(d_bra[2]));
        }
    }
    if (n_nonzero == 0) {
        printf("   forbidden\n");
    }
    else{
        printf("\n");
    }
}


// "accessor" function
// returns matrix element mat[i,j], indices=[i,j]
double complex get_dm(double complex *mat, int *indices)
{
    size_t i = indices[0];
    size_t j = indices[1];

    double complex s = mat[i * nspinors + j];

    return mat[i * nspinors + j];
}
