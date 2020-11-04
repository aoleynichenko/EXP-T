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
#include "symmetry.h"
#include "spinors.h"
#include "options.h"

double complex get_dm(double complex *mat, int *indices);

void tdms_ground_to_excited_1h1p(double complex **dip_mat,
                                 char *rep_name, size_t ms_size, size_t nroots, slater_det_t *dets,
                                 double complex *eigval, double *energy_cm,
                                 double complex *vl, double complex *vr);

double complex contract_prop_with_dm(int sect_h, int sect_p, size_t dim_dm, double complex *dm, double complex *prp);


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
                prp_ij = 0.0 + 0.0 * I;   // ???????? no!
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

    read_model_vectors_unformatted(sect_h, sect_p, NULL, &nrep, mv_blocks);

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
                                       prp_slater[ms_size_2 * a + b];
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


void construct_ms_density_matrix(int sect_h, int sect_p,
                                 int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
                                 int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
                                 double complex *dm, size_t *dim_dm);

void model_space_tdms(int sect_h, int sect_p, double complex **dip_mat,
                      char *bra_rep_name, size_t bra_ms_size, slater_det_t *bra_dets, size_t bra_nroots, double complex *bra_vecs, double *bra_energies_cm,
                      char *ket_rep_name, size_t ket_ms_size, slater_det_t *ket_dets, size_t ket_nroots, double complex *ket_vecs, double *ket_energies_cm
);


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
    printf("\n");
    printf("\t\t\t\t\t*************************\n");
    printf("\t\t\t\t\t**  MS-TDM SUBROUTINE  **\n");
    printf("\t\t\t\t\t*************************\n");
    printf("\n");

    if (!io_file_exists("XDIPLEN") ||
        !io_file_exists("YDIPLEN") ||
        !io_file_exists("ZDIPLEN")) {
        printf(" error: no files with dipole moment integrals (*DIPLEN)\n");
        printf(" MS-TDM calculation will be skipped\n");
        return;
    }

    int nrep = 0;
    struct mv_block mv_blocks[64];
    char mvcoef_file_name[64];
    sprintf(mvcoef_file_name, "MVCOEF%d%d", sect_h, sect_p);
    printf(" Reading model vectors from the %s file...\n", mvcoef_file_name);
    read_model_vectors_unformatted(sect_h, sect_p, NULL, &nrep, mv_blocks);

    printf(" reading dipole moment integrals in spinor basis...\n");
    double complex *d_spinor[3];
    d_spinor[0] = zzeros(nspinors, nspinors);
    d_spinor[1] = zzeros(nspinors, nspinors);
    d_spinor[2] = zzeros(nspinors, nspinors);
    read_property_formatted("XDIPLEN", nspinors, d_spinor[0], 0);
    read_property_formatted("YDIPLEN", nspinors, d_spinor[1], 0);
    read_property_formatted("ZDIPLEN", nspinors, d_spinor[2], 0);
    printf(" done\n");

    printf(" NOTE: only model vectors are used to estimate TDMs\n");
    printf(" When using the MS-TDM technique, please, cite the following papers:\n");
    printf("  [1] A. Hehn, L. Visscher (to be published)\n");
    printf("  [2] A. Zaitsevskii, A. V. Oleynichenko, E. Eliav.\n");
    printf("      Finite-field calculations of transition properties by the Fock space\n");
    printf("      relativistic coupled cluster method: transitions between different\n");
    printf("      Fock space sectors.\n");
    printf("      Preprints 2020, 2020100442, doi: 10.20944/preprints202010.0442.v1\n");
    printf("\n");
    printf(" -----------------------------------------------------------------------------------------------------------\n");
    printf("                                        Transition dipole moments\n");
    printf("                                             < %dh%dp | d | %dh%dp >\n", sect_h, sect_p, sect_h, sect_p);
    printf(" -----------------------------------------------------------------------------------------------------------\n");
    printf("\n");

    for (int irep1 = 0; irep1 < nrep; irep1++) {
        for (int irep2 = 0; irep2 < nrep; irep2++) {
            struct mv_block *block1 = &mv_blocks[irep1];
            struct mv_block *block2 = &mv_blocks[irep2];

            model_space_tdms(sect_h, sect_p, d_spinor,
                    block1->rep_name, block1->ms_size, block1->dets, block1->nroots, block1->vl, block1->energy_cm,
                    block2->rep_name, block2->ms_size, block2->dets, block2->nroots, block2->vr, block2->energy_cm
            );

        }
    }

    if (sect_h == 1 && sect_p == 1 && cc_opts->mixed == 0) {

        // without restoration of intermediate normalization

        printf("\n\n");
        printf(" -----------------------------------------------------------------------------------------------------------\n");
        printf("                           transition moments: ground state -> excited states\n");
        printf("                                             < vac | d | 1h1p>   \n");
        printf(" -----------------------------------------------------------------------------------------------------------\n");
        printf("\n");

        slater_det_t vac_det;
        vac_det.indices[0] = 0;
        vac_det.indices[1] = 0;
        double complex vac_coef = 1.0 + 0.0*I;
        double energy_0 = 0.0;

        for (size_t irep = 0; irep < nrep; irep++){
            struct mv_block *b = &mv_blocks[irep];

            /*tdms_ground_to_excited_1h1p(d_spinor,
                                             b->rep_name, b->ms_size, b->nroots, b->dets,
                                             b->eigval, b->energy_cm,
                                             b->vl, b->vr);*/

            model_space_tdms(1, 1, d_spinor,
                             get_irrep_name(get_vacuum_irrep()), 1, &vac_det, 1, &vac_coef, &energy_0,
                             b->rep_name, b->ms_size, b->dets, b->nroots, b->vr, b->energy_cm
            );

            model_space_tdms(1, 1, d_spinor,
                             b->rep_name, b->ms_size, b->dets, b->nroots, b->vl, b->energy_cm,
                             get_irrep_name(get_vacuum_irrep()), 1, &vac_det, 1, &vac_coef, &energy_0
            );
        }

        // with transformed vectors, 0h0p+1h1p are mixed

        printf("\n\n");
        printf(" -----------------------------------------------------------------------------------------------------------\n");
        printf("                            transition moments: ground state -> excited states\n");
        printf("                                        < 0h0p+1h1p | d | 0h0p+1h1p>          \n");
        printf("                                    (restored intermediate normalization)     \n");
        printf(" -----------------------------------------------------------------------------------------------------------\n");
        printf("\n");


        int nrep_0011;
        struct mv_block mv_blocks_0011[64];

        read_model_vectors_unformatted(1, 1, "MVCOEF0011", &nrep_0011, mv_blocks_0011);
        struct mv_block *b0011 = &mv_blocks_0011[0];

        // the same irrep
        // ground -> excited
        model_space_tdms(1, 1, d_spinor,
                         b0011->rep_name, b0011->ms_size, b0011->dets, 1, b0011->vl, b0011->energy_cm,
                         b0011->rep_name, b0011->ms_size, b0011->dets, b0011->nroots-1, b0011->vr+b0011->ms_size, b0011->energy_cm+1
        );
        // excited -> ground
        model_space_tdms(1, 1, d_spinor,
                         b0011->rep_name, b0011->ms_size, b0011->dets, b0011->nroots-1, b0011->vl+b0011->ms_size, b0011->energy_cm+1,
                         b0011->rep_name, b0011->ms_size, b0011->dets, 1, b0011->vr, b0011->energy_cm
        );

        // other irreps
        for (size_t irep = 0; irep < nrep; irep++){
            struct mv_block *b = &mv_blocks[irep];
            if (strcmp(get_irrep_name(get_vacuum_irrep()), b->rep_name) == 0) {
                continue;
            }
            // ground -> excited
            model_space_tdms(1, 1, d_spinor,
                             b0011->rep_name, b0011->ms_size, b0011->dets, 1, b0011->vl, b0011->energy_cm,
                             b->rep_name, b->ms_size, b->dets, b->nroots, b->vr, b->energy_cm
            );
            // excited -> ground
            model_space_tdms(1, 1, d_spinor,
                             b->rep_name, b->ms_size, b->dets, b->nroots, b->vl, b->energy_cm,
                             b0011->rep_name, b0011->ms_size, b0011->dets, 1, b0011->vr, b0011->energy_cm
            );
        }
    }

    // TODO: destroy mvcoef block

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

    printf("\n");
    printf("\t\t\t\t********************************\n");
    printf("\t\t\t\t**  END OF MS-TDM SUBROUTINE  **\n");
    printf("\t\t\t\t********************************\n");
    printf("\n");
}


void model_space_tdms(int sect_h, int sect_p, double complex **dip_mat,
        char *bra_rep_name, size_t bra_ms_size, slater_det_t *bra_dets, size_t bra_nroots, double complex *bra_vecs, double *bra_energies_cm,
        char *ket_rep_name, size_t ket_ms_size, slater_det_t *ket_dets, size_t ket_nroots, double complex *ket_vecs, double *ket_energies_cm
        )
{
    double complex *dm;
    size_t dim_dm;
    double const AU2CM = 219474.6313702;
    int non_zero = 0;

    dm = zzeros(nspinors, nspinors);

    printf(" < %4s | d | %4s >", bra_rep_name, ket_rep_name);

    for (int i = 0; i < bra_nroots; i++) {
        for (int f = 0; f < ket_nroots; f++) {
            double e_i = bra_energies_cm[i];
            double e_f = ket_energies_cm[f];

            if (fabs(e_i - e_f) < 1e-3) { // skip transitions between degenerate states
                continue;
            }
            /*if ((strcmp(bra_rep_name, ket_rep_name) == 0) && i == f) { // skip diagonal elements
                continue;
            }*/

            double complex *coef_i = &bra_vecs[bra_ms_size * i];
            double complex *coef_f = &ket_vecs[ket_ms_size * f];

            construct_ms_density_matrix(sect_h, sect_p, bra_ms_size, coef_i, bra_dets, ket_ms_size, coef_f, ket_dets, dm, &dim_dm);
            double complex dx = contract_prop_with_dm(sect_h, sect_p, dim_dm, dm, dip_mat[0]);
            double complex dy = contract_prop_with_dm(sect_h, sect_p, dim_dm, dm, dip_mat[1]);
            double complex dz = contract_prop_with_dm(sect_h, sect_p, dim_dm, dm, dip_mat[2]);
            double abs_dx = cabs(dx);
            double abs_dy = cabs(dy);
            double abs_dz = cabs(dz);
            double d2 = abs_dx*abs_dx + abs_dy*abs_dy + abs_dz*abs_dz;
            double d = sqrt(d2);
            // intensity ? osc str I = 2/3 * |d|^2 * excitation energy
            // in atomic units
            double osc_str = 2.0 / 3.0 * d2 * (fabs(e_f - e_i) / AU2CM);
            //printf("d = %f    %d\n", d, (d < 1e-6));

            if (d > 1e-6) {
                non_zero++;
                if (non_zero == 1) {
                    printf("     E1(cm-1)  E2(cm-1)       |d|       osc str          |dx|        |dy|        |dz|\n");
                }
                printf("%2d (%4s) -> %2d (%4s) %10.2f%10.2f", i + 1, bra_rep_name, f + 1, ket_rep_name, e_i, e_f);
                printf("%12.6f%12.6f    %12.6f%12.6f%12.6f\n", d, osc_str, abs_dx, abs_dy, abs_dz);
            }
        }
    }

    if (non_zero == 0) {
        printf("     forbidden\n");
    }
    printf("\n");

    cc_free(dm);
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
        //printf("\n");
        for (icoord = 0; icoord < 3; icoord++){
            for (idet = 0; idet < ms_size; idet++){
                slater_det_t *det = &dets[idet];
                i = det->indices[0];
                a = det->indices[1];
                //double complex c = coef_l[idet];
                //double complex d = dip_mat[icoord][nspinors * a + i];
                d_bra[icoord] += conj(coef_l[idet]) * dip_mat[icoord][nspinors * i + a];
                d_ket[icoord] += coef_r[idet] * dip_mat[icoord][nspinors * a + i];
                if (icoord == 0) {
                    //printf("+ (%f,%f) x (%f,%f) = (%f,%f)\n", creal(c), cimag(c), creal(d), cimag(d), creal(c*d), cimag(c*d));
                }
            }
            if (icoord == 0) {
                //printf(" d  = (%f,%f)\n", creal(d_bra[icoord]), cimag(d_bra[icoord]));
                //printf("|d| = %f\n", cabs(d_bra[icoord]));
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
