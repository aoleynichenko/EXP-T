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
 * Model-space estimation of transition moments.
 */

#include "model_space_tdm.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "cc_properies.h"
#include "io.h"
#include "heff.h"
#include "linalg.h"
#include "../heff/mvcoef.h"
#include "../heff/slater_rules.h"
#include "symmetry.h"
#include "spinors.h"


void model_space_tdms(int sect_h, int sect_p, double complex **dip_mat,
                      char *bra_rep_name, size_t bra_ms_size, slater_det_t *bra_dets, size_t bra_nroots,
                      double complex *bra_vecs, double *bra_energies_cm,
                      char *ket_rep_name, size_t ket_ms_size, slater_det_t *ket_dets, size_t ket_nroots,
                      double complex *ket_vecs, double *ket_energies_cm
);

int read_prop_single_file(int nspinors, char *prop_name, double complex *prop_mat);

int read_prop_two_files(int nspinors, char *file_re, char *file_im, double complex *prop_mat);

void construct_model_space_density_matrix(int sect_h, int sect_p,
                                          int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
                                          int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
                                          double complex *dm, int *dim_dm);

double complex contract_prop_with_dm(int sect_h, int sect_p, size_t dim_dm, double complex *dm, double complex *prp);


/**
 * calculates dipole-length TDMs
 * uses disk files:
 * MVCOEF**   coefficients of model vectors
 * *DIPLEN    matrix elements of the dipole moment operator (electronic part only)
 */
void dipole_length_tdms(int sect_h, int sect_p)
{
    int nspinors = get_num_spinors();

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
    struct mv_block mv_blocks[CC_MAX_NUM_IRREPS];
    char mvcoef_file_name[CC_MAX_NUM_IRREPS];
    sprintf(mvcoef_file_name, "MVCOEF%d%d", sect_h, sect_p);
    printf(" Reading model vectors from the %s file...\n", mvcoef_file_name);
    mvcoef_read_vectors_unformatted(sect_h, sect_p, NULL, &nrep, mv_blocks);

    printf(" reading dipole moment integrals in spinor basis...\n");
    double complex *d_spinor[3];
    d_spinor[0] = z_zeros(nspinors, nspinors);
    d_spinor[1] = z_zeros(nspinors, nspinors);
    d_spinor[2] = z_zeros(nspinors, nspinors);
    read_prop_single_file(nspinors, "XDIPLEN", d_spinor[0]);
    read_prop_single_file(nspinors, "YDIPLEN", d_spinor[1]);
    read_prop_single_file(nspinors, "ZDIPLEN", d_spinor[2]);
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
    printf(
        " -----------------------------------------------------------------------------------------------------------\n");
    printf("                                        Transition dipole moments\n");
    printf("                                             < %dh%dp | d | %dh%dp >\n", sect_h, sect_p, sect_h, sect_p);
    printf(
        " -----------------------------------------------------------------------------------------------------------\n");
    printf("\n");

    for (int irep1 = 0; irep1 < nrep; irep1++) {
        for (int irep2 = 0; irep2 < nrep; irep2++) {
            struct mv_block *block1 = &mv_blocks[irep1];
            struct mv_block *block2 = &mv_blocks[irep2];

            model_space_tdms(sect_h, sect_p, d_spinor,
                             block1->rep_name, block1->ms_size, block1->dets, block1->nroots, block1->vl,
                             block1->energy_cm,
                             block2->rep_name, block2->ms_size, block2->dets, block2->nroots, block2->vr,
                             block2->energy_cm
            );

        }
    }

    if (sect_h == 1 && sect_p == 1) {

        // without restoration of intermediate normalization

        printf("\n\n");
        printf(
            " -----------------------------------------------------------------------------------------------------------\n");
        printf("                           transition moments: ground state -> excited states\n");
        printf("                                             < vac | d | 1h1p>   \n");
        printf(
            " -----------------------------------------------------------------------------------------------------------\n");
        printf("\n");

        slater_det_t vac_det;
        vac_det.indices[0] = 0;
        vac_det.indices[1] = 0;
        double complex vac_coef = 1.0 + 0.0 * I;
        double energy_0 = 0.0;

        for (size_t irep = 0; irep < nrep; irep++) {
            struct mv_block *b = &mv_blocks[irep];

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
        printf(
            " -----------------------------------------------------------------------------------------------------------\n");
        printf("                            transition moments: ground state -> excited states\n");
        printf("                                        < 0h0p+1h1p | d | 0h0p+1h1p>          \n");
        printf("                                    (restored intermediate normalization)     \n");
        printf(
            " -----------------------------------------------------------------------------------------------------------\n");
        printf("\n");

        int nrep_0011;
        struct mv_block mv_blocks_0011[CC_MAX_NUM_IRREPS];

        mvcoef_read_vectors_unformatted(1, 1, "MVCOEF0011", &nrep_0011, mv_blocks_0011);
        struct mv_block *b0011 = &mv_blocks_0011[0];

        // the same irrep
        // ground -> excited
        model_space_tdms(1, 1, d_spinor,
                         b0011->rep_name, b0011->ms_size, b0011->dets, 1, b0011->vl, b0011->energy_cm,
                         b0011->rep_name, b0011->ms_size, b0011->dets, b0011->nroots - 1, b0011->vr + b0011->ms_size,
                         b0011->energy_cm + 1
        );
        // excited -> ground
        model_space_tdms(1, 1, d_spinor,
                         b0011->rep_name, b0011->ms_size, b0011->dets, b0011->nroots - 1, b0011->vl + b0011->ms_size,
                         b0011->energy_cm + 1,
                         b0011->rep_name, b0011->ms_size, b0011->dets, 1, b0011->vr, b0011->energy_cm
        );

        // other irreps
        for (size_t irep = 0; irep < nrep; irep++) {
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
    for (size_t irep = 0; irep < nrep; irep++) {
        mvblock_free(mv_blocks + irep);
    }

    printf("\n");
    printf("\t\t\t\t********************************\n");
    printf("\t\t\t\t**  END OF MS-TDM SUBROUTINE  **\n");
    printf("\t\t\t\t********************************\n");
    printf("\n");
}


void model_space_tdms(
    int sect_h, int sect_p, double complex **dip_mat,
    char *bra_rep_name, size_t bra_ms_size, slater_det_t *bra_dets, size_t bra_nroots,
    double complex *bra_vecs, double *bra_energies_cm,
    char *ket_rep_name, size_t ket_ms_size, slater_det_t *ket_dets, size_t ket_nroots,
    double complex *ket_vecs, double *ket_energies_cm
)
{
    double complex *dm;
    int dim_dm;
    double const AU2CM = 219474.6313702;
    int non_zero = 0;
    int nspinors = get_num_spinors();

    dm = z_zeros(nspinors, nspinors);

    printf(" < %4s | d | %4s >", bra_rep_name, ket_rep_name);

    for (int i = 0; i < bra_nroots; i++) {
        for (int f = 0; f < ket_nroots; f++) {
            double e_i = bra_energies_cm[i];
            double e_f = ket_energies_cm[f];

            if (fabs(e_i - e_f) < 1e-3) { // skip transitions between degenerate states
                continue;
            }

            double complex *coef_i = &bra_vecs[bra_ms_size * i];
            double complex *coef_f = &ket_vecs[ket_ms_size * f];

            construct_model_space_density_matrix(
                sect_h, sect_p,
                bra_ms_size, coef_i, bra_dets,
                ket_ms_size, coef_f, ket_dets,
                dm, &dim_dm
            );
            double complex dx = contract_prop_with_dm(sect_h, sect_p, dim_dm, dm, dip_mat[0]);
            double complex dy = contract_prop_with_dm(sect_h, sect_p, dim_dm, dm, dip_mat[1]);
            double complex dz = contract_prop_with_dm(sect_h, sect_p, dim_dm, dm, dip_mat[2]);
            double abs_dx = cabs(dx);
            double abs_dy = cabs(dy);
            double abs_dz = cabs(dz);
            double d2 = abs_dx * abs_dx + abs_dy * abs_dy + abs_dz * abs_dz;
            double d = sqrt(d2);
            // intensity ? osc str I = 2/3 * |d|^2 * excitation energy
            // in atomic units
            double osc_str = 2.0 / 3.0 * d2 * (fabs(e_f - e_i) / AU2CM);

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


double complex contract_prop_with_dm(int sect_h, int sect_p, size_t dim_dm, double complex *dm, double complex *prp)
{
    int n_active = 0;
    int active_spinors[CC_MAX_SPINORS]; // local -> global spinor index mapping
    int nspinors = get_num_spinors();

    // TODO: refactor, separate utility function (move to spinors.c)
    for (int i = 0; i < nspinors; i++) {
        if ((is_act_hole(i) && sect_h > 0) || (is_act_particle(i) && sect_p > 0)) {
            active_spinors[n_active] = i;
            n_active++;
        }
    }

    double complex sum = 0.0 + 0.0 * I;
    for (int i = 0; i < n_active; i++) {
        for (int j = 0; j < n_active; j++) {
            int p = active_spinors[i];
            int q = active_spinors[j];
            sum += dm[i * n_active + j] * prp[p * nspinors + q];
        }
    }

    return sum;
}
