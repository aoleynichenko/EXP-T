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

#include "model_vectors.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "eigenvalues.h"
#include "eff_config.h"
#include "options.h"
#include "mvcoef.h"
#include "model_space.h"
#include "symmetry.h"
#include "slater_det.h"


/**
 * Prints model vectors to stdout (formatted output).
 * Also prints effective configurations of active spinors calculated via the
 * quasi-natural spinors approach.
 *
 * If the 'nroots' directive is given in the input file, only the
 * vectors specified by this directive will be printed and analyzed.
 */
void print_model_vectors_stdout(int sector_h, int sector_p, slater_det_t **det_list, size_t *block_dims,
                                double complex **eigvalues, double complex **coef_left, double complex **coef_right)
{
    double const COEF_THRESH = 1e-4;   // threshold for printing model vec-s coeff-s
    int first_irrep = first_nonzero_irrep(block_dims);
    double eff_config[CC_MAX_SPINORS];

    double overall_eff_conf[CC_MAX_SPINORS];
    memset(overall_eff_conf, 0, sizeof(overall_eff_conf));

    if (cc_opts->print_model_vectors == 0 &&
        cc_opts->print_eff_config == 0) {
        return;
    }

    int n_active = 0;
    int active_spinors[CC_MAX_SPINORS];
    get_active_space(sector_h, sector_p, &n_active, active_spinors);

    size_t nroots[CC_MAX_NUM_IRREPS];
    get_nroots(block_dims, eigvalues, nroots);

    printf("\n Sector (%dh,%dp) -- analysis of model vectors (right vectors)\n", sector_h, sector_p);
    printf(" first line : irrep, state number, total energy, eigenvalue\n");
    printf(" other lines: coefficients of contributing determinants (above a threshold of %.1e)\n",
           COEF_THRESH);

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {

        size_t dim = block_dims[irrep];
        size_t nroots_irrep = nroots[irrep];
        slater_det_t *det_basis = det_list[irrep];

        /*
         * print energy for the vacuum state in the special case of the 1h1p sector
         */
        if (sector_h == 1 && sector_p == 1 && irrep == get_vacuum_irrep()) {
            printf("\n Irrep %d (%4s) State %d Energy %23.15f Eigenvalue %14.8f%14.4E\n vacuum state\n",
                    irrep - first_irrep + 1, get_irrep_name(irrep), 0, cc_opts->eref, 0.0, 0.0);
        }

        if (dim == 0 || nroots_irrep == 0) {
            continue;
        }

        // get max number of roots for this irrep
        char *irrep_name = get_irrep_name(irrep);

        for (size_t i = 0; i < nroots_irrep; i++) {
            double complex *left_vector = coef_left[irrep] + dim * i;
            double complex *right_vector = coef_right[irrep] + dim * i;

            print_model_vector(stdout, sector_h, sector_p,
                               irrep - first_irrep + 1, irrep_name, i, eigvalues[irrep][i],
                               dim, coef_right[irrep] + dim * i, det_basis, COEF_THRESH);

            /*
             * calculate effective configurations for the target electronic states
             */
            if (cc_opts->print_eff_config) {

                get_eff_configuration(sector_h, sector_p, dim, det_basis, right_vector, eff_config);

                printf(" effective configuration:\n");
                for (size_t j = 0; j < n_active; j++) {

                    int ispinor = active_spinors[j];

                    // accumulate eff configuration over all electronic states
                    overall_eff_conf[ispinor] += eff_config[j];

                    if (fabs(eff_config[j]) < 1e-2) {
                        continue;
                    }

                    printf("  %12.6f", eff_config[j]);
                    printf("  %-6s #%4d (%12.6f)", get_irrep_name(spinor_info[ispinor].repno), ispinor + 1,
                           spinor_info[ispinor].eps);
                    if (cc_opts->spinor_labels[ispinor] != NULL) {
                        printf("  %s", cc_opts->spinor_labels[ispinor]);
                    }
                    printf("\n");
                }
            }

            /*
             * In case of simple Intermediate Hamiltonian:
             * calculate fraction (%) of main space determinants
             */
            if (cc_opts->intham_imms_opts.sectors[sector_h][sector_p]) {
                double fraction_main = get_fraction_of_main_space_determinants(sector_h, sector_p, dim, det_basis,
                                                                               right_vector);
                printf(" fraction of main space determinants: %.1f%%\n", fraction_main * 100);
            }
        }
    }


    if (cc_opts->print_eff_config) {
        printf("\n\n");
        printf(" sum of effective configurations over the electronic states considered\n");
        printf(" (this output can be used to choose a main model space in IH-IMMS calculations)\n");
        printf("\n");

        for (int i = 0; i < get_num_spinors(); i++) {
            if (!is_active(i)) {
                continue;
            }

            printf("  %12.6f", overall_eff_conf[i]);
            printf("  %-6s  #%4d (%12.6f)", get_irrep_name(spinor_info[i].repno), i + 1, spinor_info[i].eps);
            if (cc_opts->spinor_labels[i] != NULL) {
                printf("  %s", cc_opts->spinor_labels[i]);
            }
            printf("\n");
        }
    }


    printf("\n\n");
}


void write_model_vectors_unformatted(int sector_h, int sector_p, slater_det_t **det_list, size_t *block_dims,
                                     double complex **eigvalues, double complex **coef_left,
                                     double complex **coef_right)
{
    // open unformatted file with model vectors
    int f_mvcoef = mvcoef_open(sector_h, sector_p);

    double lowest_root = get_lowest_eigenvalue(block_dims, eigvalues);

    size_t nroots[CC_MAX_NUM_IRREPS];
    get_nroots(block_dims, eigvalues, nroots);

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {

        size_t dim = block_dims[irrep];
        size_t nroots_irrep = nroots[irrep];

        if (dim == 0 || nroots_irrep == 0) {
            continue;
        }

        slater_det_t *det_basis = det_list[irrep];
        char *irrep_name = get_irrep_name(irrep);

        // print model vectors to the unformatted file MVCOEF**
        mvcoef_write_vectors_unformatted(f_mvcoef, irrep_name, nroots_irrep, dim, det_basis,
                                         eigvalues[irrep], coef_left[irrep], coef_right[irrep]);
    }

    mvcoef_close(f_mvcoef, lowest_root);
}


/**
 * Prints model vectors.
 * coef_thresh -- threshold for printing model vec-s coeff-s
 */
void print_model_vector(FILE *f_out, int sect_h, int sect_p,
                        int rep_no, char *rep_name, int state_no, double complex energy,
                        size_t len, double complex *coeffs, slater_det_t *det_list, double coef_thresh)
{
    fprintf(f_out, "\n");
    fprintf(f_out, " Irrep %d (%4s) State %d Energy %23.15f Eigenvalue %14.8f%14.4E\n",
            rep_no, rep_name, state_no + 1, cc_opts->eref + creal(energy), creal(energy), cimag(energy));

    if (cc_opts->print_model_vectors == 0) {
        return;
    }

    for (size_t j = 0; j < len; j++) {
        double complex coef = coeffs[j];
        if (cabs(coef) < coef_thresh) {
            continue;
        }
        fprintf(f_out, " %10.5f%10.5f ", creal(coef), cimag(coef));
        print_slater_det(f_out, sect_h, sect_p, det_list + j);
    }
}
