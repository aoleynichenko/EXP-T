/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2023 The EXP-T developers.
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
 * Intermediate Hamiltonian (IH) for incomplete main model spaces (IMMS).
 *
 * For details, see:
 * A. V. Zaitsevskii, N. S. Mosyagin, A. V. Oleynichenko, E. Eliav.
 * Generalized relativistic small-core pseudopotentials accounting for
 * quantum electrodynamic effects: construction and pilot applications.
 * Int. J. Quantum Chem. V. 123(8), e27077 (2023)
 */

#include "intham_imms.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "comdef.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"
#include "../heff/model_space.h"


static int spinor_subspaces[CC_MAX_SPINORS];

void intham_imms_print_options(ih_imms_options_t *ih_opts);

void intham_imms_print_determinant_shifts(int sect_h, int sect_p, int shift_type, double frontier_energy);

int intham_imms_is_main_space_det(int sect_h, int sect_p, slater_det_t *det);

void get_model_determinant_for_amplitude(int exc_rank, int *exc_indices, int *valence_labels, slater_det_t *model_det);

double intham_imms_calculate_frontier_energy_auto(int sect_h, int sect_p);

double intham_imms_get_model_det_shift(int sect_h, int sect_p, slater_det_t *det, double frontier_energy);


/**
 * This function is used to calculate energy ranges for the IH-IMMS spinor subspaces.
 * It is activated if only the TOTAL number of spinors in each subspace is given by the user.
 *
 * @param emin_values
 * @param emax_values
 */
void intham_get_emin_emax(double *emin_values, double *emax_values)
{
    /*
     * 1. get energies of ACTIVE spinors
     * 2. and sort these energies
     */
    double *spinor_energies = (double *) calloc(get_num_spinors(), sizeof(double));
    int count = 0;

    for (int i = 0; i < get_num_spinors(); i++) {
        if (is_active(i)) {
            spinor_energies[count++] = get_eps(i);
        }
    }

    qsort(spinor_energies, count, sizeof(double), cmpfunc_double);

    /*
     * for each subspace of active spinors: get lower and upper bound
     */
    ih_imms_options_t *ih_opts = &cc_opts->intham_imms_opts;
    count = 0;
    for (int isub = 0; isub < ih_opts->n_spinor_subspaces; isub++) {

        int nspinors = ih_opts->subspace_total_nspinors[isub];
        double emin = spinor_energies[count] - 1e-8;
        double emax = spinor_energies[count + nspinors - 1] + 1e-8;
        count += nspinors;

        emin_values[isub] = emin;
        emax_values[isub] = emax;
    }

    free(spinor_energies);
}


/**
 * calculates shifts which will be used in the simplest IH-like technique
 */
void intham_imms_setup(int sect_h, int sect_p)
{
    if (!intham_imms_in_sector(sect_h, sect_p)) {
        return;
    }

    ih_imms_options_t *ih_opts = &cc_opts->intham_imms_opts;

    /*
     * special case: subspaces are defined by the total number of spinors in each subspace
     */
    if (ih_opts->subspaces_definition == IH_IMMS_SUBSPACES_DEF_NTOTAL) {
        double emin_values[IH_IMMS_MAX_SPINOR_SUBSPACES];
        double emax_values[IH_IMMS_MAX_SPINOR_SUBSPACES];

        intham_get_emin_emax(emin_values, emax_values);

        for (int isub = 0; isub < ih_opts->n_spinor_subspaces; isub++) {
            ih_opts->subspace_energy_ranges[isub][0] = emin_values[isub];
            ih_opts->subspace_energy_ranges[isub][1] = emax_values[isub];
        }
    }

    /*
     * 'spinor_subspaces' is a global variable
     */
    memset(spinor_subspaces, 0, sizeof(spinor_subspaces));

    for (int i = 0; i < get_num_spinors(); i++) {
        double eps = get_eps(i);

        if (is_active(i)) {

            // to which subspace this spinor belongs?
            for (int isub = 0; isub < ih_opts->n_spinor_subspaces; isub++) {
                double emin = ih_opts->subspace_energy_ranges[isub][0];
                double emax = ih_opts->subspace_energy_ranges[isub][1];
                if (emin <= eps && eps <= emax) {
                    spinor_subspaces[i] = isub + 1;
                }
            }

            // if the subspace was not found in the list => "the remaining subspace"
            if (spinor_subspaces[i] == 0) {
                spinor_subspaces[i] = ih_opts->n_spinor_subspaces + 1;
            }
        }
    }

    /*
     * print settings of the simplest IH
     */
    intham_imms_print_options(ih_opts);

    /*
     * automatic determination of the frontier energy (if needed)
     */
    if (ih_opts->det_shift_auto[sect_h][sect_p]) {
        double frontier_energy = intham_imms_calculate_frontier_energy_auto(sect_h, sect_p);
        ih_opts->det_shift_to[sect_h][sect_p] = frontier_energy;
    }

    /*
     * print table with determinantal shifts
     */
    intham_imms_print_determinant_shifts(sect_h, sect_p, ih_opts->shift_type,
                                         ih_opts->det_shift_to[sect_h][sect_p]);
}


/**
 * Checks if IH-IMMS is enabled in the given FS sector.
 */
int intham_imms_in_sector(int sect_h, int sect_p)
{
    return cc_opts->intham_imms_opts.sectors[sect_h][sect_p];
}


/**
 * prints settings of the simplest intermediate Hamiltonian-like technique.
 */
void intham_imms_print_options(ih_imms_options_t *ih_opts)
{
    printf("\n");
    printf(" Intermediate Hamiltonian-like shifts\n");

    /*
     * list of target sectors
     */
    printf("   target sectors: ");
    for (int h = 0; h < MAX_SECTOR_RANK; h++) {
        for (int p = 0; p < MAX_SECTOR_RANK; p++) {
            if (ih_opts->sectors[h][p]) {
                printf("%dh%dp ", h, p);
            }
        }
    }
    printf("\n");

    /*
     * list spinors of subspaces
     */
    printf("   spinor subspaces:\n");
    printf("    %16s%16s\n", "emin", "emax");
    for (int i = 0; i < ih_opts->n_spinor_subspaces; i++) {
        double emin = ih_opts->subspace_energy_ranges[i][0];
        double emax = ih_opts->subspace_energy_ranges[i][1];
        printf("  %2d%16.8f%16.8f\n", i, emin, emax);
    }

    /*
     * list types of model determinants which are considered as model ones
     */
    printf("   main subspaces of model determinants:\n");
    printf("       emin");
    for (int ispace = 0; ispace < ih_opts->n_spinor_subspaces; ispace++) {
        printf("%16.8f", ih_opts->subspace_energy_ranges[ispace][0]);
    }
    printf("\n");
    printf("       emax");
    for (int ispace = 0; ispace < ih_opts->n_spinor_subspaces; ispace++) {
        printf("%16.8f", ih_opts->subspace_energy_ranges[ispace][1]);
    }
    printf("\n");

    for (int imain = 0; imain < ih_opts->n_main_subspaces; imain++) {
        int npart = 0;
        for (int ispace = 0; ispace < ih_opts->n_spinor_subspaces; ispace++) {
            npart += ih_opts->main_occ[imain][ispace];
        }
        printf("   %2d  %1dh%1dp", imain + 1, 0, npart);
        for (int ispace = 0; ispace < ih_opts->n_spinor_subspaces; ispace++) {
            printf("%16d", ih_opts->main_occ[imain][ispace]);
        }
        printf("\n");
    }

    /*
     * frontier energies used to calculate total "many-electron" shifts S_K
     */
    for (int h = 0; h < MAX_SECTOR_RANK; h++) {
        for (int p = 0; p < MAX_SECTOR_RANK; p++) {
            if (ih_opts->sectors[h][p]) {
                if (ih_opts->det_shift_auto[h][p]) {
                    printf("       %dh%dp  determinant frontier energy = auto\n", h, p);
                }
                else {
                    printf("       %dh%dp  determinant frontier energy = %g a.u.\n",
                           h, p, ih_opts->det_shift_to[h][p]);
                }
            }
        }
    }

    /*
     * shifting formula parameters
     */
    printf("   formula for denominator shifts: ");
    if (ih_opts->shift_type == CC_SHIFT_REAL) {
        printf("real\n");
    }
    else if (ih_opts->shift_type == CC_SHIFT_REALIMAG) {
        printf("realimag\n");
    }
    else if (ih_opts->shift_type == CC_SHIFT_IMAG) {
        printf("imag\n");
    }
    else if (ih_opts->shift_type == CC_SHIFT_TAYLOR) {
        printf("taylor\n");
    }
    else if (ih_opts->shift_type == CC_SHIFT_NONE) {
        printf("none\n");
    }
    else {
        printf("unknown\n");
    }
    printf("   attenuation parameter (npower): %d\n", ih_opts->npower);
    printf("   scale shift by: %g\n", ih_opts->scale_shift);
}


/**
 * automatic determination of frontier energy for the given main model space.
 * returns frontier energy.
 */
double intham_imms_calculate_frontier_energy_auto(int sect_h, int sect_p)
{
    size_t block_dims[CC_MAX_NUM_IRREPS];
    ih_imms_options_t *ih1_opts = &cc_opts->intham_imms_opts;

    double frontier_energy = -1e9;

    slater_det_t **det_basis = construct_model_space(sect_h, sect_p, block_dims);
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        for (int i = 0; i < block_dims[irrep]; i++) {

            slater_det_t *det = det_basis[irrep] + i;
            double det_energy = model_det_energy(sect_h, sect_p, det);

            /*
             * checks if this determinant belongs to any of the subspaces declared as the main ones
             */
            int is_main = 0;
            int det_occ_signature[IH_IMMS_MAX_SPINOR_SUBSPACES];
            memset(det_occ_signature, 0, sizeof(det_occ_signature));

            for (int i = 0; i < sect_h + sect_p; i++) {
                int ispinor = det->indices[i];
                int subspace = spinor_subspaces[ispinor] - 1;
                det_occ_signature[subspace]++;
            }

            for (int imain = 0; imain < ih1_opts->n_main_subspaces; imain++) {
                if (intcmp(ih1_opts->n_spinor_subspaces + 1, det_occ_signature, ih1_opts->main_occ[imain]) == 0) {
                    is_main = 1;
                    break;
                }
            }

            if (is_main && det_energy > frontier_energy) {
                frontier_energy = det_energy;
            }
        }
    }

    /*
     * we add small delta E in order to ensure that the upper main space determinant
     * will be treated as the main one
     */
    return frontier_energy + 1e-8;
}


/**
 * prints beautiful table with shifts for all model space determinants
 */
void intham_imms_print_determinant_shifts(int sect_h, int sect_p, int shift_type, double frontier_energy)
{
    size_t block_dims[CC_MAX_NUM_IRREPS];
    char buf[1024];
    int count = 1;

    printf("   shift for model space determinants\n");
    printf("   frontier energy for this sector: %g\n\n", frontier_energy);

    slater_det_t **det_basis = construct_model_space(sect_h, sect_p, block_dims);
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        for (int i = 0; i < block_dims[irrep]; i++) {

            slater_det_t *det = det_basis[irrep] + i;

            double det_energy = model_det_energy(sect_h, sect_p, det);
            double det_shift = intham_imms_get_model_det_shift(sect_h, sect_p, det, frontier_energy);
            int is_main = intham_imms_is_main_space_det(sect_h, sect_p, det);

            str_print_slater_det(buf, sect_h, sect_p, det);

            printf("   %6d  %s  %-4s    Ene=%12.8f  Shift=%12.8f",
                   count, buf, is_main ? "main" : "int", det_energy, det_shift);

            /*
             * print imaginary unit in the special case of imaginary shifts
             */
            if (shift_type == CC_SHIFT_REALIMAG || shift_type == CC_SHIFT_IMAG) {
                printf("i");
            }

            printf("\n");
            count++;
        }
    }
    printf("\n");
}


/**
 * calculate shift for the given excitation amplitude.
 *
 * the shifting will be applied only to the determinants
 * which are not allowed in the "main space"
 */
double intham_imms_get_shift_value(int sect_h, int sect_p, int rank, int *indices, int *valence_labels)
{
    /*
     * extract model determinant corresponding to the cluster amplitude
     */
    slater_det_t model_det;
    get_model_determinant_for_amplitude(rank, indices, valence_labels, &model_det);

    /*
     * if the determinant belongs to the main space => no shifts
     * else (intermediate space determinant) => calculate shift
     */
    int is_main = intham_imms_is_main_space_det(sect_h, sect_p, &model_det);
    if (is_main) {
        return 0.0;
    }

    /*
     * Shift = R - E(det)
     * E(det)  zero-order energy of a given determinant
     * R       frontier energy
     * TODO: works for particle sectors only!
     */
    double frontier_energy = cc_opts->intham_imms_opts.det_shift_to[sect_h][sect_p];
    double shift = intham_imms_get_model_det_shift(sect_h, sect_p, &model_det, frontier_energy);
    double scaling_factor = cc_opts->intham_imms_opts.scale_shift;

    return scaling_factor * shift;
}


/**
 * checks if the determinant belongs to the main subspace of the model space.
 */
int intham_imms_is_main_space_det(int sect_h, int sect_p, slater_det_t *det)
{
    ih_imms_options_t *ih_opts = &cc_opts->intham_imms_opts;

    if (ih_opts->sectors[sect_h][sect_p] == 0) {
        return 0;
    }

    double det_energy = model_det_energy(sect_h, sect_p, det);
    double frontier_energy = cc_opts->intham_imms_opts.det_shift_to[sect_h][sect_p];

    if (det_energy <= frontier_energy) {
        return 1;
    }

    int det_occ_signature[IH_IMMS_MAX_SPINOR_SUBSPACES];
    memset(det_occ_signature, 0, sizeof(det_occ_signature));

    for (int i = 0; i < sect_h + sect_p; i++) {
        int ispinor = det->indices[i];
        int subspace = spinor_subspaces[ispinor] - 1;
        det_occ_signature[subspace]++;
    }

    for (int imain = 0; imain < ih_opts->n_main_subspaces; imain++) {
        if (intcmp(ih_opts->n_spinor_subspaces + 1, det_occ_signature, ih_opts->main_occ[imain]) == 0) {
            return 1;
        }
    }

    return 0;
}


/**
 * for the given excitation constructs model space determinant
 * on which this excitation acts.
 * hole indices first, then particle indices.
 */
void get_model_determinant_for_amplitude(int exc_rank, int *exc_indices, int *valence_labels, slater_det_t *model_det)
{
    int j = 0;

    // holes
    for (int i = 0; i < exc_rank; i++) {
        if (valence_labels[i] == 1 && is_hole(exc_indices[i])) {
            model_det->indices[j] = exc_indices[i];
            j++;
        }
    }

    // particles
    for (int i = 0; i < exc_rank; i++) {
        if (valence_labels[i] == 1 && is_particle(exc_indices[i])) {
            model_det->indices[j] = exc_indices[i];
            j++;
        }
    }
}


/**
 * calculates shift for the whole determinant:
 * shift = frontier_energy - energy_of_determinant
 */
double intham_imms_get_model_det_shift(int sect_h, int sect_p, slater_det_t *det, double frontier_energy)
{
    double det_energy = model_det_energy(sect_h, sect_p, det);

    // for main model space determinants: shift is ALWAYS zero
    if (intham_imms_is_main_space_det(sect_h, sect_p, det)) {
        return 0.0;
    }

    double shift = 0.0;
    if (det_energy > frontier_energy) {
        shift = frontier_energy - det_energy;
    }

    return shift;
}


/**
 * calculates fraction of main space determinants in the given model vector.
 * to be used with the simple Intermediate Hamiltonian technique.
 */
double get_fraction_of_main_space_determinants(
        int sector_h,
        int sector_p,
        size_t dim,
        slater_det_t *det_list,
        const double complex *model_vector
)
{
    double fraction_main = 0.0;

    for (size_t j = 0; j < dim; j++) {
        if (intham_imms_is_main_space_det(sector_h, sector_p, det_list + j)) {
            double complex coef = model_vector[j];
            fraction_main += cabs(coef) * cabs(coef);
        }
    }

    return fraction_main;
}

