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

#include "eigenvalues.h"

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "codata.h"
#include "memory.h"
#include "options.h"
#include "symmetry.h"
#include "utils.h"

typedef struct {
    double complex eigval;
    int repno;
    double percent_main;
} eigval_t;

int eigval_cmp(const void *aa, const void *bb);


double max_energy_of_required_roots(size_t *block_dims, double complex **eigvalues)
{
    double max_energy = -1.0e12;

    size_t nroots[CC_MAX_NUM_IRREPS];
    get_nroots(block_dims, eigvalues, nroots);

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {

        if (block_dims[irrep] == 0 || nroots[irrep] == 0) {
            continue;
        }

        // get max number of roots for this irrep

        double complex *ev = eigvalues[irrep];
        if (nroots[irrep] > 0) {
            for (size_t i = 0; i < nroots[irrep]; i++) {
                if (creal(ev[i]) > max_energy) {
                    max_energy = creal(ev[i]);
                }
            }
        }
    }

    return max_energy;
}


double get_lowest_eigenvalue(const size_t *block_dims, double complex **eigvalues)
{
    double min_eig_value = 1e9;

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        if (block_dims[irrep] == 0) {
            continue;
        }
        if (creal(eigvalues[irrep][0]) < min_eig_value) {
            min_eig_value = creal(eigvalues[irrep][0]);
        }
    }

    /*
     * special case of the 1h1p sector:
     * ground state corresponds to the vacuum determinant
     */
    if (cc_opts->sector_h == 1 && cc_opts->sector_p == 1) {
        min_eig_value = 0.0;
    }

    /*
     * special case of the 1h2p sector:
     * ground state can belong to the 0h1p sector
     */
    if (cc_opts->sector_h == 1 && cc_opts->sector_p == 2) {
        if (cc_opts->ground_energy_0h1p < min_eig_value) {
            min_eig_value = cc_opts->ground_energy_0h1p;
        }
    }

    /*
     * special case of the 2h1p sector:
     * ground state can belong to the 1h0p sector
     */
    if (cc_opts->sector_h == 2 && cc_opts->sector_p == 1) {
        if (cc_opts->ground_energy_1h0p < min_eig_value) {
            min_eig_value = cc_opts->ground_energy_1h0p;
        }
    }

    return min_eig_value;
}


void print_eigenvalues_table(int sect_h, int sect_p, size_t *block_dims, slater_det_t **det_list,
                             double complex **eigvalues, double degen_thresh, double complex **coef_right)
{
    /*
     * Allocate table for eigenvalues
     */
    size_t total_heff_dim = 0;
    for (int i = 0; i < get_num_irreps(); i++) {
        total_heff_dim += block_dims[i];
    }
    if (sect_h == 1 && sect_p == 1) {
        total_heff_dim += 1;  // for the reference energy
    }
    eigval_t *eigenvalues = (eigval_t *) cc_malloc(total_heff_dim * sizeof(eigval_t));
    size_t n_eigenvalues = 0;

    /*
     * Copy eigenvalues to the table.
     * Vacuum det is always 'main'.
     */
    if (sect_h == 1 && sect_p == 1) {
        eigenvalues[0].eigval = 0.0;
        eigenvalues[0].repno = get_vacuum_irrep();
        eigenvalues[0].percent_main = 100.0;
        n_eigenvalues++;
    }
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        size_t block_dim = block_dims[irrep];
        for (size_t i = 0; i < block_dim; i++) {
            eigenvalues[n_eigenvalues + i].eigval = eigvalues[irrep][i];
            eigenvalues[n_eigenvalues + i].repno = irrep;

            /*
             * In case of simple Intermediate Hamiltonian:
             * calculate fraction (%) of main space determinants
             */
            double complex *right_vector = coef_right[irrep] + block_dim * i;
            if (cc_opts->intham_imms_opts.sectors[sect_h][sect_p]) {
                double fraction_main = get_fraction_of_main_space_determinants(sect_h, sect_p, block_dim,
                                                                               det_list[irrep],
                                                                               right_vector);
                eigenvalues[n_eigenvalues + i].percent_main = fraction_main * 100;
            }
            else {
                eigenvalues[n_eigenvalues + i].percent_main = 0.0;
            }
        }
        n_eigenvalues += block_dim;
    }

    /*
     * Sort table in ascending order
     */
    qsort(eigenvalues, n_eigenvalues, sizeof(eigval_t), eigval_cmp);

    /*
     * Print table with energy levels
     */
    double complex e0 = eigenvalues[0].eigval;
    double reference_energy = cc_opts->eref;   // (see options.h)
    double max_energy = max_energy_of_required_roots(block_dims, eigvalues);
    if (sect_h == 0 && sect_p == 1) {
        cc_opts->ground_energy_0h1p = creal(e0);
    }
    if (sect_h == 0 && sect_p == 2) {
        cc_opts->ground_energy_0h2p = creal(e0);
    }
    if (sect_h == 0 && sect_p == 3) {
        cc_opts->ground_energy_0h3p = creal(e0);
    }
    if (sect_h == 1 && sect_p == 0) {
        cc_opts->ground_energy_1h0p = creal(e0);
    }
    if (sect_h == 2 && sect_p == 0) {
        cc_opts->ground_energy_2h0p = creal(e0);
    }
    if (sect_h == 3 && sect_p == 0) {
        cc_opts->ground_energy_3h0p = creal(e0);
    }
    if (sect_h == 1 && sect_p == 2) {
        e0 = cc_opts->ground_energy_0h1p + 0.0 * I;
    }
    if (sect_h == 2 && sect_p == 1) {
        e0 = cc_opts->ground_energy_1h0p + 0.0 * I;
    }

    int intham = 0;
    if (cc_opts->intham_imms_opts.sectors[sect_h][sect_p]) {
        intham = 1;
    }

    printf("\n Heff eigenvalues:\n (degeneracy threshold = %.1e a.u.)\n\n", degen_thresh);

    if (intham == 0) {
        printf(" Level  Re(eigenvalue)  Im(eigv)               Abs energy  Rel eigenvalue    Rel eigv, eV  Rel eigv, cm-1  deg  symmetry\n");
        printf(" ------------------------------------------------------------------------------------------------------------------------\n");
    }
    else {
        printf(" Level  Re(eigenvalue)  Im(eigv)               Abs energy  Rel eigenvalue    Rel eigv, eV  Rel eigv, cm-1  %% main  deg  symmetry\n");
        printf(" --------------------------------------------------------------------------------------------------------------------------------\n");
    }

    int ilevel = 1;
    for (size_t i = 0; i < n_eigenvalues; i++) {

        // only energy levels with eigenvalue Ei <= max_energy will be printed
        if (creal(eigenvalues[i].eigval) > max_energy + degen_thresh) {
            break;
        }

        // calculate degeneracy
        int deg = 0;
        for (size_t j = i; j < n_eigenvalues; j++) {
            if (cabs(eigenvalues[i].eigval - eigenvalues[j].eigval) > degen_thresh) {
                break;
            }
            else {
                deg++;
            }
        }

        // print info about (maybe degenerate) level
        double complex ei = eigenvalues[i].eigval;
        double abs_energy = reference_energy + creal(ei);
        double rel_energy = creal(ei) - creal(e0);
        double rel_energy_cm = rel_energy * CODATA_AU_TO_CM;
        double rel_energy_ev = rel_energy * CODATA_AU_TO_EV;

        printf("@%5d%16.10f%10.2e%25.17f%16.10f%16.10f%16.6f  ",
               ilevel, creal(ei), cimag(ei), abs_energy, rel_energy,
               rel_energy_ev, rel_energy_cm);

        // percent of main model space determinants (in case of IH-FSCC)
        if (intham) {
            printf("%6.1f  ", eigenvalues[i].percent_main);
        }

        printf("%2d  ", deg);
        for (size_t irrep = 0; irrep < get_num_irreps(); irrep++) {
            int rep_deg = 0;
            for (size_t j = 0; j < deg; j++) {
                if (eigenvalues[i + j].repno == irrep) {
                    rep_deg++;
                }
            }
            if (rep_deg == 0) { continue; }
            printf(" %s", get_irrep_name(irrep));
            if (rep_deg > 1) {
                printf("(%d)", rep_deg);
            }
        }

        printf("\n");

        // go to the next (maybe degenerate) energy level
        i += deg - 1;
        ilevel++;
    }

    printf("\n");

    cc_free(eigenvalues);
}


/* calculates how many roots in each irrep should be processed.
 * number of roots is given by the two options in the input file:
 * 'nroots' - number of lowest roots for each irrep
 * 'roots_cutoff' - lowest roots are limited by some max energy
 * finally, number of roots = min(nroots,roots_cutoff)
 */
void get_nroots(size_t *block_dims, double complex **eigvalues, size_t *nroots)
{
    double lowest_root = get_lowest_eigenvalue(block_dims, eigvalues);

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {

        size_t nroots_irrep = block_dims[irrep];
        size_t nroots_cutoff = block_dims[irrep];

        // from the 'nroots' option
        if (cc_opts->nroots_specified) {
            char *irrep_name = get_irrep_name(irrep);
            nroots_irrep = get_nroots_for_irrep(irrep_name);
        }
        else if (cc_opts->roots_cutoff_specified) {
            nroots_cutoff = get_nroots_from_cutoff(nroots_irrep, eigvalues[irrep], lowest_root, cc_opts->roots_cutoff);
        }

        nroots[irrep] = MIN(nroots_irrep, nroots_cutoff);
    }
}


size_t get_nroots_from_cutoff(size_t irrep_dim, double complex *eigvalues, double emin, double cutoff)
{
    size_t nroots = 0;

    for (size_t i = 0; i < irrep_dim; i++) {
        double e_i = creal(eigvalues[i]);
        if (e_i - emin <= cutoff) {
            nroots += 1;
        }
    }

    return nroots;
}


int get_nroots_for_irrep(char *irrep_name)
{
    int nroots_irep = 0;

    for (int ii = 0; ii < get_num_irreps(); ii++) {
        if (strcmp(irrep_name, cc_opts->nroots_specs.rep_names[ii]) == 0) {
            nroots_irep = cc_opts->nroots_specs.dim[ii];
            break;
        }
    }

    return nroots_irep;
}


size_t first_nonzero_irrep(const size_t *block_dims)
{
    for (size_t i = 0; i < get_num_irreps(); i++) {
        if (block_dims[i] != 0) {
            return i;
        }
    }

    return 0;
}


/**
 * Comparator for eigenvalues (to be used in qsort)
 */
int eigval_cmp(const void *aa, const void *bb)
{
    const eigval_t *a = aa, *b = bb;
    return (creal(a->eigval) < creal(b->eigval)) ? -1 : (creal(a->eigval) > creal(b->eigval));
}
