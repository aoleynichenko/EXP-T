/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2021 The EXP-T developers.
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
 * intham1.c
 * =========
 *
 * Simple version of the intermediate Hamiltonian-like technique:
 * (1) active space of spinors is divided into main and intermediate spaces;
 * (2) energies of intermediate spinors are shifted to some upper (for holes)
 *     and lower (for particles) bounds;
 * (3) total S shift is calculated within one of two approaches:
 *     (3.1) S = \sum_{in} s_k - \sum_{out} s_k ("box")
 *     (3.2) S = sum over quasiparticle annihilation lines ("line")
 * (4) shifted denominator is calculated via one of the formulas for dynamic
 *     shifts: real, real simulation of imaginary, imaginary shifts.
 *
 * Input section for this version of "intermediate Hamiltonian" is organized
 * as follows (example):
 *
 * intham1
 *   sectors <list-of-sector-labels>
 *   shift_type [real|realimag|imag|...]
 *   main <int n_main_holes> <int n_main_particles>
 *   shift_to <upper for holes> <lower for particles>
 *   [box|line]
 * end
 *
 * example: IH in the 1h1p sector
 *
 * intham1
 *   sector 1h1p
 *   shift_type real
 *   main 2 6
 *   shift_to -1.5 0.5
 *   box
 * end
 *
 * 2021 Alexander Oleynichenko
 ******************************************************************************/

#include "intham1.h"

#include <stdio.h>
#include <string.h>

#include "comdef.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"

static double orb_shifts[CC_MAX_SPINORS];

void print_spinor_shifts(double *spinor_shifts);


void intham1_get_main_bounds(double *lowest_hole_energy, double *highest_particle_energy)
{
    // scan for the lowest hole energy
    for (int i = 0; i < get_num_spinors(); i++) {
        if (is_active_main(i) && is_hole(i)) {
            *lowest_hole_energy = get_eps(i);
            break;
        }
    }

    // scan for the highest particle energy
    for (int i = get_num_spinors()-1; i >= 0; i--) {
        if (is_active_main(i) && is_particle(i)) {
            *highest_particle_energy = get_eps(i);
            break;
        }
    }
}


/*
 * barycenters are calculated for the main active spinors only
 * (separately for holes and particles).
 */
void intham1_get_main_barycenters(double *hole_barycenter, double *particle_barycenter)
{
    // barycenter of active main holes
    int n_main_holes = 0;
    *hole_barycenter = 0.0;
    for (int i = 0; i < get_num_spinors(); i++) {
        if (is_active_main(i) && is_hole(i)) {
            *hole_barycenter += get_eps(i);
            n_main_holes++;
        }
    }
    if (n_main_holes > 0) {
        *hole_barycenter = *hole_barycenter / n_main_holes;
    }

    // barycenter of active main particles
    int n_main_particles = 0;
    *particle_barycenter = 0.0;
    for (int i = 0; i < get_num_spinors(); i++) {
        if (is_active_main(i) && is_particle(i)) {
            *particle_barycenter += get_eps(i);
            n_main_particles++;
        }
    }
    if (n_main_particles > 0) {
        *particle_barycenter = *particle_barycenter / n_main_particles;
    }
}


void intham1_calculate_shifts(int sect_h, int sect_p)
{
    ih1_options_t *ih1_opts;

    ih1_opts = &cc_opts->ih1_opts;
    if (ih1_opts->sectors[sect_h][sect_p] == 0) {
        return;
    }

    if (ih1_opts->shift_to_formula == IH1_SHIFT_TO_MAIN_BOUNDS) {
        intham1_get_main_bounds(&ih1_opts->holes_shift_to, &ih1_opts->particles_shift_to);
    }
    else if (ih1_opts->shift_to_formula == IH1_SHIFT_TO_BARYCENTER) {
        intham1_get_main_barycenters(&ih1_opts->holes_shift_to, &ih1_opts->particles_shift_to);
    }

    memset(orb_shifts, 0, sizeof(orb_shifts));

    for (int i = 0; i < get_num_spinors(); i++) {
        if (is_active_intermediate(i)) {
            if (is_hole(i)) {
                orb_shifts[i] = ih1_opts->holes_shift_to - get_eps(i);
            }
            else { // particle
                orb_shifts[i] = ih1_opts->particles_shift_to - get_eps(i);
            }
        }
    }

    printf("\n");
    printf(" Intermediate Hamiltonian-like shifts\n");
    printf("   target sectors: ");
    for (int h = 0; h < MAX_SECTOR_RANK; h++) {
        for (int p = 0; p < MAX_SECTOR_RANK; p++) {
            if (ih1_opts->sectors[h][p]) {
                printf("%dh%dp ", h, p);
            }
        }
    }
    printf("\n");
    printf("   main subspace: %d hole spinors, %d particle spinors\n", ih1_opts->main_holes, ih1_opts->main_particles);
    if (ih1_opts->shift_to_formula == IH1_SHIFT_TO_MAIN_BOUNDS) {
        printf("   shift levels: automatically to the bounds of the main subspace\n");
    }
    else if (ih1_opts->shift_to_formula == IH1_SHIFT_TO_BARYCENTER) {
        printf("   shift levels: automatically to the main subspace barycenters\n");
    }
    else {
        printf("   shift levels: custom\n");
    }
    printf("   shift hole levels to: %g a.u.\n", ih1_opts->holes_shift_to);
    printf("   shift particle levels to: %g a.u.\n", ih1_opts->particles_shift_to);
    printf("   formula for denominator shifts: ");
    if (ih1_opts->shift_type == CC_SHIFT_REAL) {
        printf("real\n");
    }
    else if (ih1_opts->shift_type == CC_SHIFT_REALIMAG) {
        printf("realimag\n");
    }
    else if (ih1_opts->shift_type == CC_SHIFT_IMAG) {
        printf("imag\n");
    }
    else if (ih1_opts->shift_type == CC_SHIFT_TAYLOR) {
        printf("taylor\n");
    }
    else if (ih1_opts->shift_type == CC_SHIFT_NONE) {
        printf("none\n");
    }
    else {
        printf("unknown\n");
    }
    printf("   attenuation parameter (npower): %d\n", ih1_opts->npower);
    printf("   formula for summation of shifts: ");
    if (ih1_opts->sum_formula == IH1_SUM_FORMULA_BOX) {
        printf("box\n");
    }
    else if (ih1_opts->sum_formula == IH1_SUM_FORMULA_LINE) {
        printf("line\n");
    }
    else {
        printf("unknown\n");
    }

    // print table with shifts
    print_spinor_shifts(orb_shifts);
}


void print_spinor_shifts(double *spinor_shifts)
{
    printf("   shifts for active spinors:\n");
    printf("   ----------------------------------------------------------------------------------------\n");
    printf("    no    rep         occ        one-el energy          one-el shift        shifted energy\n");
    printf("   ----------------------------------------------------------------------------------------\n");
    for (int i = 0; i < get_num_spinors(); i++) {
        if (!is_active(i)) {
            continue;
        }
        printf("   %4d%4d \"%-6s\"%4d%22.12f%22.12f%22.12f\n",
               i + 1, spinor_info[i].repno,
               get_irrep_name(spinor_info[i].repno), is_hole(i),
               spinor_info[i].eps, spinor_shifts[i], spinor_info[i].eps+spinor_shifts[i]);
    }
    printf("   ----------------------------------------------------------------------------------------\n");

    printf("\n");
}


double intham1_get_shift_value(int sect_h, int sect_p, int rank, int *indices, int *valence_labels)
{
    /*
     * Box:
     * S = \sum_{in} s_i - \sum_{out} s_i
     */
    if (cc_opts->ih1_opts.sum_formula == IH1_SUM_FORMULA_BOX) {
        double shift = 0.0;
        for (int i = 0; i < rank/2; i++) {
            shift += orb_shifts[indices[i]];
            //printf(" + %4d(%10.6f) ", indices[i]+1, orb_shifts[indices[i]]);
        }
        for (int i = rank/2; i < rank; i++) {
            shift -= orb_shifts[indices[i]];
            //printf(" - %4d(%10.6f) ", indices[i]+1, orb_shifts[indices[i]]);
        }
        return shift;
    }
    /*
     * Line:
     * S = \sum_{in,valence} s_i - \sum_{out,valence} s_i
     */
    else if (cc_opts->ih1_opts.sum_formula == IH1_SUM_FORMULA_LINE) {
        double shift = 0.0;
        for (int i = 0; i < rank/2; i++) {
            if (valence_labels[i] == 1) {
                shift += orb_shifts[indices[i]];
                //printf(" + %4d(%10.6f) ", indices[i]+1, orb_shifts[indices[i]]);
            }
        }
        for (int i = rank/2; i < rank; i++) {
            if (valence_labels[i] == 1) {
                shift -= orb_shifts[indices[i]];
                //printf(" - %4d(%10.6f) ", indices[i]+1, orb_shifts[indices[i]]);
            }
        }
        return shift;
    }
}

