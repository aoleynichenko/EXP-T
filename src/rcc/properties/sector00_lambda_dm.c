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
 * Analytic density matrices in the 0h0p Fock space sector.
 *
 * Diagrams for lambda equations and density matrix can be found in:
 * I. Shavitt, R. J. Bartlett, Many-Body Methods in Chemistry and Physics,
 * Cambridge University Press (2009)
 *
 * For details on lambda equations in the relativistic framework, see:
 * A. Shee, L. Visscher, T. Saue.
 * Analytic one-electron properties at the 4-component relativistic
 * coupled cluster level with inclusion of spin-orbit coupling.
 * J. Chem. Phys. 145, 184107 (2016)
 * doi: 10.1063/1.4966643
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "engine.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"


void sector_0h0p_analytic_density_matrix_block_hh();

void sector_0h0p_analytic_density_matrix_block_pp();

void sector_0h0p_analytic_density_matrix_block_hp();

void sector_0h0p_analytic_density_matrix_block_ph();

double complex sector_0h0p_calculate_overlap(int nspinors, double complex *dm);

void save_density_matrix(char *path, int nspinors, double complex *dm);

void diagram_clear_off_diagonal(char *diagram_name);

void calculate_natural_spinors_and_occ_numbers(char *path_to_dm_file, char *path_to_nat_spinors_file);


/**
 * constructs analytic one-particle density matrix in the 0h0p Fock space sector
 */
void sector_0h0p_analytic_density_matrix_lambda()
{
    double t_start = abs_time();
    printf(" begin construction of density matrix\n");

    int nspinors = get_num_spinors();
    double complex *dm = z_zeros(nspinors, nspinors);
    memset(dm, 0, sizeof(double complex) * nspinors * nspinors);

    // intermediate diagram - "cross"
    /*tmplt("lambda_t", "pp", "11", "12", NOT_PERM_UNIQUE);
    dg_stack_pos_t pos = get_stack_pos();
    reorder("L01_2c", "r1", "3412");
    mult("s2c", "r1", "r2", 3);
    update("lambda_t", 0.5, "r2");
    restore_stack_pos(pos);*/

    /*
     * construct matrices for the hp, hh and pp blocks of the DM
     */
    sector_0h0p_analytic_density_matrix_block_hp();   // => "dm00_hp"
    sector_0h0p_analytic_density_matrix_block_ph();   // => "dm00_ph"
    sector_0h0p_analytic_density_matrix_block_hh();   // => "dm00_hh"
    sector_0h0p_analytic_density_matrix_block_pp();   // => "dm00_pp"

    /*
     * extract matrix elements from diagrams and save them into matrix
     */
    diagram_t *dg_dm00_hp = diagram_stack_find("dm00_hp");
    diagram_t *dg_dm00_ph = diagram_stack_find("dm00_ph");
    diagram_t *dg_dm00_hh = diagram_stack_find("dm00_hh");
    diagram_t *dg_dm00_pp = diagram_stack_find("dm00_pp");

    for (int p = 0; p < nspinors; p++) {
        for (int q = 0; q < nspinors; q++) {
            double complex gamma_pq = 0.0 + 0.0 * I;

            if (is_hole(p) && is_hole(q)) { // h-h
                gamma_pq = (p == q) + diagram_get_2(dg_dm00_hh, p, q);
            }
            else if (is_hole(p) && is_particle(q)) { // h-p
                gamma_pq = diagram_get_2(dg_dm00_hp, p, q);
            }
            else if (is_particle(p) && is_hole(q)) { // p-h
                gamma_pq = diagram_get_2(dg_dm00_ph, p, q);
            }
            else { // p-p
                gamma_pq = diagram_get_2(dg_dm00_pp, p, q);
            }

            dm[p * nspinors + q] = gamma_pq;
        }
    }

    /*
     * save density matrix to disk
     */
    save_density_matrix("density_0h0p.txt", nspinors, dm);

    /*
     * calculate overlap integral <psi|psi> and wavefunction norm
     */
    double complex trace = ztrace(nspinors, dm);
    double complex norm2 = sector_0h0p_calculate_overlap(nspinors, dm);
    double norm = sqrt(creal(norm2));

    printf(" density matrix constructed in %.2f seconds\n", abs_time() - t_start);
    printf(" density matrix is written to file 'density_0h0p.txt'\n");
    printf(" sum of diagonal elements   = %20.12f %20.12f\n", creal(trace), cimag(trace));
    printf(" overlap integral <psi|psi> = %20.12f %20.12f\n", creal(norm2), cimag(norm2));
    printf(" norm                       = %20.12f\n", norm);
    printf("\n");

    calculate_natural_spinors_and_occ_numbers("density_0h0p.txt", "natural_spinors_0h0p.txt");

    // clean up
    cc_free(dm);
}


/**
 * Hole-particle block of the one-particle density matrix.
 */
void sector_0h0p_analytic_density_matrix_block_hp()
{
    tmplt("dm00_hp", "hp", "00", "12", NOT_PERM_UNIQUE);

    dg_stack_pos_t pos = get_stack_pos();

    // DM-1
    update("dm00_hp", 1.0, "t1c");
    restore_stack_pos(pos);

    // DM-2
    reorder("t2c", "r1", "3142");
    mult("r1", "L00_1c", "r2", 2);
    reorder("r2", "r3", "21");
    update("dm00_hp", 1.0, "r3");
    restore_stack_pos(pos);

    // DM-3
    reorder("t1c", "r1", "21");
    mult("r1", "L00_1c", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    update("dm00_hp", -1.0, "r3");
    restore_stack_pos(pos);

    // DM-5
    reorder("L00_2c", "r1", "3412");
    reorder("t1c", "r3", "21");
    mult("t2c", "r1", "r2", 3);
    mult("r2", "r3", "r4", 1);
    update("dm00_hp", -0.5, "r4");
    restore_stack_pos(pos);

    // DM-6
    reorder("t2c", "r1", "3412");
    mult("r1", "L00_2c", "r2", 3);
    mult("t1c", "r2", "r3", 1);
    update("dm00_hp", -0.5, "r3");
    restore_stack_pos(pos);

    /*
     * contribution from perturbative triples (CCSD(T))
     * (diagram DM-7)
     */
    if (cc_opts->cc_model == CC_MODEL_CCSD_T3) {
        diagram_conjugate("t2c", "t2c+");
        reorder("t3c", "r1", "145623");
        mult("r1", "t2c+", "r2", 4);
        update("dm00_hp", 0.25, "r2");
        restore_stack_pos(pos);
    }

    if (cc_opts->sector_h == 0 && cc_opts->sector_p == 1) {
        // DM_01_2
        copy("dm1", "cross");
        tmplt("cross_pv", "pp", "01", "12", NOT_PERM_UNIQUE);
        expand_diagram("cross", "cross_pv");
        reorder("s2c", "r1", "2431");
        mult("r1", "cross_pv", "r2", 2);
        update("dm00_hp", 1.0, "r2");
        restore_stack_pos(pos);

        // DM_01_5
        /*reorder("L01_2c", "r1", "3412");
        mult("t2c", "r1", "r2", 3);
        tmplt("r3", "hp", "00", "12", NOT_PERM_UNIQUE);
        expand_diagram("r2", "r3");
        update("dm00_hp", -0.5, "r3");
        restore_stack_pos(pos);*/
    }
}


/**
 * Particle-hole block of the one-particle density matrix.
 */
void sector_0h0p_analytic_density_matrix_block_ph()
{
    tmplt("dm00_ph", "ph", "00", "12", NOT_PERM_UNIQUE);

    // DM-1'
    update("dm00_ph", 1.0, "L00_1c");
}


/**
 * Hole-hole block of the one-particle density matrix.
 */
void sector_0h0p_analytic_density_matrix_block_hh()
{
    tmplt("dm00_hh", "hh", "00", "12", NOT_PERM_UNIQUE);

    dg_stack_pos_t pos = get_stack_pos();

    // DM-20
    reorder("L00_1c", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    update("dm00_hh", -1.0, "r2");
    restore_stack_pos(pos);

    // DM-21
    reorder("L00_2c", "r1", "3412");
    mult("t2c", "r1", "r2", 3);
    update("dm00_hh", -0.5, "r2");
    restore_stack_pos(pos);

    /*
     * contribution from perturbative triples (CCSD(T))
     * (diagram DM-25)
     */
    if (cc_opts->cc_model == CC_MODEL_CCSD_T3) {

        reorder("t3c+", "r1", "456123");
        mult("t3c", "r1", "r2", 5);
        diagram_clear_off_diagonal("r2");
        update("dm00_hh", -1.0 / 12, "r2");
        restore_stack_pos(pos);

        reorder("t3(d)+", "r1", "456123");
        mult("t3c", "r1", "r2", 5);
        diagram_clear_off_diagonal("r2");
        update("dm00_hh", -1.0 / 12, "r2");
        restore_stack_pos(pos);
    }

    if (cc_opts->sector_h == 0 && cc_opts->sector_p == 1) {
        // DM_01_3
        /*reorder("L01_2c", "r1", "4123");
        reorder("s2c", "r2", "2341");
        mult("r2", "r1", "r3", 3);
        update("dm00_hh", 0.5, "r3");
        restore_stack_pos(pos);*/
    }
}


/**
 * Particle-particle block of the one-particle density matrix.
 */
void sector_0h0p_analytic_density_matrix_block_pp()
{
    tmplt("dm00_pp", "pp", "00", "12", NOT_PERM_UNIQUE);

    dg_stack_pos_t pos = get_stack_pos();

    // DM-20
    reorder("t1c", "r1", "21");
    mult("L00_1c", "r1", "r2", 1);
    update("dm00_pp", 1.0, "r2");
    restore_stack_pos(pos);

    // DM-21
    reorder("t2c", "r1", "3412");
    mult("L00_2c", "r1", "r2", 3);
    update("dm00_pp", 0.5, "r2");
    restore_stack_pos(pos);

    /*
     * contribution from perturbative triples (CCSD(T))
     * (diagram DM-25)
     */
    if (cc_opts->cc_model == CC_MODEL_CCSD_T3) {

        reorder("t3c", "r1", "456123");
        mult("t3c+", "r1", "r2", 5);
        diagram_clear_off_diagonal("r2");
        update("dm00_pp", +1.0 / 12, "r2");
        restore_stack_pos(pos);

        reorder("t3c", "r1", "456123");
        mult("t3(d)+", "r1", "r2", 5);
        diagram_clear_off_diagonal("r2");
        update("dm00_pp", +1.0 / 12, "r2");
        restore_stack_pos(pos);
    }

    if (cc_opts->sector_h == 0 && cc_opts->sector_p == 1) {
        // DM_01_1
        copy("dm1", "cross");
        tmplt("cross_pp", "pp", "00", "12", NOT_PERM_UNIQUE);
        expand_diagram("cross", "cross_pp");
        update("dm00_pp", 1.0, "cross_pp");
        restore_stack_pos(pos);

        // DM_01_4
        /*reorder("L01_2c", "r1", "2341");
        reorder("s2c", "r2", "4123");
        mult("r1", "r2", "r3", 3);
        update("dm00_pp", -1.0, "r3");
        restore_stack_pos(pos);*/
    }
}


/**
 * sets off-diagonal matrix elements of the 1-particle diagram to zero
 */
void diagram_clear_off_diagonal(char *diagram_name)
{
    diagram_t *diag = diagram_stack_find(diagram_name);
    if (diag == NULL) {
        errquit("diagram_clear_off_diagonal(): diagram '%s' not found");
    }
    if (diag->rank != 2) {
        summary(diagram_name);
        errquit("diagram_clear_off_diagonal(): only rank 2 diagrams are allowed");
    }

    for (int p = 0; p < get_num_spinors(); p++) {
        for (int q = 0; q < get_num_spinors(); q++) {
            if (p != q) {
                int idx[2];
                idx[0] = p;
                idx[1] = q;

                diagram_set(diag, 0.0 + 0.0 * I, idx);
            }
        }
    }
}


/**
 * Prints density matrix to the formatted file
 */
void save_density_matrix(char *path, int nspinors, double complex *dm)
{
    FILE *f = fopen(path, "w");
    if (f == NULL) {
        errquit("cannot write density matrix to disk");
    }

    for (int i = 0; i < nspinors; i++) {
        for (int j = 0; j < nspinors; j++) {
            fprintf(f, "%10d%10d%25.16E%25.16E\n", i + 1, j + 1, creal(dm[i * nspinors + j]),
                    cimag(dm[i * nspinors + j]));
        }
    }

    fclose(f);
}

