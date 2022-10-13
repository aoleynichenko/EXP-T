/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2022 The EXP-T developers.
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
 * Direct calculation of norms and overlaps of correlated wavefunctions
 * in different Fock space sector.
 * The formula <psi_I|psi_J> = < phi_I | {e^T^\dagger} {e^T} | phi_J > is used,
 * where |psi> are exact wavefunctions and |phi> are model vectors.
 * Infinite expansion is truncated at quadratic terms (T^2).
 *
 * 2022 Alexander Oleynichenko
 */

#include "cc_wfn_overlap.h"

#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "cc_properies.h"
#include "datamodel.h"
#include "engine.h"
#include "mvcoef.h"
#include "spinors.h"
#include "symmetry.h"
#include "ms_prop.h"


double complex sector_0h0p_overlap(int approximation);

void sector_1h0p_overlap(char *result_name, int approximation);

void sector_0h1p_overlap(char *result_name, int approximation);

void sector_0h2p_overlap(char *result_name, int approximation);


void msprop_transform_slater_to_model(size_t nroots_i, size_t nroots_f,
                                      size_t ms_size_i, size_t ms_size_f,
                                      double complex *coefs_bra, double complex *coefs_ket,
                                      double complex *prop_slater, double complex *prop_model);

void construct_property_matrix_slater(
        int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p,
        size_t n_det_bra, slater_det_t *bra_dets, size_t n_det_ket, slater_det_t *ket_dets,
        double complex *prp_spinor, double complex *prp_slater
);

void print_property_matrix(
        int nroots_i, int nroots_f,
        double *energies_i, double *energies_f,
        char *rep_name_i, char *rep_name_f,
        double complex *prop_matrix);

void conjugate_cluster_amplitudes_0h0p();

void conjugate_cluster_amplitudes_1h0p();

void conjugate_cluster_amplitudes_0h1p();

void conjugate_cluster_amplitudes_0h2p();

void construct_cc_wfn_overlap_matrix(int sect_h, int sect_p, int dim, slater_det_t *det_list,
                                     double complex *overlap_slater, double scalar_term,
                                     int n_overlap_diagrams, char **overlap_diagrams);


void calculate_wavefunction_norms_and_overlaps(int sect_h, int sect_p)
{
    int n_overlap_diagrams = 0;
    char *overlap_diagrams[10];

    if (cc_opts->calc_overlap[sect_h][sect_p] == 0) {
        return;
    }

    /*
     * get complex conjugate for all cluster operators
     */
    conjugate_cluster_amplitudes_0h0p();
    if (sect_h >= 1) {
        conjugate_cluster_amplitudes_1h0p();
    }
    if (sect_p >= 1) {
        conjugate_cluster_amplitudes_0h1p();
    }
    if (sect_p >= 2) {
        conjugate_cluster_amplitudes_0h2p();
    }

    /*
     * scalar part:
     * 1 + fully contracted diagrams without any external lines
     */
    double complex scalar_part = sector_0h0p_overlap(CC_PROPERTIES_APPROX_QUADRATIC);

    /*
     * one-electron part
     */
    if (sect_h >= 1) {
        sector_1h0p_overlap("overlap_1h0p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_1h0p");
    }
    if (sect_p >= 1) {
        sector_0h1p_overlap("overlap_0h1p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_0h1p");
    }

    /*
     * two-electron part
     */
    if (sect_p >= 2) {
        sector_0h2p_overlap("overlap_0h2p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_0h2p");
    }

    printf("\n");
    printf(" **\n");
    printf(" ** Overlap of electronic wavefunctions\n");
    printf(" ** Approximation: quadratic\n");
    printf(" ** Scalar term: %.6f %.6fi\n", creal(scalar_part), cimag(scalar_part));
    printf(" ** Diagrams: ");
    for (int i = 0; i < n_overlap_diagrams; i++) {
        printf("%s ", overlap_diagrams[i]);
    }
    printf("\n");
    printf(" **\n");
    printf("\n");

    /*
     * read model vectors from disk
     */
    int n_irreps = 0;
    struct mv_block mv_blocks[CC_MAX_NUM_IRREPS];
    mvcoef_read_vectors_unformatted(sect_h, sect_p, NULL, &n_irreps, mv_blocks);

    for (int irrep = 0; irrep < n_irreps; irrep++) {

        struct mv_block *block = mv_blocks + irrep;

        size_t ms_size = block->ms_size;
        size_t n_roots = block->nroots;
        char *rep_name = block->rep_name;
        slater_det_t *det_list = block->dets;
        double *energies_cm = block->energy_cm;
        double complex *vectors_left = block->vl;
        double complex *vectors_right = block->vr;

        printf(" Irrep %s\n", block->rep_name);

        double complex *overlap_slater = x_zeros(CC_COMPLEX, ms_size, ms_size);
        double complex *overlap = x_zeros(CC_COMPLEX, n_roots, n_roots);

        // obtain overlap matrix in the basis of Slater determinants
        construct_cc_wfn_overlap_matrix(sect_h, sect_p, ms_size, det_list, overlap_slater,
                                        scalar_part, n_overlap_diagrams, overlap_diagrams);

        // construct property matrix in the basis of model vectors
        msprop_transform_slater_to_model(
                n_roots, n_roots, ms_size, ms_size,
                vectors_right, vectors_right, overlap_slater, overlap
        );

        // print smart table of matrix elements
        print_property_matrix(
                n_roots, n_roots, energies_cm, energies_cm,
                rep_name, rep_name, overlap
        );

        cc_free(overlap_slater);
        cc_free(overlap);
    }

    /*
     * clean up
     */
    for (size_t irep = 0; irep < n_irreps; irep++) {
        mvblock_free(mv_blocks + irep);
    }
}


void construct_cc_wfn_overlap_matrix(int sect_h, int sect_p, int dim, slater_det_t *det_list,
                                     double complex *overlap_slater, double scalar_term,
                                     int n_diagrams, char **list_diagrams)
{
    memset(overlap_slater, 0, sizeof(double complex) * dim * dim);

    for (int i = 0; i < dim; i++) {
        overlap_slater[i * dim + i] = scalar_term;
    }

    for (int k = 0; k < n_diagrams; k++) {
        char *dg_name = list_diagrams[k];
        diagram_t *dg_norm = diagram_stack_find(dg_name);
        setup_slater(dg_norm, (matrix_getter_fun) diagram_get, sect_h, sect_p, sect_h, sect_p, rank(dg_name) / 2);

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                slater_det_t *bra = det_list + i;
                slater_det_t *ket = det_list + j;
                overlap_slater[i * dim + j] += slater_rule(bra, ket);
            }
        }
    }
}


double complex sector_0h0p_overlap(int approximation)
{
    double overlap = 1.0 + 0.0 * I;

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return overlap;
    }

    double t1t1 = scalar_product("C", "N", "t1c", "t1c");
    double t2t2 = 0.25 * scalar_product("C", "N", "t2c", "t2c");
    overlap += t1t1 + t2t2;

    if (triples_enabled()) {
        double t3t3 = (1.0 / 36.0) * scalar_product("C", "N", "t3c", "t3c");
        overlap += t3t3;
    }

    return overlap;
}


void sector_0h1p_overlap(char *result_name, int approximation)
{
    tmplt(result_name, "pp", "11", "12", NOT_PERM_UNIQUE);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    dg_stack_pos_t pos = get_stack_pos();

    // Q1
    reorder("t1c", "r1", "21");
    mult("t1c+", "r1", "r2", 1);
    closed("r2", "r3");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q2
    reorder("t2c", "r1", "3412");
    mult("t2c+", "r1", "r2", 3);
    closed("r2", "r3");
    update(result_name, -0.5, "r3");
    restore_stack_pos(pos);

    // Q3
    reorder("s1c+", "r1", "21");
    mult("s1c", "r1", "r2", 1);
    closed("r2", "r3");
    update(result_name, +1.0, "r3");
    restore_stack_pos(pos);

    // Q4
    reorder("s2c", "r1", "1342");
    mult("r1", "t1c+", "r2", 2);
    closed("r2", "r3");
    update(result_name, +1.0, "r3");
    restore_stack_pos(pos);

    // Q5
    reorder("s2c+", "r1", "1342");
    mult("r1", "t1c", "r2", 2);
    closed("r2", "r3");
    update(result_name, +1.0, "r3");
    restore_stack_pos(pos);

    // Q6
    reorder("s2c+", "r1", "3412");
    mult("s2c", "r1", "r2", 3);
    closed("r2", "r3");
    update(result_name, +0.5, "r3");
    restore_stack_pos(pos);

    if (triples_enabled()) {

        // Q7
        reorder("t3c", "r1", "456123");
        mult("t3c+", "r1", "r2", 5);
        closed("r2", "r3");
        update(result_name, -1.0 / 12.0, "r3");
        restore_stack_pos(pos);

        // Q8
        reorder("s3c", "r1", "145623");
        mult("r1", "t2c+", "r2", 4);
        closed("r2", "r3");
        update(result_name, 0.25, "r3");
        restore_stack_pos(pos);

        // Q9
        reorder("s3c+", "r1", "145623");
        mult("r1", "t2c", "r2", 4);
        closed("r2", "r3");
        update(result_name, 0.25, "r3");
        restore_stack_pos(pos);

        // Q10
        reorder("s3c+", "r1", "456123");
        mult("s3c", "r1", "r2", 5);
        closed("r2", "r3");
        update(result_name, +1.0 / 12.0, "r3");
        restore_stack_pos(pos);
    }
}


void sector_1h0p_overlap(char *result_name, int approximation)
{
    tmplt(result_name, "hh", "11", "12", NOT_PERM_UNIQUE);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    dg_stack_pos_t pos = get_stack_pos();

    // Q1
    reorder("t1c+", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    closed("r2", "r3");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q2
    reorder("t2c+", "r1", "3412");
    mult("t2c", "r1", "r2", 3);
    closed("r2", "r3");
    update(result_name, 0.5, "r3");
    restore_stack_pos(pos);

    // Q3
    reorder("h1c", "r1", "21");
    mult("h1c+", "r1", "r2", 1);
    closed("r2", "r3");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q4
    reorder("h2c", "r1", "1342");
    mult("r1", "t1c+", "r2", 2);
    closed("r2", "r3");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q5
    reorder("h2c+", "r1", "1342");
    mult("r1", "t1c", "r2", 2);
    closed("r2", "r3");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q6
    reorder("h2c", "r1", "3412");
    mult("h2c+", "r1", "r2", 3);
    closed("r2", "r3");
    update(result_name, -0.5, "r3");
    restore_stack_pos(pos);
}


void sector_0h2p_overlap(char *result_name, int approximation)
{
    tmplt(result_name, "pppp", "1111", "1234", NOT_PERM_UNIQUE);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    dg_stack_pos_t pos = get_stack_pos();

    // Q1
    reorder("t2c", "r1", "3412");
    mult("t2c+", "r1", "r2", 2);
    closed("r2", "r3");
    update(result_name, +0.5, "r3");
    restore_stack_pos(pos);

    // Q2
    reorder("s2c", "r1", "1342");
    mult("r1", "t1c+", "r2", 1);
    reorder("r2", "r3", "1423");
    closed("r3", "r4");
    perm("r4", "(12)");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q3
    reorder("t1c", "r1", "21");
    mult("s2c+", "r1", "r2", 1);
    closed("r2", "r3");
    perm("r3", "(34)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q4
    reorder("s2c+", "s2c+_21", "2143"); // interchange electrons 1 <-> 2
    reorder("s2c+_21", "r1", "2431");
    reorder("s2c", "r2", "1324");
    mult("r2", "r1", "r3", 2);
    reorder("r3", "r4", "1324");
    closed("r4", "r5");
    perm("r5", "(12|34)");
    update(result_name, +1.0, "r5");
    restore_stack_pos(pos);

    // Q5
    reorder("s1c+", "r1", "21");
    mult("x2c", "r1", "r2", 1);
    closed("r2", "r3");
    perm("r3", "(34)");
    update(result_name, +1.0, "r3");
    restore_stack_pos(pos);

    // Q6
    reorder("x2c+", "r1", "3412");
    mult("r1", "s1c", "r2", 1);
    reorder("r2", "r3", "3412");
    closed("r3", "r4");
    perm("r4", "(12)");
    update(result_name, +1.0, "r4");
    restore_stack_pos(pos);

    // Q7
    reorder("x2c+", "r1", "3412");
    mult("x2c", "r1", "r2", 2);
    closed("r2", "r3");
    update(result_name, +0.5, "r3");
    restore_stack_pos(pos);

    if (triples_enabled()) {

        // Q8
        reorder("x3c", "r1", "124563");
        mult("r1", "t1c+", "r2", 2);
        closed("r2", "r3");
        update(result_name, +1.0, "r3");
        restore_stack_pos(pos);

        // Q9
        reorder("x3c+", "r1", "124563");
        mult("r1", "t1c", "r2", 2);
        closed("r2", "r3");
        update(result_name, +1.0, "r3");
        restore_stack_pos(pos);

        // Q10
        reorder("s2c+", "r1", "3412");
        reorder("x3c", "r2", "124356");
        mult("r2", "r1", "r3", 3);
        closed("r3", "r4");
        perm("r4", "(34)");
        update(result_name, +0.5, "r4");
        restore_stack_pos(pos);

        // Q11
        reorder("x3c+", "r1", "145623");
        mult("r1", "s2c", "r2", 3);
        reorder("r2", "r3", "1423");
        closed("r3", "r4");
        perm("r4", "(12)");
        update(result_name, +0.5, "r4");
        restore_stack_pos(pos);

        // Q12
        reorder("s3c", "r1", "145623");
        mult("r1", "t2c+", "r2", 3);
        reorder("r2", "r3", "1423");
        closed("r3", "r4");
        perm("r4", "(12)");
        update(result_name, -0.5, "r4");
        restore_stack_pos(pos);

        // Q13
        reorder("s3c+", "r1", "124356");
        reorder("t2c", "r2", "3412");
        mult("r1", "r2", "r3", 3);
        closed("r3", "r4");
        perm("r4", "(34)");
        update(result_name, -0.5, "r4");
        restore_stack_pos(pos);

        // Q14
        reorder("t3c", "r1", "456123");
        mult("t3c+", "r1", "r2", 4);
        closed("r2", "r3");
        update(result_name, 1.0 / 6.0, "r3");
        restore_stack_pos(pos);

        // Q15
        reorder("s3c+", "r1", "245613");
        reorder("s3c", "r2", "152346");
        mult("r2", "r1", "r3", 4);
        reorder("r3", "r4", "1342");
        closed("r4", "r5");
        perm("r5", "(12|34)");
        update(result_name, -0.25, "r5");
        restore_stack_pos(pos);

        // Q16
        reorder("x3c+", "r1", "456123");
        mult("x3c", "r1", "r2", 4);
        closed("r2", "r3");
        update(result_name, 1.0 / 6.0, "r3");
        restore_stack_pos(pos);
    }
}
