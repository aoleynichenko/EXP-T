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
 * intham.c
 * ========
 *
 * Intermediate Hamiltonian implementation.
 *
 * 2021 Alexander Oleynichenko
 ******************************************************************************/

#include "intham.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "codata.h"
#include "datamodel.h"
#include "heff.h"
#include "model_space.h"
#include "mvcoef.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"

static const double ZERO_THRESH = 1e-14;

void construct_int_space_projector(size_t *dims, double complex **vec, double complex **P_i);
void construct_int_space_biorth_projector(size_t *dims, double complex **left_vec, double complex **right_vec, double complex **P_i);
void get_eff_configuration(int sect_h, int sect_p, int ms_size, slater_det_t *det_list,
                           double complex *coef_left, double complex *coef_right, double *config);
void ih_diveps_block_rank_4(block_t *block, double *spinor_shifts, size_t *block_dims, slater_det_t **det_basis,
                            double complex **P_i);
void ih_spinor_shifts(int sect_h, int sect_p, size_t *block_dims, slater_det_t **det_basis,
                      double complex **eigvalues, double complex **left_vectors,
                      double complex **right_vectors, double *spinor_shifts);
void print_spinor_shifts(double *spinor_shifts);
double upper_bound_by_occ(int sect_h, int sect_p, int n_active, int *active_spinors, double *max_eff_config);
void ih_diveps(char *name, double *spinor_shifts, size_t *block_dims, slater_det_t **det_basis,
               double complex **P_i);
void mult_by_shift_rank_4(block_t *src_block, block_t *dst_block, double *spinor_shifts);
void mult_by_shift(char *src_name, char *dst_name, double *spinor_shifts);
double complex get_P_ij(int i1, int i2, int j1, int j2, size_t *block_dims, slater_det_t **det_basis, double complex **P_i);
void matrix_to_tensor(char *tensor_name, size_t *block_dims, slater_det_t **det_basis, double complex **matrix);
void get_off_diagonal(size_t dim, double complex *src, double complex *off_diag);


void intham(int sect_h, int sect_p)
{
    double shifts[CC_MAX_SPINORS];

    printf("\nIntermediate Hamiltonian\n");

    // construct and diagonalize effective Hamiltonian, obtain model vectors

    double complex *heff[CC_MAX_NUM_IRREPS];
    double complex *eigvalues[CC_MAX_NUM_IRREPS];
    double complex *coef_left[CC_MAX_NUM_IRREPS];
    double complex *coef_right[CC_MAX_NUM_IRREPS];
    size_t block_dims[CC_MAX_NUM_IRREPS];

    /*
     * Construct model space
     */
    slater_det_t **det_basis = construct_model_space(sect_h, sect_p, block_dims);
    int n_diagrams = 2;
    char *diagram_names[] = {"veff01", "veff02", NULL};
    construct_heff(sect_h, sect_p, det_basis, block_dims, heff, n_diagrams, diagram_names);

    /*
     * Diagonalization
     */
    diagonalize_heff(sect_h, sect_p, block_dims, heff, eigvalues, coef_left, coef_right);

    /*
     * Calculate shifts for each spinor
     */
    ih_spinor_shifts(sect_h, sect_p, block_dims, det_basis, eigvalues, coef_left, coef_right, shifts);

    /*
     * Construct projector onto intermediate states
     */
    double complex *P_i[CC_MAX_NUM_IRREPS];
    double complex *P_i_offdiag[CC_MAX_NUM_IRREPS];
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        size_t dim = block_dims[irrep];
        if (dim == 0) {
            P_i[irrep] = NULL;
            P_i_offdiag[irrep] = NULL;
        }
        else {
            P_i[irrep] = z_zeros(dim, dim);
            P_i_offdiag[irrep] = z_zeros(dim, dim);
        }
    }
    //construct_int_space_biorth_projector(block_dims, coef_left, coef_right, P_i);
    construct_int_space_projector(block_dims, coef_right, P_i);
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        get_off_diagonal(block_dims[irrep], P_i[irrep], P_i_offdiag[irrep]);
    }

    /*for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        if (P_i[irrep] == NULL) continue;
        //cc_free(P_i[irrep]);
        for (int i = 0; i < block_dims[irrep]; i++) {
            print_slater_det(stdout, 0, 2, det_basis[irrep]+i);
        }
        xprimat(CC_COMPLEX, P_i[irrep], block_dims[irrep], block_dims[irrep], "P_i");
        xprimat(CC_COMPLEX, P_i_offdiag[irrep], block_dims[irrep], block_dims[irrep], "P_i_offdiag");
    }*/

    /*printf("Diagonal of P_i:\n");
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        printf("Irrep %s\n", get_irrep_name(irrep));
        for (int i = 0; i < block_dims[irrep]; i++) {
            printf("%15.8e%15.8e\n", creal(P_i[irrep][i*block_dims[irrep]+i]), cimag(P_i[irrep][i*block_dims[irrep]+i]));
        }
    }*/


    dg_stack_pos_t pos = get_stack_pos();
    tmplt("P_i", "pppp", "1111", "1234", IS_PERM_UNIQUE);
    matrix_to_tensor("P_i", block_dims, det_basis, P_i_offdiag);

    mult_by_shift("x2c", "tksk", shifts);

    reorder("x2c", "r1", "3412");
    mult("P_i", "r1", "r2", 2);
    update("x2nw", -1.0, "r2");

    restore_stack_pos(pos);

    ih_diveps("x2nw", shifts, block_dims, det_basis, P_i);


    /*
     * Cleanup
     */
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        if (block_dims[irrep] == 0) {
            continue;
        }
        cc_free(heff[irrep]);
        cc_free(eigvalues[irrep]);
        cc_free(coef_left[irrep]);
        cc_free(coef_right[irrep]);
        cc_free(det_basis[irrep]);
    }
    cc_free(det_basis);
}


int get_intham_main_subspace_for_irrep(char *irrep_name)
{
    int nroots_irep = 0;

    for (int ii = 0; ii < get_num_irreps(); ii++) {
        if (strcmp(irrep_name, cc_opts->intham_params.main_space.rep_names[ii]) == 0) {
            nroots_irep = cc_opts->intham_params.main_space.dim[ii];
            break;
        }
    }

    return nroots_irep;
}


/*
 * Calculates one-electron (spinor) shifts.
 */
void ih_spinor_shifts(int sect_h, int sect_p, size_t *block_dims, slater_det_t **det_basis,
                      double complex **eigvalues, double complex **left_vectors,
                      double complex **right_vectors, double *spinor_shifts)
{
    // construct list of indices of active spinors
    int n_active = 0;
    int active_spinors[CC_MAX_SPINORS]; // local -> global spinor index mapping
    get_active_space(sect_h, sect_p, &n_active, active_spinors);
    double *max_eff_config = d_zeros(n_active, 1);

    // print list of target (main) states
    printf("Main (target) states:\n");

    printf("                  ");
    for (int j = 0; j < n_active; j++) {
        int ispinor = active_spinors[j];
        printf("%10d", ispinor+1);
    }
    printf("\n");
    printf("                  ");
    for (int j = 0; j < n_active; j++) {
        int ispinor = active_spinors[j];
        char *spinor_irrep = get_irrep_name(spinor_info[ispinor].repno);
        printf("%10s", spinor_irrep);
    }
    printf("\n");
    printf("                  ");
    for (int j = 0; j < n_active; j++) {
        int ispinor = active_spinors[j];
        double eps = spinor_info[ispinor].eps;
        printf("%10.3f", eps);
    }
    printf("\n");

    double lowest_eigvalue = get_lowest_eigenvalue(block_dims, eigvalues);

    for (size_t irep = 0; irep < get_num_irreps(); irep++) {

        char *irrep_name = get_irrep_name(irep);
        int n_main_states = get_intham_main_subspace_for_irrep(irrep_name);
        size_t ms_size = block_dims[irep];
        if (n_main_states == 0) {
            continue;
        }

        printf("  Irrep %s\n", irrep_name);
        for (int i = 0; i < n_main_states; i++) {
            double e_i = CODATA_AU_TO_CM * (eigvalues[irep][i] - lowest_eigvalue);
            printf("    [%d] %10.3f", i+1, e_i);

            double config[CC_MAX_SPINORS];
            double complex *left_vector = left_vectors[irep] + ms_size * i;
            double complex *right_vector = right_vectors[irep] + ms_size * i;
            get_eff_configuration(sect_h, sect_p, ms_size, det_basis[irep], left_vector, right_vector, config);
            for (int j = 0; j < n_active; j++) {
                if (config[j] > max_eff_config[j]) {
                    max_eff_config[j] = config[j];
                }
            }

            // print configuration
            for (int j = 0; j < n_active; j++) {
                printf("%10.6f", config[j]);
            }
            printf("\n");
        }
    }

    printf("\nmax eff config =  ");
    for (int j = 0; j < n_active; j++) {
        printf("%10.6f", max_eff_config[j]);
    }
    printf("\n");

    double upper_bound = upper_bound_by_occ(sect_h, sect_p, n_active, active_spinors, max_eff_config);
    printf("upper bound of spinor energy = %6f\n", upper_bound);

    for (int i = 0; i < get_num_spinors(); i++) {
        if (!is_active(i)) {
            spinor_shifts[i] = 0.0;
            continue;
        }
        double eps = spinor_info[i].eps;
        if (eps <= upper_bound + 1e-10) {
            spinor_shifts[i] = 0.0;
        }
        else {
            spinor_shifts[i] = upper_bound - eps;
        }
    }

    print_spinor_shifts(spinor_shifts);

    cc_free(max_eff_config);
}


double upper_bound_by_occ(int sect_h, int sect_p, int n_active, int *active_spinors, double *max_eff_config)
{
    const double OCC_THRESH = 0.1;

    double upper_bound = -1e10;
    for (int i = 0; i < n_active; i++) {
        int ispinor = active_spinors[i];
        double eps = spinor_info[ispinor].eps;
        if (eps > upper_bound && max_eff_config[i] > OCC_THRESH) {
            upper_bound = eps;
        }
    }

    return upper_bound;
}


/**
 * Constructs projector onto intermediate states in the basis of Slater determinants.
 * Projector is stored as a blocked square matrix, number of blocks is equal to
 * the number of irreps.
 *
 * Projector P_i is calculated using the formula:
 * P_i = P - P_m
 * where
 * P = identity matrix
 * P_m = projector onto main (target) states
 *
 * @param dims  dimensions of the model space for each irrep
 * @param vec   model vectors
 * @param P_i   (output) projector onto intermediate states, array of square matrices.
 *              (must be pre-allocated)
 */
void construct_int_space_projector(size_t *dims, double complex **vec, double complex **P_i)
{
    //printf("Construct P_i\n");
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        char *irrep_name = get_irrep_name(irrep);
        size_t dim = dims[irrep];
        if (dim == 0) {
            continue;
        }

        // P = identity
        double complex *P = z_identity(dim);

        //xprimat(CC_COMPLEX, P, dim, dim, "Identity");

        // P_m (will be stored in P_i)
        int n_main = get_intham_main_subspace_for_irrep(irrep_name);
        double complex *P_m = z_zeros(dim, dim);
        construct_projector(dims[irrep], n_main, vec[irrep], P_m);

        // P_i = P - P_m
        for (size_t i = 0; i < dim*dim; i++) {
            P_i[irrep][i] = P[i] - P_m[i];
        }
        //xaxpy(CC_COMPLEX, dim*dim, -1.0, P_i[irrep], P);

        //xprimat(CC_COMPLEX, P_i[irrep], dim, dim, "P_i[irrep]");

        cc_free(P);
        cc_free(P_m);
    }
    //printf("End of construct P_i\n");
    //exit(0);
}


void construct_int_space_biorth_projector(size_t *dims, double complex **left_vec, double complex **right_vec, double complex **P_i)
{
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        char *irrep_name = get_irrep_name(irrep);
        size_t dim = dims[irrep];
        if (dim == 0) {
            continue;
        }

        // P = identity
        double complex *P = z_identity(dim);

        // P_m (will be stored in P_i)
        int n_main = get_intham_main_subspace_for_irrep(irrep_name);
        double complex *P_m = z_zeros(dim, dim);
        construct_biorth_projector(dims[irrep], n_main, left_vec[irrep], right_vec[irrep], P_m);

        // P_i = P - P_m
        for (size_t i = 0; i < dim*dim; i++) {
            P_i[irrep][i] = P[i] - P_m[i];
        }

        cc_free(P);
        cc_free(P_m);
    }
}



/**
 * Division of the diagram by the energy denominators.
 * Shifted denominators are used, shifts are multiplied by the diagonal matrix elements
 * of the projector onto the intermediate states:
 * D_K^L' = D_K^L + S_K^L * Pi_K;K
 */
void ih_diveps(char *name, double *spinor_shifts, size_t *block_dims, slater_det_t **det_basis,
               double complex **P_i)
{
    diagram_t *dg;


    /*for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        if (block_dims[irrep] == 0) {
            continue;
        }
        for (int i = 0; i < block_dims[irrep]; i++) {
            print_slater_det(stdout, 0, 2, det_basis[irrep]+i);
        }
        xprimat(CC_COMPLEX, P_i[irrep], block_dims[irrep], block_dims[irrep], "P_i");
    }*/

    timer_new_entry("intham_diveps", "IH energy denominators");
    timer_start("intham_diveps");

    dg = diagram_stack_find(name);
    if (dg == NULL) {
        errquit("intham_diveps(): diagram '%s' not found", name);
    }

    for (size_t isb = 0; isb < dg->n_blocks; isb++) {
        block_t *block = dg->blocks[isb];
        symblock_load(block);
        if (block->is_unique == 0) {
            continue;
        }

        if (dg->rank == 4) {
            ih_diveps_block_rank_4(block, spinor_shifts, block_dims, det_basis, P_i);
        }

        /*else if (rank == 4) {
            diveps_block_rank_4(block);
        }
        else if (rank == 6) {
            diveps_block_rank_6(block);
        }
        else {
            diveps_block_general(block);
        }*/

        symblock_store(block);
    }

    timer_stop("intham_diveps");
}


double complex get_P_kk(int idx_0, int idx_1, size_t *block_dims, slater_det_t **det_basis,
                        double complex **P_i);


void ih_diveps_block_rank_4(block_t *block, double *spinor_shifts, size_t *block_dims, slater_det_t **det_basis,
                            double complex **P_i)
{
    int nspinors = get_num_spinors();
    int sect_h = cc_opts->curr_sector_h;
    int sect_p = cc_opts->curr_sector_p;

    int dims_0 = block->indices[0][0];
    int dims_1 = block->indices[1][0];
    int dims_2 = block->indices[2][0];
    int dims_3 = block->indices[3][0];
    int coef_0 = dims_1 * dims_2 * dims_3;
    int coef_1 = dims_2 * dims_3;
    int coef_2 = dims_3;
    //int coef_3 = 1;
    int *block_indices_0 = block->indices[0] + 1; // +1 to skip first element (length)
    int *block_indices_1 = block->indices[1] + 1;
    int *block_indices_2 = block->indices[2] + 1;
    int *block_indices_3 = block->indices[3] + 1;

    double *eps = (double *) cc_malloc(nspinors * sizeof(double));
    get_spinor_energies(nspinors, eps);

    // loop over matrix elements
    for (size_t i = 0; i < block->size; i++) {
        // if matrix element is zero we can skip it
        double complex t = carith ? block->buf[i] : ((double *) block->buf)[i] + 0.0 * I;
        if (cabs(t) < ZERO_THRESH) {
            continue;
        }
        // calculate compound index
        int index = i;
        int idx_0 = index / coef_0;
        index = index % coef_0;
        int idx_1 = index / coef_1;
        index = index % coef_1;
        int idx_2 = index / coef_2;
        index = index % coef_2;
        int idx_3 = index;    // since coef_3 == 1

        // relative index to spinor indices
        idx_0 = block_indices_0[idx_0];
        idx_1 = block_indices_1[idx_1];
        idx_2 = block_indices_2[idx_2];
        idx_3 = block_indices_3[idx_3];

        // calculate denominator and divide matrix element by it
        double denom = eps[idx_0] + eps[idx_1] - eps[idx_2] - eps[idx_3];
        double shift = spinor_shifts[idx_0] + spinor_shifts[idx_1] - spinor_shifts[idx_2] - spinor_shifts[idx_3];
        double complex P_kk = get_P_ij(idx_0, idx_1, idx_0, idx_1, block_dims, det_basis, P_i);
        //printf("%d %d  P_kk = %f %f\n", idx_0, idx_1, creal(P_kk), cimag(P_kk));
        t /= (denom + shift * creal(P_kk));

        if (carith) {
            block->buf[i] = t;
        }
        else {
            ((double *) block->buf)[i] = creal(t);
        }
    }

    cc_free(eps);
}


double complex get_P_kk(int idx_0, int idx_1, size_t *block_dims, slater_det_t **det_basis,
                         double complex **P_i)
{
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        size_t dim = block_dims[irrep];
        for (size_t i = 0; i < dim; i++) {
            slater_det_t *d = det_basis[irrep] + i;
            if ((d->indices[0] == idx_0 && d->indices[1] == idx_1) ||
                (d->indices[0] == idx_1 && d->indices[1] == idx_0)) {
                return P_i[irrep][dim * i + i];
            }
        }
    }

    return 0.0 + 0.0*I;
}


int find_determinant(int idx_0, int idx_1, size_t dim, slater_det_t *det_list)
{
    for (int i = 0; i < dim; i++) {
        slater_det_t *d = det_list + i;
        if ((d->indices[0] == idx_0 && d->indices[1] == idx_1) ||
            (d->indices[0] == idx_1 && d->indices[1] == idx_0)) {
            return i;
        }
    }

    return -1;
}


double complex get_P_ij(int i1, int i2, int j1, int j2, size_t *block_dims, slater_det_t **det_basis, double complex **P_i)
{
    // 1st determinant
    int irrep1 = mulrep2_abelian(get_spinor_irrep(i1), get_spinor_irrep(i2));
    size_t idx_1 = find_determinant(i1, i2, block_dims[irrep1], det_basis[irrep1]);

    // 2nd determinant
    int irrep2 = mulrep2_abelian(get_spinor_irrep(j1), get_spinor_irrep(j2));
    size_t idx_2 = find_determinant(j1, j2, block_dims[irrep2], det_basis[irrep2]);

    if (irrep1 != irrep2) {
        return 0.0 + 0.0*I;
    }
    if (idx_1 == -1 || idx_2 == -1) {
        return 0.0 + 0.0*I;
    }

    //printf("%d %d - %d %d\n", i1, i2, j1, j2);
    //printf("irrep=%d   irrep=%d\n", irrep1, irrep2);
    //printf("dims=%d  idx_1=%d  idx_2=%d\n", block_dims[irrep1], idx_1, idx_2);
    return P_i[irrep1][block_dims[irrep1] * idx_1 + idx_2];
}


void mult_by_shift(char *src_name, char *dst_name, double *spinor_shifts)
{
    copy(src_name, dst_name);

    timer_new_entry("mult_by_shift", "Multiplication by shifts");
    timer_start("mult_by_shift");

    diagram_t *src = diagram_stack_find(src_name);
    if (src == NULL) {
        errquit("mult_by_shift(): diagram '%s' not found", src_name);
    }
    diagram_t *dst = diagram_stack_find(dst_name);

    for (size_t isb = 0; isb < src->n_blocks; isb++) {
        block_t *src_block = src->blocks[isb];
        block_t *dst_block = dst->blocks[isb];
        symblock_load(src_block);
        symblock_load(dst_block);
        if (src_block->is_unique == 0) {
            continue;
        }

        if (src->rank == 4) {
            mult_by_shift_rank_4(src_block, dst_block, spinor_shifts);
        }

        /*else if (rank == 4) {
            diveps_block_rank_4(src_block);
        }
        else if (rank == 6) {
            diveps_block_rank_6(src_block);
        }
        else {
            diveps_block_general(src_block);
        }*/

        symblock_store(dst_block);
    }

    timer_stop("mult_by_shift");
}


void mult_by_shift_rank_4(block_t *src_block, block_t *dst_block, double *spinor_shifts)
{
    int nspinors = get_num_spinors();
    int sect_h = cc_opts->curr_sector_h;
    int sect_p = cc_opts->curr_sector_p;

    int dims_0 = src_block->indices[0][0];
    int dims_1 = src_block->indices[1][0];
    int dims_2 = src_block->indices[2][0];
    int dims_3 = src_block->indices[3][0];
    int coef_0 = dims_1 * dims_2 * dims_3;
    int coef_1 = dims_2 * dims_3;
    int coef_2 = dims_3;
    //int coef_3 = 1;
    int *block_indices_0 = src_block->indices[0] + 1; // +1 to skip first element (length)
    int *block_indices_1 = src_block->indices[1] + 1;
    int *block_indices_2 = src_block->indices[2] + 1;
    int *block_indices_3 = src_block->indices[3] + 1;


    // loop over matrix elements
    for (size_t i = 0; i < src_block->size; i++) {
        // if matrix element is zero we can skip it
        double complex t = carith ? src_block->buf[i] : ((double *) src_block->buf)[i] + 0.0 * I;

        // calculate compound index
        int index = i;
        int idx_0 = index / coef_0;
        index = index % coef_0;
        int idx_1 = index / coef_1;
        index = index % coef_1;
        int idx_2 = index / coef_2;
        index = index % coef_2;
        int idx_3 = index;    // since coef_3 == 1

        // relative index to spinor indices
        idx_0 = block_indices_0[idx_0];
        idx_1 = block_indices_1[idx_1];
        idx_2 = block_indices_2[idx_2];
        idx_3 = block_indices_3[idx_3];

        // calculate denominator and divide matrix element by it
        double shift = spinor_shifts[idx_0] + spinor_shifts[idx_1] - spinor_shifts[idx_2] - spinor_shifts[idx_3];
        t = t * shift;

        if (carith) {
            dst_block->buf[i] = t;
        }
        else {
            ((double *) dst_block->buf)[i] = creal(t);
        }
    }
}


void get_off_diagonal(size_t dim, double complex *src, double complex *off_diag)
{
    memcpy(off_diag, src, sizeof(double complex) * dim * dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            if (i == j) {
                off_diag[dim*i+j] = 0.0+0.0;
            }
        }
    }
}


void matrix_to_tensor(char *tensor_name, size_t *block_dims, slater_det_t **det_basis, double complex **matrix)
{
    diagram_t *tensor = diagram_stack_find(tensor_name);
    if (tensor == NULL) {
        errquit("in matrix_to_tensor(): diagram '%s' not found", tensor_name);
    }

    assert(tensor->rank == 4);

    for (size_t isb = 0; isb < tensor->n_blocks; isb++) {
        block_t *block = tensor->blocks[isb];
        symblock_load(block);
        if (block->is_unique == 0) {
            continue;
        }

        int dims_0 = block->indices[0][0];
        int dims_1 = block->indices[1][0];
        int dims_2 = block->indices[2][0];
        int dims_3 = block->indices[3][0];
        int coef_0 = dims_1 * dims_2 * dims_3;
        int coef_1 = dims_2 * dims_3;
        int coef_2 = dims_3;
        //int coef_3 = 1;
        int *block_indices_0 = block->indices[0] + 1; // +1 to skip first element (length)
        int *block_indices_1 = block->indices[1] + 1;
        int *block_indices_2 = block->indices[2] + 1;
        int *block_indices_3 = block->indices[3] + 1;


        // loop over matrix elements
        for (size_t i = 0; i < block->size; i++) {

            // calculate compound index
            int index = i;
            int idx_0 = index / coef_0;
            index = index % coef_0;
            int idx_1 = index / coef_1;
            index = index % coef_1;
            int idx_2 = index / coef_2;
            index = index % coef_2;
            int idx_3 = index;    // since coef_3 == 1

            // relative index to spinor indices
            idx_0 = block_indices_0[idx_0];
            idx_1 = block_indices_1[idx_1];
            idx_2 = block_indices_2[idx_2];
            idx_3 = block_indices_3[idx_3];

            double complex t = get_P_ij(idx_0, idx_1, idx_2, idx_3, block_dims, det_basis, matrix);


            if (carith) {
                block->buf[i] = t;
            }
            else {
                ((double *) block->buf)[i] = creal(t);
            }
        }

        symblock_store(block);
    }
}
