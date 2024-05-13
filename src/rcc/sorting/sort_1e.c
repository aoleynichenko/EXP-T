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

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sorting_request.h"

#include "error.h"
#include "io.h"
#include "memory.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"

void fill_block_one_elec(block_t *block, double complex *ints_matrix, int ignore_diagonal);

void reconstruct_fock(int nspinors, double complex *h_matrix, double complex *fock_matrix);

double recalculate_scf_energy(double complex *h_ints);

int read_prop_two_files(int nspinors, char *file_re, char *file_im, double complex *prop_mat);

int read_prop_single_file(int nspinors, char *prop_name, double complex *prop_mat);

void fock_add_oneprop(double complex *fock, double complex lambda, char *file_re, char *file_im);

void fock_add_mdprop(double complex *fock, double complex lambda, char *prop_name);

double max_diagonal_diff(double complex *fock);


void sort_onel()
{
    int nspinors = get_num_spinors();
    double complex *h_ints = (double complex *) x_zeros(CC_COMPLEX, nspinors, nspinors);
    double complex *f_ints = (double complex *) x_zeros(CC_COMPLEX, nspinors, nspinors);

    printf(" sorting one-electron integrals ...\n");

    // read one-electron integrals -- core Fock operator
    // ? а что если HINT нет?
    int fd = io_open("HINT", "r");
    io_read_compressed(fd, h_ints, sizeof(double complex) * nspinors * nspinors);
    io_close(fd);

    // read one-electron operators from the OneProp code by Leonid V. Skripnikov
    if (cc_opts->oneprop_on) {
        for (int ioper = 0; ioper < cc_opts->n_oneprop; ioper++) {
            fock_add_oneprop(h_ints, cc_opts->oneprop_lambda[ioper],
                             cc_opts->oneprop_file_re[ioper], cc_opts->oneprop_file_im[ioper]);
        }
    }

    if (cc_opts->n_mdprop > 0) {
        for (int ioper = 0; ioper < cc_opts->n_mdprop; ioper++) {
            fock_add_mdprop(h_ints, cc_opts->mdprop_lambda[ioper], cc_opts->mdprop_file[ioper]);
        }
    }

    // construct Fock matrix
    if (!cc_opts->x2cmmf) {
        reconstruct_fock(nspinors, h_ints, f_ints);
    }
    else {  // x2cmmf hamiltonian, no recompute fock matrix
        printf("   Fock matrix reconstruction will be skipped\n");
        memset(f_ints, 0, sizeof(double complex) * nspinors * nspinors);
        for (int i = 0; i < nspinors; i++) {
            f_ints[i * nspinors + i] = spinor_info[i].eps;
        }
    }

    // recalculate SCF energy (and update it if needed)
    if (!cc_opts->x2cmmf) {
        double new_escf = recalculate_scf_energy(h_ints);
        printf("   SCF energy (energy of reference determinant) = %20.12f a.u.\n", new_escf);
        if (fabs(cc_opts->escf - new_escf) > 1e-10) {
            printf("   SCF energy (energy of reference determinant) was updated:\n");
            printf("     old energy = %20.12f a.u.\n", cc_opts->escf);
            printf("     new energy = %20.12f a.u.\n", new_escf);
            cc_opts->escf = new_escf;
        }
    }

    // compare diagonal elements of the reconstructed Fock matrix with orbital
    // energies read from the MRCONEE file. if they don't coincide then
    // substitute eps[] with new Fock matrix diagonal elements for the subsequent
    // use in perturbation expressions
    double max_eps_diff = max_diagonal_diff(f_ints);
    if ((max_eps_diff > 1e-10) && (cc_opts->use_oe == 0)) {
        if (max_eps_diff > 1e-6) {
            printf("   The diagonal elements of the reconstructed Fock matrix don't "
                   "coincide with the orbital energies!\n");
            printf("   NOTE: The diagonal elements of the recomputed Fock matrix "
                   "(right column) are used in perturbation expressions.\n");

            printf("\n");
            printf("   no    rep         occ    active     one-el energy        recalc energy            delta    \n");
            printf("  ---------------------------------------------------------------------------------------------\n");
            for (int i = 0; i < nspinors; i++) {
                double eps = spinor_info[i].eps;
                double f_ii = creal(f_ints[i * nspinors + i]);
                printf("  %4d%4d \"%-6s\"%4d       %1s     %16.10f     %16.10f     %16.6e\n",
                       i + 1, spinor_info[i].repno,
                       get_irrep_name(spinor_info[i].repno), is_hole(i),
                       is_active(i) ? "a" : "i", eps, f_ii, f_ii - eps);
            }
            printf("  ---------------------------------------------------------------------------------------------\n");
        }
        printf(" max deviation of the diagonal elements of the reconstructed"
               " Fock matrix and orbital energies = %.6e\n", max_eps_diff);

        // recalculate energies
        for (int i = 0; i < nspinors; i++) {
            spinor_info[i].eps = creal(f_ints[i * nspinors + i]);
            f_ints[i * nspinors + i] = 0.0 + 0.0 * I;
        }
    }

    // только запросы ?
    printf(" fill 1-electron diagrams ... ");
    for (int ireq = 0; ireq < n_requests; ireq++) {
        diagram_t *dg = sorting_requests[ireq].dg;
        if (dg->rank != 2) { // only one-electron diagrams
            continue;
        }
        printf("%s ", dg->name);

        // assign values of matrix elements
        for (size_t isb = 0; isb < dg->n_blocks; isb++) {
            block_t *block = dg->blocks[isb];
            block_load(block);
            int ignore_diagonal = cc_opts->use_oe ? 0 : 1;
            fill_block_one_elec(block, f_ints, ignore_diagonal);
            block_unload(block);
        }
    }

    cc_free(h_ints);
    cc_free(f_ints);

    printf("done\n");
}


void sort_prop(int nspinors, double complex *oper_ints)
{
    printf(" sorting property integrals ...\n");

    for (int ireq = 0; ireq < n_requests; ireq++) {
        diagram_t *dg = sorting_requests[ireq].dg;
        if (dg->rank != 2) { // only one-electron diagrams
            continue;
        }
        // assign values of matrix elements
        for (size_t isb = 0; isb < dg->n_blocks; isb++) {
            block_t *block = dg->blocks[isb];
            block_load(block);
            fill_block_one_elec(block, oper_ints, 0);
            block_unload(block);
        }
    }

    // finalize:
    n_requests = 0;
}


void fock_add_oneprop(double complex *fock, double complex lambda, char *file_re, char *file_im)
{
    int nspinors = get_num_spinors();
    double complex *oper = x_zeros(CC_COMPLEX, nspinors, nspinors);

    printf("   reading one-electron property from Oneprop "
           "(lambda = %e %e, file_re = %s, file_im = %s)\n",
           creal(lambda), cimag(lambda), file_re, file_im);

    // read files with (a) real part of the prop matrix (b) imaginary part
    int status = read_prop_two_files(nspinors, file_re, file_im, oper);
    if (status == EXIT_FAILURE) {
        errquit("Cannot open OneProp files '%s' (re) '%s' (im)", file_re, file_im);
    }

    // add operator matrix to the Fock matrix:
    // F = F + lambda*Oper (lambda = perturbation parameter)
    xaxpy(CC_COMPLEX, nspinors * nspinors, lambda, oper, fock);

    // all done
    cc_free(oper);
}


void fock_add_mdprop(double complex *fock, double complex lambda, char *prop_name)
{
    int nspinors = get_num_spinors();
    double complex *oper = x_zeros(CC_COMPLEX, nspinors, nspinors);

    printf("   reading one-electron property %s from MDPROP "
           "(lambda = %e %e)\n", prop_name, creal(lambda), cimag(lambda));

    // read file containing the prop matrix
    int status = read_prop_single_file(nspinors, prop_name, oper);
    if (status == EXIT_FAILURE) {
        errquit("Cannot open property matrix file '%s'", prop_name);
    }

    // add operator matrix to the Fock matrix:
    // F = F + lambda*Oper (lambda = perturbation parameter)
    xaxpy(CC_COMPLEX, nspinors * nspinors, lambda, oper, fock);

    // all done
    cc_free(oper);
}


/**
 * compares diagonal elements of the reconstructed Fock matrix with orbital
 * energies read from the MRCONEE file;
 * returns maximum difference between them.
 */
double max_diagonal_diff(double complex *fock)
{
    double eps_diff = 0.0;
    double max_eps_diff = 0.0;
    double const PRINT_THRESH = 1e-13;
    int nspinors = get_num_spinors();

    for (int i = 0; i < nspinors; i++) {
        double complex f_ii = fock[i * nspinors + i];
        if (cimag(f_ii) > PRINT_THRESH) {
            printf("   (!) in sort_onel(): imaginary value of the [%d,%d] diagonal "
                   "element of the Fock matrix = %.8f %.3E (threshold = 1e-13)\n",
                   i, i, creal(f_ii), cimag(f_ii));
        }
        eps_diff = cabs(f_ii - spinor_info[i].eps);
        if (eps_diff > max_eps_diff) {
            max_eps_diff = eps_diff;
        }
    }

    return max_eps_diff;
}


/**
 * copy matrix elements from square matrix into the block.
 */
void fill_block_one_elec(block_t *block, double complex *ints_matrix, int ignore_diagonal)
{
    assert(block->rank == 2);

    int nspinors = get_num_spinors();

    int dims_1 = block->shape[0];
    int dims_2 = block->shape[1];
    int *block_indices_1 = block->indices[0];
    int *block_indices_2 = block->indices[1];

    double *dbuf = (double *) block->buf;

    size_t index = 0;
    for (int i1 = 0; i1 < dims_1; i1++) {
        for (int i2 = 0; i2 < dims_2; i2++) {
            int idx_1 = block_indices_1[i1];
            int idx_2 = block_indices_2[i2];

            if ((idx_1 == idx_2) && ignore_diagonal) {
                if (arith == CC_ARITH_COMPLEX) {
                    block->buf[index] = 0.0 + 0.0 * I;
                }
                else {
                    dbuf[index] = 0.0;
                }
            }
            else {
                if (arith == CC_ARITH_COMPLEX) {
                    block->buf[index] = ints_matrix[idx_1 * nspinors + idx_2];
                }
                else {
                    dbuf[index] = creal(ints_matrix[idx_1 * nspinors + idx_2]);
                }
            }
            index++;
        }
    }
}


void reconstruct_fock(int nspinors, double complex *h_matrix, double complex *fock_matrix)
{
    int idx4[4];

    // get pointer to the diagrams required for the construction of the
    // exchange-correlation part of the Fock matrix
    diagram_t *dg_hhhh = diagram_stack_find("hhhh");
    if (dg_hhhh == NULL) {
        errquit("in sort_onel(): diagram 'hhhh' is required to construct Fock matrix, but it was not found");
    }
    diagram_t *dg_hphh = diagram_stack_find("hphh");
    if (dg_hphh == NULL) {
        errquit("in sort_onel(): diagram 'hphh' is required to construct Fock matrix, but it was not found");
    }
    diagram_t *dg_hhhp = diagram_stack_find("hhhp");
    if (dg_hhhp == NULL) {
        errquit("in sort_onel(): diagram 'hhhp' is required to construct Fock matrix, but it was not found");
    }
    diagram_t *dg_hphp = diagram_stack_find("hphp");
    if (dg_hphp == NULL) {
        errquit("in sort_onel(): diagram 'hphp' is required to construct Fock matrix, but it was not found");
    }

    int nocc = get_num_electrons();
    int *occ_spinors = (int *) cc_malloc(sizeof(int) * nocc);
    get_spinor_indices_occupied(occ_spinors);

    printf("   Fock matrix reconstruction ...\n");
    for (int i = 0; i < nspinors; i++) {
        for (int j = 0; j < nspinors; j++) {
            fock_matrix[i * nspinors + j] = h_matrix[i * nspinors + j];

            for (int h = 0; h < nocc; h++) {
                int idx_h = occ_spinors[h];
                idx4[0] = idx_h;
                idx4[1] = i;
                idx4[2] = idx_h;
                idx4[3] = j;
                if (is_hole(i)) {
                    if (is_hole(j)) {
                        fock_matrix[i * nspinors + j] += diagram_get(dg_hhhh, idx4);
                    }
                    else { // 'j is particle'
                        fock_matrix[i * nspinors + j] += diagram_get(dg_hhhp, idx4);
                    }
                }
                else { // 'i' is particle
                    if (is_hole(j)) {
                        fock_matrix[i * nspinors + j] += diagram_get(dg_hphh, idx4);
                    }
                    else { // 'j is particle'
                        fock_matrix[i * nspinors + j] += diagram_get(dg_hphp, idx4);
                    }
                }
            }
        }
    }

    cc_free(occ_spinors);
}


double recalculate_scf_energy(double complex *h_ints)
{
    int idx4[4];
    double new_escf = cc_opts->enuc;
    int nspinors = get_num_spinors();

    diagram_t *dg_hhhh = diagram_stack_find("hhhh");
    if (dg_hhhh == NULL) {
        errquit("in sort_onel(): diagram 'hhhh' is required to construct Fock matrix, but it was not found");
    }

    for (int i = 0; i < nspinors; i++) {
        if (!is_hole(i)) { continue; }
        new_escf += h_ints[i * nspinors + i];
        for (int j = 0; j < nspinors; j++) {
            if (!is_hole(j)) { continue; }
            idx4[0] = i;
            idx4[1] = j;
            idx4[2] = i;
            idx4[3] = j;
            new_escf += 0.5 * diagram_get(dg_hhhh, idx4);
        }
    }

    return new_escf;
}


/**
 * Reads pair of formatted files containing real square matrices of properties.
 * NOTE: transposition is included to be able to use files from DIRAC.
 *
 * @param nspinors number of spinors
 * @param file_re real part
 * @param file_im imag part
 * @param prop_mat (is returned) NSPINORS x NSPINORS complex matrix
 */
int read_prop_two_files(int nspinors, char *file_re, char *file_im, double complex *prop_mat)
{
    FILE *f_prop;
    int idx_i, idx_j;
    double value;

    memset(prop_mat, 0, nspinors * nspinors * sizeof(double complex));

    // read files with (a) real part of the prop matrix
    f_prop = fopen(file_re, "r");
    if (f_prop == NULL) {
        return EXIT_FAILURE;
    }
    while (fscanf(f_prop, "%d%d%lf", &idx_i, &idx_j, &value) == 3) {
        prop_mat[(idx_i - 1) * nspinors + (idx_j - 1)] = value + 0.0 * I;
    }
    fclose(f_prop);

    // (b) imaginary part
    f_prop = fopen(file_im, "r");
    if (f_prop == NULL) {
        return EXIT_FAILURE;
    }
    while (fscanf(f_prop, "%d%d%lf", &idx_i, &idx_j, &value) == 3) {
        prop_mat[(idx_i - 1) * nspinors + (idx_j - 1)] += 0.0 + value * I;
    }
    fclose(f_prop);

    return EXIT_SUCCESS;
}


/**
 * Reads property matrix from file 'prop_name' (extracted from MDPROP)
 *
 * @param nspinors number of spinors
 * @param prop_name name of the property under consideration in the MDPROP file
 * @param prop_mat (is returned) NSPINORS x NSPINORS complex matrix
 */
int read_prop_single_file(int nspinors, char *prop_name, double complex *prop_mat)
{
    FILE *f_prop;
    int idx_i, idx_j;
    double re_value, im_value;

    memset(prop_mat, 0, nspinors * nspinors * sizeof(double complex));

    // read file with the prop matrix
    f_prop = fopen(prop_name, "r");
    if (f_prop == NULL) {
        return EXIT_FAILURE;
    }
    while (fscanf(f_prop, "%d%d%lf%lf", &idx_i, &idx_j, &re_value, &im_value) == 4) {
        prop_mat[(idx_i - 1) * nspinors + (idx_j - 1)] = re_value + im_value * I;
    }
    fclose(f_prop);

    return EXIT_SUCCESS;
}
