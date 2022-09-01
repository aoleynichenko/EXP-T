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
 * Restoration of the intermediate normalization of the wave operator
 * Omega = { exp(T) }.
 * This code is used primarily for the FF calculation of off-diagonal matrix
 * elements between the 0h0p ground and 1h1p excited states.
 *
 * For details, see:
 * A. Zaitsevskii, A. V. Oleynichenko, E. Eliav,
 * Finite-Field Calculations of Transition Properties by the Fock Space
 * Relativistic Coupled Cluster Method: Transitions between Different
 * Fock Space Sectors.
 * Symmetry 2020, 12(11), 1845; https://doi.org/10.3390/sym12111845
 *
 * 2020 Alexander Oleynichenko
 */

#include "renorm_omega.h"

#include <assert.h>
#include <string.h>

#include "datamodel.h"
#include "io.h"
#include "mvcoef.h"
#include "slater_det.h"
#include "../linalg/cblas_lapacke.h"


size_t first_nonzero_irrep(size_t *block_dims);
FILE *hefff_open(int sect_h, int sect_p, char *label);
void hefff_close(FILE *hefff);
void hefff_write_block(FILE *hefff, int carith, int rep_no, size_t dim, double complex *heff);
size_t first_nonzero_irrep(size_t *block_dims);
void renormalize_wave_operator_0h0p_0h1p(size_t dim, slater_det_t *det_list,
                                         double complex *heff, double complex *heff_prime, double complex *omega);
int get_nroots_for_irrep(char *irrep_name);
void get_nroots(size_t *block_dims, double complex **eigvalues, size_t *nroots);


void restore_intermediate_normalization(size_t *block_dims, slater_det_t **det_list, double complex **heff,
                                        double complex **eigvalues)
{
    if (!(cc_opts->sector_h == 1 && cc_opts->sector_p == 1)) {
        return;
    }

    int first_irrep = first_nonzero_irrep(block_dims);
    int vacuum_irrep = get_vacuum_irrep();
    char *vacuum_irrep_name = get_irrep_name(vacuum_irrep);
    slater_det_t *det_basis = det_list[vacuum_irrep];
    size_t dim = block_dims[vacuum_irrep];
    double complex *heff_block = heff[vacuum_irrep];

    size_t dim_prime = dim + 1;  // dimension of the Heff' matrix
    double complex *heff_0h0p_1h1p_prime = z_zeros(dim_prime, dim_prime);
    double complex *p_omega_p = z_zeros(dim_prime, dim_prime);

    // 0h0p-1h1p: transformation of the Heff to with the (P\OmegaP)^-1 matrix.
    // + store this 0h0p+1h1p block; it will be written to the formatted file
    // AFTER all blocks belonging to the "pure" 1h1p sector
    renormalize_wave_operator_0h0p_0h1p(dim, det_basis, heff_block, heff_0h0p_1h1p_prime, p_omega_p);

    // write the 0h0p+1h1p transformed Heff block
    // together with the transformation P \Omega P matrix
    FILE *hefff = hefff_open(1, 1, "0011");
    // heff block
    hefff_write_block(hefff, arith == CC_ARITH_COMPLEX, vacuum_irrep-first_irrep+1, dim_prime, heff_0h0p_1h1p_prime);
    // P \omega P
    for (int i = 0; i < dim_prime * dim_prime; i++) {
        if (arith == CC_ARITH_COMPLEX) {
            fprintf(hefff, "%21.12E%21.12E", creal(p_omega_p[i]), cimag(p_omega_p[i]));
            if (i > 0 && i % 2 != 0) { fprintf(hefff, "\n"); }
        }
        else {
            fprintf(hefff, "%21.12E", creal(p_omega_p[i]));
            if (i > 0 && (i + 1) % 4 == 0) { fprintf(hefff, "\n"); }
        }
    }
    if (dim * dim % 2 != 0) { fprintf(hefff, "\n"); }
    hefff_close(hefff);

    /*
     * The 0h0p-1h1p Heff block should be diagonalized since its eigenvectors
     * are required to construct properties and quasi-natural orbitals.
     * These eigenvectors are to be stored in the MVCOEF0011 unformatted file.
     */
    double complex *vl_0011 = z_zeros(dim_prime, dim_prime);
    double complex *vr_0011 = z_zeros(dim_prime, dim_prime);
    double complex *ev_0011 = z_zeros(dim_prime, 1);

    /*int nroots_irep = get_nroots_for_irrep(vacuum_irrep_name);
    size_t nroots = cc_opts->nroots_specified ? nroots_irep : dim_prime;*/
    size_t nroots_array[CC_MAX_NUM_IRREPS];
    get_nroots(block_dims, eigvalues, nroots_array);
    size_t nroots = (cc_opts->nroots_specified || cc_opts->roots_cutoff_specified) ? nroots_array[vacuum_irrep] : dim_prime;

    eigensolver(dim_prime, heff_0h0p_1h1p_prime, ev_0011, vl_0011, vr_0011);
    if (cc_opts->do_hermit == 1) {
        loewdin_orth(dim_prime, vr_0011, vr_0011, cc_opts->print_level == CC_PRINT_DEBUG ? 5 : 0);
        memcpy(vl_0011, vr_0011, dim_prime * dim_prime * sizeof(double complex));
    }

    // print model vectors to the unformatted file MVCOEF**
    slater_det_t *det_list_0011 = (slater_det_t *) cc_malloc(sizeof(slater_det_t) * dim_prime);
    det_list_0011[0].indices[0] = 0;
    det_list_0011[0].indices[1] = 0;
    det_list_0011[0].sym = vacuum_irrep;
    memcpy(det_list_0011 + 1, det_basis, sizeof(slater_det_t) * (dim_prime - 1));
    ev_0011[0] = 0.0 + 0.0 * I;

    int f_mvcoef_0011 = io_open("MVCOEF0011", "w");
    mvcoef_write_vectors_unformatted(f_mvcoef_0011, vacuum_irrep_name, nroots + 1, dim_prime, det_list_0011, ev_0011,
                                     vl_0011, vr_0011);

    //double lowest_root = get_lowest_eigenvalue(block_dims, eigvalues);
    mvcoef_close(f_mvcoef_0011, 0.0);  // lowest root == 0.0
    cc_free(vl_0011);
    cc_free(vr_0011);
    cc_free(ev_0011);

    cc_free(heff_0h0p_1h1p_prime);
    cc_free(p_omega_p);
}


/**
 * constructs matrix of the model-space projection of the wave operator \Omega
 * in the basis of 0h0p (vacuum) and 1h1p (singly excited) determinants:
 * P \Omega P, P=0h0p+1h1p
 * note that the vacuum and 1h1p determinants must be of the same symmetry.
 *
 * T1 = exc cluster oper (0h0p)
 * S1 = de-exc cluster oper (1h1p)
 */
void omega_0h0p_0h1p(size_t dim, slater_det_t *det_list_1h1p, double complex *omega,
                     char *exc_oper_name, char *deexc_oper_name)
{
    diagram_t *dg_exc = diagram_stack_find(exc_oper_name);
    diagram_t *dg_deexc = diagram_stack_find(deexc_oper_name);
    assert(dg_exc != NULL && dg_deexc != NULL);

    /*
     * < 0h0p | Omega | 0h0p > = 1
     */
    omega[0] = 1.0 + 0.0 * I;

    /*
     * < 0h0p | Omega | 1h1p > = S1[a,i]
     */
    for (size_t iket = 0; iket < dim; iket++) {
        slater_det_t *det_ia = det_list_1h1p + iket;
        int i = det_ia->indices[0];
        int a = det_ia->indices[1];
        omega[(dim + 1) * (iket + 1)] = diagram_get_2(dg_deexc, a, i);
    }

    /*
     * < 1h1p | Omega | 0h0p > = T1[i,a]
     */
    for (size_t ibra = 0; ibra < dim; ibra++) {
        slater_det_t *det_ia = det_list_1h1p + ibra;
        int i = det_ia->indices[0];
        int a = det_ia->indices[1];
        omega[ibra + 1] = diagram_get_2(dg_exc, i, a);
    }

    /*
     * < 1h1p i,a | Omega | 1h1p j,b > = \delta_ij \delta_ab + T1[i,a] * S1[j,b]
     */
    for (size_t ibra = 0; ibra < dim; ibra++) {
        for (size_t iket = 0; iket < dim; iket++) {
            slater_det_t *det_ia = det_list_1h1p + ibra;
            slater_det_t *det_jb = det_list_1h1p + iket;
            int i = det_ia->indices[0];
            int a = det_ia->indices[1];
            int j = det_jb->indices[0];
            int b = det_jb->indices[1];
            omega[(dim + 1) * (iket + 1) + ibra + 1] =
                    (i == j) * (a == b) + diagram_get_2(dg_exc, i, a) * diagram_get_2(dg_deexc, b, j);
        }
    }
}


/**
 * Performs renormalization of the wave operator \Omega in order to force
 * intermediate normalization P\OmegaP=P in the P=0h0p+1h1p model space.
 *
 * @param dim dimension of the Heff block in the 1h1p subspace
 * @param det_list list of Slater determinants; their symmetry
 *      must coincide with the symmetry of the vacuum determinant
 * @param heff Heff block in the 1h1p subspace
 * @param heff_prime Heff block in the 0h0p+1h1p subspace, transformed
 *      accordingly with renormalized wave operator \Omega
 *      dim(heff') = dim(heff)+1
 * @param omega transformation matrix (model-space part of the wave operator)
 *      dim(omega) = dim(heff)+1
 */
void renormalize_wave_operator_0h0p_0h1p(size_t dim, slater_det_t *det_list,
                                         double complex *heff, double complex *heff_prime, double complex *omega)
{
    size_t dimx = dim + 1;  // dimension of eXtended (0h0p+1h1p) matrices
    double complex alpha = 1.0 + 0.0 * I;
    double complex beta = 0.0 + 0.0 * I;

    /*
     * allocate working arrays, init with zeros
     */
    memset(omega, 0, sizeof(double complex) * dimx * dimx);
    double complex *omega_inv = z_zeros(dimx, dimx);
    double complex *heff_0h0p_1h1p = z_zeros(dimx, dimx);
    double complex *buf = z_zeros(dimx, dimx);

    /*
     * construct 0h0p+1h1p (extended) block-diagonal Heff matrix
     * Heff 1h1p block will be placed into the right lower angle.
     * Ecorr is subtracted from the whole diagonal.
     */
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            heff_0h0p_1h1p[dimx * (i + 1) + j + 1] = heff[dim * i + j];
        }
    }

    /*
     * construct P \Omega P matrix (P=0h0p+1h1p)
     */
    omega_0h0p_0h1p(dim, det_list, omega, "t1c", "e1c");

    /*
     * construct inverse matrix (P \Omega P)^{-1}
     */
    inverse_matrix(dimx, omega, omega_inv);

    /*
     * Heff' = { (P\OmegaP) * { Heff * (P\OmegaP)^-1 } }
     */
    xgemm(CC_DOUBLE_COMPLEX, "N", "N",
                dimx, dimx, dimx, &alpha, heff_0h0p_1h1p, dimx, omega_inv, dimx, &beta, buf, dimx);
    xgemm(CC_DOUBLE_COMPLEX, "N", "N",
                dimx, dimx, dimx, &alpha, omega, dimx, buf, dimx, &beta, heff_prime, dimx);

    /*
     * deallocate working arrays
     */
    cc_free(omega_inv);
    cc_free(heff_0h0p_1h1p);
    cc_free(buf);
}
