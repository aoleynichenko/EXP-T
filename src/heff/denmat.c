/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2020 The EXP-T developers.
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
 * denmat.c
 *
 * Calculation of density matrices and natural orbitals.
 * Associated header file: heff.h
 *
 * 2019 Alexander Oleynichenko
 ******************************************************************************/

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "io.h"
#include "linalg.h"
#include "memory.h"
#include "mvcoef.h"
#include "slater.h"
#include "spinors.h"
#include "symmetry.h"
#include "options.h"
#include "utils.h"

// functions used in this module
int cmp_natorb_occupations(double complex occ1, double complex occ2);

void write_NO(char *natorb_file_name, int sect_h, int sect_p,
              int dim, double *occ_numbers, double complex *nat_orbs,
              double occ_thresh);

int density_matrix_element(int sect_h, int sect_p, int p, int q, slater_det_t *bra, slater_det_t *ket);

void sort_vectors(int n, double complex *ev, double complex *vl, double complex *vr,
                  int (*cmp)(double complex, double complex));


/*******************************************************************************
 * density_matrix
 *
 * Computes density matrix for (a) selected electronic state; (b) selected pair
 * of electronic states (transition density matrix).
 * Performs diagonalization (in order to obtain NOs or NTOs).
 *
 * NOTE: density matrices computed are not true density matrices since only
 * model vectors are used to compute them. However, this do not seem to be a
 * problem, since density matrices, NOs and NTOs produced are used only for
 * semi-quantitative estimation and visualization purposes.
 *
 * Algorithm notes:
 * (1) Only 'valence' part of the density matrix is constructed =>
 *     dim = nacth * (h > 0) + nactp * (p > 0) for the (h,p) FS sector
 * (2) in all formulas for density matrices we must use global indices of
 *     spinors, however, rows and columns of DMs are enumerated with local
 *     indices; local indices are transformed into global when needed.
 * (3) the 1h1p sector is the special case:
 *     - calculation of approximate natural orbitals for the 0h0p state is
 *       meaningless (NOs = MOs)
 *     - 0h0p - 1h1p transitions ("CIS-like") require separate code for the DM
 *       evaluation (again, two cases: 0h0p->1h1p and 1h1p->0h0p, are not equal
 *       in case of non-hermitian Hamiltonian
 *     - 1h1p - 1h1p transitions can be calculated in the usual way
 *     NOTE: state numbering for the case of the 1h1p sector as a target
 *       must be quite different: state 1 == 0h0p, states from 1h1p: 2,3,4...
 *       (numbering is shifted by 1)
 *
 ******************************************************************************/
void density_matrix(int sect_h, int sect_p, int rep1, int state1, int rep2, int state2)
{
    int nrep;
    struct mv_block mv_blocks[64];
    double const occ_thresh = 1e-6;   // threshold for printing natural orbitals
    double const coef_thresh = 1e-4;   // threshold for printing model vec-s coeff-s
    // for state 1
    char *rep_name1;
    size_t ms_size1;
    slater_det_t *dets1;
    double eigval_1;
    double energy_cm_1;
    double complex *coef_left1, *coef_right1;
    // for state 2
    char *rep_name2;
    size_t ms_size2;
    slater_det_t *dets2;
    double eigval_2;
    double energy_cm_2;
    double complex *coef_left2, *coef_right2;

    rep_name1 = get_irrep_name(rep1);
    rep_name2 = get_irrep_name(rep2);

    // NOs or NTOs (transitions) ?
    int tran = !((rep1 == rep2) && (state1 == state2));

    // print header for this module
    if (!tran) {
        printf("\n\n *** DENSITY MATRIX AND NATURAL ORBITALS FOR THE IRREP %d (%s) STATE %d ***\n", rep1 + 1, rep_name1,
               state1 + 1);
    }
    else{
        printf("\n\n ***   TRANSITION DENSITY MATRIX AND NATURAL TRANSITION ORBITALS   ***\n");
        printf(" *** FOR THE IRREP %2d (%s) STATE %2d  --->  IRREP %2d (%s) STATE %2d TRANSITION ***\n",
               rep1, rep_name1, state1 + 1, rep2, rep_name2, state2 + 1);
    }

    // only 'valence' parts of DMs are be constructed
    // dim DMs = n_active x n_active
    int n_active = 0;
    int active_spinors[CC_MAX_SPINORS]; // local -> global spinor index mapping
    double complex *denmat;
    double complex *nat_occ, *natorb_left, *natorb_right;
    double *lambda;

    // TODO: refactor, separate utility function (move to spinors.c)
    for (int i = 0; i < nspinors; i++){
        if (is_act_hole(i) && sect_h > 0 || is_act_part(i) && sect_p > 0) {
            active_spinors[n_active] = i;
            n_active++;
        }
    }

    // in the 1h1p sector states:
    // 0h0p => state 1
    // 1h1p => states 2, 3, ...
    // (only in the irrep containing the 0h0p state)
    if (sect_h == 1 && sect_p == 1) {
        int vac_irrep = get_vacuum_irrep();
        if (rep1 == vac_irrep) {
            state1 -= 1;
        }
        if (rep2 == vac_irrep) {
            state2 -= 1;
        }
        if (state1 == -1 && state2 == -1) {
            printf(" Calculation of approximate natural orbitals for the vacuum state (0h0p)\n");
            printf(" is meaningless and will be ignored\n");
        }
    }

    // allocate working arrays
    denmat  = zzeros(n_active, n_active);
    nat_occ = zzeros(n_active, 1);  // NO occ numbers
    lambda  = dzeros(n_active, 1);     // NTO singular values
    natorb_left  = zzeros(n_active, n_active); // for NTOs too
    natorb_right = zzeros(n_active, n_active);

    // extract model vectors and eigenvalues from the MVCOEF* unformatted file
    read_model_vectors_unformatted(sect_h, sect_p, &nrep, mv_blocks);
    struct mv_block *mvb1 = NULL, *mvb2 = NULL;
    for (size_t ib = 0; ib < nrep; ib++){
        if (strcmp(mv_blocks[ib].rep_name, rep_name1) == 0) {
            mvb1 = mv_blocks + ib;
        }
        if (strcmp(mv_blocks[ib].rep_name, rep_name2) == 0) {
            mvb2 = mv_blocks + ib;
        }
    }

    // extract data for state 1
    if (state1 != -1) {
        ms_size1 = mvb1->ms_size;
        dets1 = mvb1->dets;
        eigval_1 = creal(mvb1->eigval[state1]);
        energy_cm_1 = mvb1->energy_cm[state1];
        coef_left1 = mvb1->vl + ms_size1 * state1;
        coef_right1 = mvb1->vr + ms_size1 * state1;
    }
    else{
        ms_size1 = 1;
        dets1 = (slater_det_t *) cc_malloc(sizeof(slater_det_t) * 1);
        set_vacuum_det(dets1);
        eigval_1 = 0.0 + 0.0 * I;
        energy_cm_1 = 0.0;
        coef_left1 = zzeros(1, 1);
        coef_left1[0] = 1.0 + 0.0 * I;
    }
    // extract data for state 2 (if tran == 0 => the same as for state 1)
    if (state2 != -1) {
        ms_size2 = mvb2->ms_size;
        dets2 = mvb2->dets;
        eigval_2 = creal(mvb2->eigval[state2]);
        energy_cm_2 = mvb2->energy_cm[state2];
        coef_left2 = mvb2->vl + ms_size2 * state2;
        coef_right2 = mvb2->vr + ms_size2 * state2;
    }
    else{
        ms_size2 = 1;
        dets1 = (slater_det_t *) cc_malloc(sizeof(slater_det_t) * 1);
        set_vacuum_det(dets1);
        eigval_2 = 0.0 + 0.0 * I;
        energy_cm_2 = 0.0;
        coef_right2 = zzeros(1, 1);
        coef_right2[0] = 1.0 + 0.0 * I;
    }

    // density matrix (DM) construction
    // the same code for both cases of DMs and transition DMs
    // TODO: separate routine ?
    for (int p = 0; p < n_active; p++){
        for (int q = 0; q < n_active; q++){
            double complex d_pq = 0.0 + 0.0 * I;
            // loops over model vectors
            for (size_t i = 0; i < ms_size1; i++){
                for (size_t j = 0; j < ms_size2; j++){
                    slater_det_t *bra = &dets1[i];
                    slater_det_t *ket = &dets2[j];
                    int bra_d_ket = density_matrix_element(sect_h, sect_p, active_spinors[p], active_spinors[q], bra,
                                                           ket);
                    d_pq += conj(coef_left1[i]) * coef_right2[j] * bra_d_ket;
                }
            }
            denmat[p * n_active + q] = d_pq;
        }
    }

    // just in order to check correctness: calculate transition dipole moment
    double complex tdm[] = {0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I};
    char *file_names[] = {"XDIPLEN", "YDIPLEN", "ZDIPLEN"};
    for (int icoord = 0; icoord < 3; icoord++){
        double complex *d_spinor = zzeros(nspinors, nspinors);
        FILE *f = fopen(file_names[icoord], "r");
        if (f == NULL) {
            printf("Cannot calculate transition dipole moment!\n");
            printf("File '%s' does not exist\n", file_names[icoord]);
            break;
        }
        int idx1, idx2;
        double re, im;
        while (fscanf(f, "%d%d%lf%lf", &idx1, &idx2, &re, &im) == 4) {
            d_spinor[(idx2 - 1) * nspinors + (idx1 - 1)] = re + I * im;  // together with transposition
        }
        fclose(f);
        tdm[icoord] = 0.0 + 0.0 * I;
        for (int i = 0; i < n_active; i++){
            for (int j = 0; j < n_active; j++){
                idx1 = active_spinors[i];
                idx2 = active_spinors[j];
                tdm[icoord] += denmat[i * n_active + j] * d_spinor[idx1 * nspinors + idx2];
            }
        }
        cc_free(d_spinor);
    }
    double complex dx = tdm[0];
    double complex dy = tdm[1];
    double complex dz = tdm[2];
    double tranmom = sqrt(cabs(dx) * cabs(dx) + cabs(dy) * cabs(dy) + cabs(dz) * cabs(dz));


    // calculate natural transition spinors (NTOs)
    // SVD of transition density matrix
    if (tran) {
        // SVD: T = U*S*VH
        // U = eigenvectors of T*TH
        // V = eigenvectors of TH*T
        // S = sqrt(eigenvalues of U or V)
        // pairs (eigenvalue,eigenvector) must be sorted in descending order (by eigenvalue)
        svd(n_active, denmat, lambda, natorb_left, natorb_right);
    }
        // calculate natural spinors (NOs)
        // diagonalize density matrix
    else{
        // тут происходит очень неустойчивая диагонализация матрицы с кучей нулевых собственных значений
        // но почему все нормально на неортогонализованных векторах?
        // потому что диагонализуемая матрица плотности неэрмитова!
        eig(n_active, denmat, nat_occ, natorb_left, natorb_right);
        sort_vectors(n_active, nat_occ, natorb_left, natorb_right, cmp_natorb_occupations);
    }

    // print occ-s and NOs to stdout
    if (tran) {
        printf(" Symmetry: %s --> %s\n", rep_name1, rep_name2);
        printf(" State 1: %10.2f cm^-1 (%s)\n", energy_cm_1, rep_name1);
        printf(" State 2: %10.2f cm^-1 (%s)\n", energy_cm_2, rep_name2);
        printf(" (dipole-length) transition dipole moment (a.u.):\n");
        printf("  dx = %10.6f%10.6f     |dx| = %10.6f\n", creal(dx), cimag(dx), cabs(dx));
        printf("  dy = %10.6f%10.6f     |dy| = %10.6f\n", creal(dy), cimag(dy), cabs(dy));
        printf("  dz = %10.6f%10.6f     |dz| = %10.6f\n", creal(dz), cimag(dz), cabs(dz));
        printf("                                |d|  = %10.6f\n", tranmom);
        double sum_lambda = 0.0;  // trace
        for (int i = 0; i < n_active; i++){
            sum_lambda += lambda[i];
            if (lambda[i] < occ_thresh) {
                continue;
            }

            printf(" Singular value = %.6f\n", lambda[i]);
            printf("    from:\n");
            for (size_t j = 0; j < n_active; j++){
                double complex coef = natorb_left[n_active * i + j];
                if (cabs(coef) < coef_thresh) {
                    continue;
                }
                printf("  %12.6f%12.6f ", creal(coef), cimag(coef));
                int ispinor = active_spinors[j];
                printf("  %4s #%4d (%12.6f)\n", rep_names[spinor_info[ispinor].repno], ispinor + 1,
                       spinor_info[ispinor].eps);
            }
            printf("    to:\n");
            for (size_t j = 0; j < n_active; j++){
                double complex coef = natorb_right[n_active * i + j];
                if (cabs(coef) < coef_thresh) {
                    continue;
                }
                printf("  %12.6f%12.6f ", creal(coef), cimag(coef));
                int ispinor = active_spinors[j];
                printf("  %4s #%4d (%12.6f)\n", rep_names[spinor_info[ispinor].repno], ispinor + 1,
                       spinor_info[ispinor].eps);
            }
        }
        printf(" Sum of singular values lambda = %.6f\n", sum_lambda);

        // flush full information about NOs to the formatted file
        char natorb_file_name[256];
        sprintf(natorb_file_name, "NTO_%dh%dp_%s:%d-%s:%d_FROM.dat", sect_h, sect_p, rep_name1, state1 + 1, rep_name2,
                state2 + 1);
        str_replace(natorb_file_name, '/', '|');  // for some irrep names with '/' in order to prevent bugs
        write_NO(natorb_file_name, sect_h, sect_p, n_active, lambda, natorb_left, occ_thresh);
        sprintf(natorb_file_name, "NTO_%dh%dp_%s:%d-%s:%d_TO.dat", sect_h, sect_p, rep_name1, state1 + 1, rep_name2,
                state2 + 1);
        str_replace(natorb_file_name, '/', '|');
        write_NO(natorb_file_name, sect_h, sect_p, n_active, lambda, natorb_right, occ_thresh);
    }
    else{
        printf(" Symmetry of this electronic state: %s\n", rep_name1);
        printf(" Eigenvalue   = %.8f a.u.\n", eigval_1);
        printf(" Energy level = %.2f cm^-1\n", energy_cm_1);
        double complex sum_occ = 0.0 + 0.0 * I;  // trace
        for (int i = 0; i < n_active; i++){
            sum_occ += nat_occ[i];
            if (cimag(nat_occ[i]) > 1e-13) {
                printf(" Imaginary occupation number: %.3e %.3e\n", creal(nat_occ[i]), cimag(nat_occ[i]));
            }
            if (fabs(creal(nat_occ[i])) < occ_thresh) {
                continue;
            }

            printf(" Occ number = %.6f\n", creal(nat_occ[i]));
            for (size_t j = 0; j < n_active; j++){
                double complex coef = natorb_right[n_active * i + j];
                if (cabs(coef) < coef_thresh) {
                    continue;
                }
                printf("  %12.6f%12.6f ", creal(coef), cimag(coef));
                int ispinor = active_spinors[j];
                printf("  %4s #%4d (%12.6f)\n", rep_names[spinor_info[ispinor].repno], ispinor + 1,
                       spinor_info[ispinor].eps);
            }
        }
        printf(" Sum of occ numbers = %.6f (must be = %d)\n", creal(sum_occ), -sect_h + sect_p);

        // flush full information about NOs to the formatted file
        char natorb_file_name[256];
        sprintf(natorb_file_name, "NATORB_%dh%dp_%s:%d.dat", sect_h, sect_p, rep_name1, state1 + 1);
        str_replace(natorb_file_name, '/', '|');  // for some irrep names with '/' in order to prevent bugs
        for (int i = 0; i < n_active; i++){
            lambda[i] = creal(nat_occ[i]);
        }
        write_NO(natorb_file_name, sect_h, sect_p, n_active, lambda, natorb_right, occ_thresh);
    }


    // cleanup
    cc_free(denmat);
    cc_free(nat_occ);
    cc_free(lambda);
    cc_free(natorb_left);
    cc_free(natorb_right);
    for (size_t irep = 0; irep < nrep; irep++){
        struct mv_block *b = &mv_blocks[irep];
        cc_free(b->dets);
        cc_free(b->eigval);
        cc_free(b->energy_cm);
        cc_free(b->vl);
        cc_free(b->vr);
    }

    if (tran) {
        printf(" *** end of TRANSITION DENSITY MATRIX AND NATURAL TRANSITION ORBITALS module ***\n");
    }
    else{
        printf(" *** end of DENSITY MATRIX AND NATURAL ORBITALS module ***\n");
    }
}


/*******************************************************************************
 * write_NO
 *
 * Write natural orbital (spinor) expansions to the formatted file with
 * the standard name NATORB_[H]h[P]p_[REP]:[STATE].dat
 ******************************************************************************/
void write_NO(char *natorb_file_name, int sect_h, int sect_p,
              int dim, double *occ_numbers, double complex *nat_orbs,
              double occ_thresh)
{
    FILE *f;
    int n_natorb;
    int n_active;
    int active_spinors[CC_MAX_SPINORS]; // local -> global spinor index mapping

    // construct 'local indices -> global indices' mapping
    n_active = 0;
    for (int i = 0; i < nspinors; i++){
        if (is_act_hole(i) && sect_h > 0 || is_act_part(i) && sect_p > 0) {
            active_spinors[n_active] = i;
            n_active++;
        }
    }

    // header contains brief info about spinors with numbering
    f = fopen(natorb_file_name, "w");
    fprintf(f, "dim %d\n", n_active);
    fprintf(f, "spinor info:\n");
    for (int i = 0; i < n_active; i++){
        int ispinor = active_spinors[i];
        fprintf(f, "%4d%20.12f%10s\n", ispinor + 1, spinor_info[ispinor].eps, rep_names[spinor_info[ispinor].repno]);
    }

    // calculate number of NOs to be printed (with occ_no >= thresh)
    n_natorb = 0;
    for (int i = 0; i < n_active; i++){
        if (fabs(creal(occ_numbers[i])) >= occ_thresh) {
            n_natorb++;
        }
    }
    fprintf(f, "nspinors %d\n", n_natorb);

    // write NOs
    for (int i = 0; i < n_active; i++){
        double occ = occ_numbers[i];
        if (fabs(occ) < occ_thresh) {
            continue;
        }

        fprintf(f, "occ %.6f\n", occ);

        // write coefficients (all)
        for (size_t j = 0; j < n_active; j++){
            double complex coef = nat_orbs[n_active * i + j];
            int ispinor = active_spinors[j];
            fprintf(f, "%4d  %20.12e%20.12e\n", ispinor + 1, creal(coef), cimag(coef));
        }
    }
    printf(" written to file: %s\n", natorb_file_name);

    fclose(f);
}


/*******************************************************************************
 * cmp_natorb_occupations
 *
 * Comparator function for rearrangement of natural orbitals by their eigenvalues
 * -- occupation numbers
 *
 * Rules of comparation:
 * (1) negative, then positive
 * (2) negative occ numbers in ascending order: -1 --> 0
 * (3) positive occ numbers in descending order: +1 --> 0
 ******************************************************************************/
int cmp_natorb_occupations(double complex occ1, double complex occ2)
{
    double occ1_re = creal(occ1);
    double occ2_re = creal(occ2);

    if (occ1_re * occ2_re < 0) {
        return sgn(occ1_re - occ2_re);  // ascending
    }
    else{
        return sgn(fabs(occ2_re) - fabs(occ1_re));  // descending
    }
}


/*******************************************************************************
 * density_matrix_element*
 *
 * matrix elements of the a_p^+ a_q operator (in the determinantal basis)
 * these matrix elements are simply integers (since all formulas contain only
 * Kronecker delta symbols)
 *
 * encapsulates code for different sectors
 ******************************************************************************/
int density_matrix_element_0h1p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_1h0p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_1h1p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_0h2p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_2h0p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_0h0p_1h1p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_1h1p_0h0p(int p, int q, slater_det_t *bra, slater_det_t *ket);


int density_matrix_element(int sect_h, int sect_p, int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    if (sect_h == 0 && sect_p == 1) {
        return density_matrix_element_0h1p(p, q, bra, ket);
    }
    else if (sect_h == 1 && sect_p == 0) {
        return density_matrix_element_1h0p(p, q, bra, ket);
    }
    else if (sect_h == 1 && sect_p == 1) {
        if (is_vacuum_det(bra)) {
            return density_matrix_element_0h0p_1h1p(p, q, bra, ket);
        }
        else if (is_vacuum_det(ket)) {
            return density_matrix_element_1h1p_0h0p(p, q, bra, ket);
        }
        else{
            return density_matrix_element_1h1p(p, q, bra, ket);
        }
    }
    else if (sect_h == 0 && sect_p == 2) {
        return density_matrix_element_0h2p(p, q, bra, ket);
    }
    else if (sect_h == 2 && sect_p == 0) {
        return density_matrix_element_2h0p(p, q, bra, ket);
    }
    else{
        errquit("no density matrix elements implementation for the %dh%dp sector", sect_h, sect_p);
    }
}


int density_matrix_element_0h1p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int a = bra->indices[0];
    int b = ket->indices[0];

    // d_ap d_qb
    return (a == p) * (q == b);
}


int density_matrix_element_1h0p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int i = bra->indices[0];
    int j = ket->indices[0];

    // - d_iq d_jp
    return -(j == p) * (i == q);
}


int density_matrix_element_1h1p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int i = bra->indices[0];
    int a = bra->indices[1];
    int j = ket->indices[0];
    int b = ket->indices[1];

    // - d_iq d_ab d_jp
    // + d_ij d_ap d_bq

    return -(i == q) * (a == b) * (j == p) + (i == j) * (a == p) * (b == q);
}


int density_matrix_element_0h0p_1h1p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int j = ket->indices[0];
    int b = ket->indices[1];

    // {p+ q} {b+ j}

    return (p == j) * (q == b);
}


int density_matrix_element_1h1p_0h0p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int i = bra->indices[0];
    int a = bra->indices[1];

    // {i+ a} {p+ q}

    return (p == a) * (q == i);
}


int density_matrix_element_0h2p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int a = bra->indices[0];
    int b = bra->indices[1];
    int c = ket->indices[0];
    int d = ket->indices[1];

    /*
    + d_bp d_ac d_dq
    - d_bp d_ad d_cq
    - d_bc d_ap d_dq
    + d_bd d_ap d_cq
    */

    return (b == p) * (a == c) * (d == q) - (b == p) * (a == d) * (c == q) - (b == c) * (a == p) * (d == q) +
           (b == d) * (a == p) * (c == q);
}


int density_matrix_element_2h0p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int i = bra->indices[0];
    int j = bra->indices[1];
    int k = ket->indices[0];
    int l = ket->indices[1];

    /*
    - d_jq d_ik d_lp
    + d_jq d_il d_kp
    + d_jk d_iq d_lp
    - d_jl d_iq d_kp
    */

    return
            -(j == q) * (i == k) * (l == p)
            + (j == q) * (i == l) * (k == p)
            + (j == k) * (i == q) * (l == p)
            - (j == l) * (i == q) * (k == p);
}