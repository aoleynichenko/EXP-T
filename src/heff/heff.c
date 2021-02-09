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
 * heff.c
 * ======
 *
 * (1) construction and diagonalization of the effective Hamiltonian matrix;
 * (2) analysis of its eigenvalues and eigenvectors (model vectors);
 * (3) write matrix and vectors to disk.
 * NOTE: ALL THIS CODE IS DESIGNED FOR ABELIAN GROUPS ONLY!
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#include "heff.h"

#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "codata.h"
#include "io.h"
#include "engine.h"
#include "datamodel.h"
#include "linalg.h"
#include "mvcoef.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"

void zero_order_heff(int sect_h, int sect_p, size_t dim, slater_det_t *det_list, double complex *heff);
void omega_0h0p_0h1p(size_t dim, slater_det_t *det_list_1h1p, double complex *omega,
                     char *exc_oper_name, char *deexc_oper_name);
FILE *hefff_open(int sect_h, int sect_p);
void hefff_close(FILE *hefff);
void hefff_write_block(FILE *hefff, int carith, int rep_no, size_t dim, double complex *heff);
void print_model_vector(
        FILE *f_out,
        int sect_h, int sect_p,
        int rep_no, char *rep_name, int state_no,
        double complex energy,
        size_t len, double complex *coeffs, slater_det_t *det_list
);
int get_nroots_for_irrep(char *irrep_name);
void renormalize_wave_operator_0h0p_0h1p(size_t dim, slater_det_t *det_list,
                                         double complex *heff, double complex *heff_prime, double complex *omega);


/*******************************************************************************
 * get_model_space_size
 *
 * Returns total number of model-space determinants for the given FS sector
 ******************************************************************************/
size_t get_model_space_size(int sect_h, int sect_p, int nspinors, spinor_attr_t *spinor_info)
{
    int nacth, nactp;
    size_t ndeth, ndetp;

    get_active_space_size(&nacth, &nactp);

    // count number of independent determinants
    // 1 added electron:  Ndet = nactp
    // 2 added electrons: Ndet = nactp * (nactp-1)/2
    // ... (independent off-diagonal elements of n-dim matrix)

    ndeth = 1;
    ndetp = 1;

    // contribution from holes
    if (sect_h == 1) {
        ndeth = nacth;
    }
    else if (sect_h == 2) {
        ndeth = nacth * (nacth - 1) / 2;
    }
    else if (sect_h == 3) { // ?? / n! ?
        ndeth = nacth * (nacth - 1) * (nacth - 2) / 6;
    }
    else if (sect_h == 4) {
        ndeth = nacth * (nacth - 1) * (nacth - 2) * (nacth - 3) / 24;
    }

    // contibution from particles
    if (sect_p == 1) {
        ndetp = nactp;
    }
    else if (sect_p == 2) {
        ndetp = nactp * (nactp - 1) / 2;
    }
    else if (sect_p == 3) { // ?? / n! ?
        ndetp = nactp * (nactp - 1) * (nactp - 2) / 6;
    }
    else if (sect_p == 4) {
        ndetp = nactp * (nactp - 1) * (nactp - 2) * (nactp - 3) / 24;
    }

    return ndeth * ndetp;
}


/*******************************************************************************
 * detcmp
 *
 * Comparator for Slater determinants; the lower the number of irrep =>
 * the older the determinant
 ******************************************************************************/
int detcmp(const void *_d1, const void *_d2)
{
    const slater_det_t *d1 = (const slater_det_t *) _d1;
    const slater_det_t *d2 = (const slater_det_t *) _d2;
    if (d1->sym != d2->sym) {
        return abs(d1->sym) - abs(d2->sym);
    }
    else {
        return d2->sym - d1->sym;
    }
}


/*******************************************************************************
 * create_model_dets
 *
 * Returns list of model-space Slater determinants and number of determinants
 * belonging to each symmetry (via the 'ms_rep_sizes' array).
 * NOTE: detlist[] array MUST be pre-allocated
 *
 * global variables used: nspinors, nsym
 ******************************************************************************/
void create_model_dets(int sect_h, int sect_p, size_t *ms_rep_sizes, slater_det_t *detlist)
{
    moindex_t *active_holes_indices;
    moindex_t *active_parts_indices;
    int nacth, nactp;
    size_t ndet;

    // set counters to zero
    for (int irep = 0; irep < nsym; irep++) {
        ms_rep_sizes[irep] = 0;
    }

    active_holes_indices = (moindex_t *) cc_malloc(nspinors * sizeof(moindex_t));
    active_parts_indices = (moindex_t *) cc_malloc(nspinors * sizeof(moindex_t));

    get_active_space(&nacth, &nactp, active_holes_indices, active_parts_indices);

    // in order to avoid dynamically-nested loops we shall consider
    // different Fock-space sectors simply as different cases

    if (sect_h == 0 && sect_p == 0) {
        return;
    }
    else if (sect_h == 1 && sect_p == 0) {
        ndet = 0;
        for (int k = 0; k < nacth; k++) {
            moindex_t idx_k = active_holes_indices[k];
            int rep_k = spinor_info[idx_k].repno;
            detlist[ndet].indices[0] = idx_k;
            detlist[ndet].sym = rep_k;
            ms_rep_sizes[rep_k]++;
            ndet++;
        }
    }
    else if (sect_h == 0 && sect_p == 1) {
        ndet = 0;
        for (int a = 0; a < nactp; a++) {
            moindex_t idx_a = active_parts_indices[a];
            int rep_a = spinor_info[idx_a].repno;
            detlist[ndet].indices[0] = idx_a;
            detlist[ndet].sym = rep_a;
            ms_rep_sizes[rep_a]++;
            ndet++;
        }
    }
    else if (sect_h == 0 && sect_p == 2) {
        ndet = 0;
        for (int a = 0; a < nactp; a++) {
            for (int b = a + 1; b < nactp; b++) {
                moindex_t idx_a = active_parts_indices[a];
                moindex_t idx_b = active_parts_indices[b];
                detlist[ndet].indices[0] = idx_a;
                detlist[ndet].indices[1] = idx_b;

                int rep_a = spinor_info[idx_a].repno;
                int rep_b = spinor_info[idx_b].repno;
                int rep_ab = mulrep2_abelian(rep_a, rep_b);
                detlist[ndet].sym = rep_ab;
                ms_rep_sizes[rep_ab]++;
                ndet++;
            }
        }
    }
    else if (sect_h == 2 && sect_p == 0) {
        ndet = 0;
        for (int k = 0; k < nacth; k++) {
            for (int m = k + 1; m < nacth; m++) {
                moindex_t idx_k = active_holes_indices[k];
                moindex_t idx_m = active_holes_indices[m];
                detlist[ndet].indices[0] = idx_k;
                detlist[ndet].indices[1] = idx_m;

                int rep_k = spinor_info[idx_k].repno;
                int rep_m = spinor_info[idx_m].repno;
                int rep_km = mulrep2_abelian(rep_k, rep_m);
                detlist[ndet].sym = rep_km;
                ms_rep_sizes[rep_km]++;
                ndet++;
            }
        }
    }
    else if (sect_h == 1 && sect_p == 1) {
        ndet = 0;
        // add the 0h0p vacuum determinant -- if required
        if (cc_opts->mixed) {
            int rep_00 = get_totally_symmetric_irrep(); // is valid only for closed-shell references!
            detlist[0].indices[0] = 0;
            detlist[0].indices[1] = 0;
            detlist[0].sym = rep_00;
            ms_rep_sizes[rep_00]++;
            ndet++;
        }
        // 1h1p determinants
        for (int i = 0; i < nacth; i++) {
            for (int a = 0; a < nactp; a++) {
                moindex_t idx_i = active_holes_indices[i];
                moindex_t idx_a = active_parts_indices[a];
                detlist[ndet].indices[0] = idx_i;
                detlist[ndet].indices[1] = idx_a;

                int rep_i = inverse_irrep_abelian(spinor_info[idx_i].repno);
                //if (rep_i < 4) rep_i += 4;
                //else rep_i -= 4;
                int rep_a = spinor_info[idx_a].repno;
                int rep_ia = mulrep2_abelian(rep_i, rep_a);
                detlist[ndet].sym = rep_ia;
                ms_rep_sizes[rep_ia]++;
                ndet++;
            }
        }
    }
    else if (sect_h == 0 && sect_p == 3) {
        ndet = 0;
        for (int a = 0; a < nactp; a++) {
            for (int b = a + 1; b < nactp; b++) {
                for (int c = b + 1; c < nactp; c++) {
                    moindex_t idx_a = active_parts_indices[a];
                    moindex_t idx_b = active_parts_indices[b];
                    moindex_t idx_c = active_parts_indices[c];
                    detlist[ndet].indices[0] = idx_a;
                    detlist[ndet].indices[1] = idx_b;
                    detlist[ndet].indices[2] = idx_c;

                    int rep_a = spinor_info[idx_a].repno;
                    int rep_b = spinor_info[idx_b].repno;
                    int rep_c = spinor_info[idx_c].repno;
                    int rep_ab = mulrep2_abelian(rep_a, rep_b);
                    int rep_abc = mulrep2_abelian(rep_ab, rep_c);
                    detlist[ndet].sym = rep_abc;

                    ms_rep_sizes[rep_abc]++;
                    ndet++;
                }
            }
        }
    }

    // sort determinants by irrep (in ascending order)
    qsort(detlist, ndet, sizeof(slater_det_t), detcmp);

    cc_free(active_holes_indices);
    cc_free(active_parts_indices);
}


/*******************************************************************************
 * print_slater_det
 *
 * Prints Slater determinant 'det' to the stream 'f'.
 ******************************************************************************/
void print_slater_det(FILE *f, int sect_h, int sect_p, slater_det_t *det)
{
    // здесь нужно сделать исключение -- печать детерминанта из сектора 00, который у нас считается вакуумным
    // детектировать его очень просто -- det->indices[0] == 0 && det->indices[1] == 0
    // нужна отдельная функция, говорящая нам, это вакуумный детерминант или нет
    // и еще отдельная функция, создающая вакуумный детерминант
    if (is_vacuum_det(det)) {
        fprintf(f, "| vac > (%4s)\n", rep_names[det->sym]);
        return;
    }

    fprintf(f, "|");
    for (int ih = 0; ih < sect_h; ih++) {
        int idx = det->indices[ih];
        int rep = spinor_info[idx].repno;
        fprintf(f, " %4s #%4d (%12.6f)", rep_names[rep], idx + 1, spinor_info[idx].eps);
    }
    if (sect_h != 0 && sect_p != 0) {
        fprintf(f, " -> ");
    }
    for (int ip = 0; ip < sect_p; ip++) {
        int idx = det->indices[sect_h + ip];
        int rep = spinor_info[idx].repno;
        fprintf(f, " %4s #%4d (%12.6f)", rep_names[rep], idx + 1, spinor_info[idx].eps);
    }
    fprintf(f, " > (%4s)\n", rep_names[det->sym]);
}


/*******************************************************************************
 * diag_heff
 *
 * Constructs and diagonalizes effective Hamiltonian.
 * Writes Heff and right eigenvectors to the formatted files.
 * Agruments:
 *   sect_h, sect_p   number of the Fock-space sector
 *   ...              list of names of the diagrams containing 1-particle,
 *                    2-particle etc parts of the effective interaction operator
 ******************************************************************************/
void diag_heff(int sect_h, int sect_p, ...)
{
    int nacth, nactp;
    const int MAX_REPS_NUMBER = 64;
    size_t ms_size;
    size_t ms_rep_sizes[MAX_REPS_NUMBER];
    slater_det_t *dets;
    va_list vargs;
    double const coef_thresh = 1e-4;   // threshold for printing model vec-s coeff-s
    double max_energy = -1.0e9;
    int tot_sym_rep = get_totally_symmetric_irrep();
    int vac_det_rep = get_totally_symmetric_irrep();

    // for transition moments 0h0p -> 1h1p
    size_t heff_0h0p_1h1p_prime_dim;
    double complex *heff_0h0p_1h1p_prime = NULL;
    double complex *p_omega_p = NULL;
    slater_det_t *rep_dets_0011 = NULL;

    printf("\n");
    printf(" Effective Hamiltonian analysis\n");

    get_active_space_size(&nacth, &nactp);
    ms_size = get_model_space_size(sect_h, sect_p, nspinors, spinor_info);
    if (cc_opts->mixed && sect_h == 1 && sect_p == 1) {
        ms_size += 1;  // alloc space for the 0h0p determinant
    }
    printf(" Active space: %d holes, %d particles\n", nacth, nactp);
    printf(" Model space size: %d determinants (total)\n", ms_size);

    dets = (slater_det_t *) cc_malloc(ms_size * sizeof(slater_det_t));

    create_model_dets(sect_h, sect_p, ms_rep_sizes, dets);
    size_t max_heff_size = size_t_max(nsym, ms_rep_sizes);

    // TODO: to the other place!
    // for mixed-sector model: include the 0h0p determinant
    // it is supposed to belong to the totally symmetric irrep
    // (not so for open-shell references, however).
    // here we just allocate place for the 0h0p determinant

    if (cc_opts->print_level >= CC_PRINT_MEDIUM) {
        printf("\n List of model space determinants:\n");
        for (size_t i = 0; i < ms_size; i++) {
            printf(" %3d  ", i);
            print_slater_det(stdout, sect_h, sect_p, &dets[i]);
        }
    }

    printf(" Dimensions of symmetry blocks of Heff:\n");
    for (int irep = 0; irep < nsym; irep++) {
        if (ms_rep_sizes[irep] == 0) {
            continue;
        }
        printf("%4s [%d]  ", rep_names[irep], ms_rep_sizes[irep]);
    }
    printf(" (max %d)", max_heff_size);
    printf("\n");

    // allocate working arrays
    size_t nbytes = max_heff_size * max_heff_size * sizeof(double complex);
    double complex *heff = zzeros(max_heff_size, max_heff_size);
    double complex *vl = zzeros(max_heff_size, max_heff_size);
    double complex *vr = zzeros(max_heff_size, max_heff_size);
    double complex *ev = zzeros(max_heff_size, 1);
    eigval_t *eigenvalues = (eigval_t *) cc_malloc(
            (ms_size + 1) * sizeof(eigval_t)); // +1 for the reference energy in case of 1h1p
    size_t n_eigenvalues = 0;

    // not needed for the mixed-sector model 0h0p-1h1p
    if (sect_h == 1 && sect_p == 1 && cc_opts->mixed == 0) {
        eigenvalues[0].eigval = 0.0;
        eigenvalues[0].repno = get_vacuum_irrep();
        n_eigenvalues++;
    }

    // open file with formatted Heff
    FILE *hefff = hefff_open(sect_h, sect_p);

    // open unformatted file with model vectors
    int f_mvcoef = mvcoef_open(sect_h, sect_p);

    printf("\n Sector (%dh,%dp) -- analysis of model vectors (right vectors)\n", sect_h, sect_p);
    printf(" first line : irrep, state number, total energy, eigenvalue\n");
    printf(" other lines: coefficients of contributing determinants (above a threshold of %.1e)\n",
           coef_thresh);

    // for each irrep: construct Heff, diagonalize it, perform analysis of model vectors
    size_t rep_offset = 0;
    int rep0 = -1;
    for (int irep = 0; irep < nsym; irep++) {
        slater_det_t *rep_dets;

        if (ms_rep_sizes[irep] == 0) {
            continue;
        }
        if (rep0 == -1) {
            rep0 = irep;
        }

        // тут нужно проверить, не в это ли представление входит вакуумный детерминант

        // clear workspace
        memset(heff, 0, nbytes);
        memset(vl, 0, nbytes);
        memset(vr, 0, nbytes);
        memset(ev, 0, max_heff_size * sizeof(double complex));

        // get pointer to the current block of det-s with symmetry 'irep'
        rep_dets = &dets[rep_offset];
        rep_offset += ms_rep_sizes[irep];  // go to the next block of dets

        // size of the Heff subblock
        size_t ms_size = ms_rep_sizes[irep];

        // construct zero-order Hamiltonian
        zero_order_heff(sect_h, sect_p, ms_size, rep_dets, heff);

        int n_var_arg = (sect_h + 1) * (sect_p + 1) - 1; // rectangle - (0,0) sector
        va_start(vargs, sect_p);
        // loop over n-particle contibutions (n=1,2,...) to the Veff operator
        for (int ioper = 0; ioper < n_var_arg; ioper++) {
            char *dg_name = va_arg(vargs, char *);
            diagram_t *dg_veff = diagram_stack_find(dg_name);
            if (dg_veff == NULL) {
                errquit("in diag_heff(): diagram '%s' not found", dg_name);
            }
            setup_slater(dg_veff, (matrix_getter_fun) diagram_get, sect_h, sect_p, sect_h, sect_p, rank(dg_name) / 2);

            // construct matrix of the effective interaction (Veff) in the basis
            // of Slater determinants
            //
            // NOTE: here we exclude vacuum determinants!
            // They arise in the mixed-sector model, but matrix elements between
            // something and |vac> must be treated separately
            for (int i = 0; i < ms_size; i++) {
                slater_det_t *bra = &rep_dets[i];
                if (is_vacuum_det(bra)) {
                    continue;
                }
                for (int j = 0; j < ms_size; j++) {
                    slater_det_t *ket = &rep_dets[j];
                    if (is_vacuum_det(ket)) {
                        continue;
                    }
                    heff[i * ms_size + j] += slater(bra, ket);
                }
            }
        }
        va_end(vargs);

        // 0h0p-1h1p: transformation of the Heff to with the (P\OmegaP)^-1 matrix.
        // + store this 0h0p+1h1p block; it will be written to the formatted file
        // AFTER all blocks belonging to the "pure" 1h1p sector
        if (sect_h == 1 && sect_p == 1 && irep == vac_det_rep) {
            heff_0h0p_1h1p_prime_dim = ms_size + 1;
            heff_0h0p_1h1p_prime = zzeros(ms_size + 1, ms_size + 1);
            p_omega_p = zzeros(ms_size + 1, ms_size + 1);
            renormalize_wave_operator_0h0p_0h1p(ms_size, rep_dets, heff, heff_0h0p_1h1p_prime, p_omega_p);
            rep_dets_0011 = rep_dets;
        }

        // write this block of the effective Hamiltonian to the formatted file 'HEFF'
        hefff_write_block(hefff, carith, irep + 1 - rep0, ms_size, heff);

        // diagonalization:
        // 1. eigenvalues are sorted in ascending order
        // 2. eigenvectors are resorted too (and biorthonormalized)
        // 3. HEFF matrix will be destroyed
        eig(ms_size, heff, ev, vl, vr);

        // perform Loewdin orthogonalization (if not prohibited)
        if (cc_opts->do_hermit == 1) {
            loewdin_orth(ms_size, vr, vr, cc_opts->print_level == CC_PRINT_DEBUG ? 5 : 0);
            memcpy(vl, vr, ms_size * ms_size * sizeof(double complex));
        }
        // else: no orthogonalization => biorthogonal vectors => TDM_if != TDM_fi

        // get max number of roots for this irrep
        char *irrep_name = get_irrep_name(irep);
        int nroots_irep = get_nroots_for_irrep(irrep_name);
        size_t nroots = cc_opts->nroots_specified ? nroots_irep : ms_size;

        // print & save model vectors
        if (nroots > 0) {
            // print 'nroots' model vectors to stdout (formatted)
            for (size_t i = 0; i < nroots; i++) {
                print_model_vector(stdout, sect_h, sect_p,
                                   irep - rep0 + 1, irrep_name, i, ev[i],
                                   ms_size, vr + ms_size * i, rep_dets);
                // update max energy (used in eigenvalues_table() to cut unused roots)
                if (creal(ev[i]) > max_energy) {
                    max_energy = creal(ev[i]);
                }
            }
            // print model vectors to the unformatted file MVCOEF**
            mvcoef_write_vectors_unformatted(f_mvcoef, irrep_name, nroots, ms_size, rep_dets, ev, vl, vr);
        }

        // save eigenvalues
        for (size_t i = 0; i < ms_size; i++) {
            eigenvalues[n_eigenvalues + i].eigval = ev[i];
            eigenvalues[n_eigenvalues + i].repno = irep;
        }
        n_eigenvalues += ms_size;
    }

    qsort(eigenvalues, n_eigenvalues, sizeof(eigval_t), eigval_cmp);
    eigenvalues_table(n_eigenvalues, eigenvalues, cc_opts->degen_thresh, max_energy);

    double eigval_0 = creal(eigenvalues[0].eigval);
    if (sect_h != sect_p) {
        printf("\n Ionization potential wrt reference state = %18.12f a.u. = %8.4f eV = %10.2f cm^-1\n",
               -eigval_0, -eigval_0 * 27.21138602, -eigval_0 * 219474.6313702);
    }

    // write the 0h0p+1h1p transformed Heff block
    // together with the transformation P \Omega P matrix
    if (sect_h == 1 && sect_p == 1) {
        fprintf(hefff, "0011         # sector\n");
        // heff block
        hefff_write_block(hefff, carith, vac_det_rep + 1 - rep0, heff_0h0p_1h1p_prime_dim, heff_0h0p_1h1p_prime);
        // P \omega P
        size_t dim = heff_0h0p_1h1p_prime_dim;
        for (int i = 0; i < dim * dim; i++) {
            if (carith) {
                fprintf(hefff, "%21.12E%21.12E", creal(p_omega_p[i]), cimag(p_omega_p[i]));
                if (i > 0 && i % 2 != 0) { fprintf(hefff, "\n"); }
            }
            else {
                fprintf(hefff, "%21.12E", creal(p_omega_p[i]));
                if (i > 0 && (i + 1) % 4 == 0) { fprintf(hefff, "\n"); }
            }
        }
        if (dim * dim % 2 != 0) { fprintf(hefff, "\n"); }

        // diagonalize the 00-11 block
        double complex *vl_0011 = zzeros(dim, dim);
        double complex *vr_0011 = zzeros(dim, dim);
        double complex *ev_0011 = zzeros(dim, 1);
        int irep = vac_det_rep;
        char *irrep_name = get_irrep_name(irep);
        int nroots_irep = get_nroots_for_irrep(irrep_name);
        size_t nroots = cc_opts->nroots_specified ? nroots_irep : dim;

        eig(dim, heff_0h0p_1h1p_prime, ev_0011, vl_0011, vr_0011);
        if (cc_opts->do_hermit == 1) {
            loewdin_orth(dim, vr_0011, vr_0011, cc_opts->print_level == CC_PRINT_DEBUG ? 5 : 0);
            memcpy(vl_0011, vr_0011, dim * dim * sizeof(double complex));
        }

        // print model vectors to the unformatted file MVCOEF**
        slater_det_t *det_list_0011 = (slater_det_t *) cc_malloc(sizeof(slater_det_t) * dim);
        det_list_0011[0].indices[0] = 0;
        det_list_0011[0].indices[1] = 0;
        det_list_0011[0].sym = vac_det_rep;
        memcpy(det_list_0011 + 1, rep_dets_0011, sizeof(slater_det_t) * (dim - 1));
        ev_0011[0] = 0.0 + 0.0 * I;

        int f_mvcoef_0011 = io_open("MVCOEF0011", "w");
        mvcoef_write_vectors_unformatted(f_mvcoef_0011, rep_names[irep], nroots + 1, dim, det_list_0011, ev_0011,
                                         vl_0011, vr_0011);

        mvcoef_close(f_mvcoef_0011, eigval_0);
        cc_free(vl_0011);
        cc_free(vr_0011);
        cc_free(ev_0011);
    }

    mvcoef_close(f_mvcoef, eigval_0);
    hefff_close(hefff);
    cc_free(eigenvalues);
    cc_free(heff);
    cc_free(ev);
    cc_free(vr);
    cc_free(vl);
    cc_free(dets);
    if (heff_0h0p_1h1p_prime != NULL) {
        cc_free(heff_0h0p_1h1p_prime);
    }
    if (p_omega_p != NULL) {
        cc_free(p_omega_p);
    }

    // calculation of density matrices and natural orbitals
    // for the target sector only
    if (cc_opts->sector_h == sect_h && cc_opts->sector_p == sect_p &&
        cc_opts->n_denmat) {
        for (int ipair = 0; ipair < cc_opts->n_denmat; ipair++) {
            int sect2_h = cc_opts->denmat_query[ipair].sect2[0];
            int sect2_p = cc_opts->denmat_query[ipair].sect2[1];
            if (!(sect2_h == sect_h && sect2_p == sect_p)) {
                continue;
            }
            int rep1 = get_rep_number(cc_opts->denmat_query[ipair].rep1_name);
            int rep2 = get_rep_number(cc_opts->denmat_query[ipair].rep2_name);
            int state1 = cc_opts->denmat_query[ipair].state1;
            int state2 = cc_opts->denmat_query[ipair].state2;
            // single code for NOs and NTOs:
            // rep1 == rep2 && state1 == state2  =>  single state (NO)
            // else: pair of states (NTO)
            density_matrix(sect_h, sect_p, rep1, state1, rep2, state2);
        }
    }

    // transition dipole moments via the DL-TDM techniques
    // for the target sector only
    if (cc_opts->sector_h == sect_h && cc_opts->sector_p == sect_p &&
        cc_opts->do_diplen_tdm) {
        dipole_length_tdms(sect_h, sect_p);
    }

    // estimate properties: only for the target sector
    if (cc_opts->sector_h == sect_h && cc_opts->sector_p == sect_p &&
        cc_opts->n_ms_props > 0) {
        for (int i = 0; i < cc_opts->n_ms_props; i++) {
            model_space_property(cc_opts->prop_queries + i);
        }
    }
}


/**
 * construct zero-order effective Hamiltonian (H0) matrix:
 *  - only the block belonging to one irrep;
 *  - only one-electron spinor energies are used.
 * It is a diagonal matrix: heff[i,i] = sum{particles} eps - sum{holes} eps
 *
 * @param sect_h, sect_p  Fock space sector signature
 * @param dim dimension of this Heff block (number of model det-s of one symmetry)
 * @param det_list list of model space determinants (with one symmetry)
 * @param heff Heff matrix, of size dim x dim
 */
void zero_order_heff(int sect_h, int sect_p, size_t dim, slater_det_t *det_list, double complex *heff)
{
    for (size_t i = 0; i < dim; i++) {
        slater_det_t *d = det_list + i;

        double H0 = 0.0;
        for (int h = 0; h < sect_h; h++) {
            H0 -= spinor_info[d->indices[h]].eps;
        }
        for (int p = 0; p < sect_p; p++) {
            H0 += spinor_info[d->indices[sect_h + p]].eps;
        }
        heff[i * dim + i] = H0 + 0.0 * I;
    }
}


/*******************************************************************************
 * eigval_cmp
 *
 * Comparator for eigenvalues (for qsort)
 ******************************************************************************/
int eigval_cmp(const void *aa, const void *bb)
{
    const eigval_t *a = aa, *b = bb;
    return (creal(a->eigval) < creal(b->eigval)) ? -1 : (creal(a->eigval) > creal(b->eigval));
}


/*******************************************************************************
 * eigenvalues_table
 *
 * Prints nice table with all information about eigenvalues and symmetries of
 * solutions
 ******************************************************************************/
void eigenvalues_table(size_t n_eigenvalues, eigval_t *eigenvalues, double degen_thresh, double max_energy)
{
    size_t i, j;
    size_t ilevel;
    int irep;
    int rep_deg;
    int n_nz_reps;
    int nz_reps[64]; // irreps with non-zero number of electronic states
    double complex e0 = eigenvalues[0].eigval;
    double complex ei;
    double rel_energy;
    double abs_energy;
    double rel_energy_ev;
    double rel_energy_cm;
    double eref = cc_opts->eref;   // (see options.h)

    // find reps with nonzero Heff and their number
    n_nz_reps = 0;
    for (i = 0; i < nsym; i++) {
        nz_reps[i] = 0;
    }
    for (i = 0; i < n_eigenvalues; i++) {
        nz_reps[eigenvalues[i].repno]++;
    }
    for (i = 0; i < nsym; i++) {
        if (nz_reps[i] != 0) {
            nz_reps[n_nz_reps++] = i;
        }
    }

    printf("\n Heff eigenvalues:\n (degeneracy threshold = %.1e a.u.)\n\n", degen_thresh);

    // print table header (see symmetry.h for 'rep_names[]')
    printf(" Level  Re(eigenvalue)  Im(eigv)               Abs energy  Rel eigenvalue    Rel eigv, eV  Rel eigv, cm-1  deg  symmetry\n");
    /*for (irep = 0; irep < n_nz_reps; irep++) {
        printf("%4s|", rep_names[nz_reps[irep]]);
    }*/
    printf(" ------------------------------------------------------------------------------------------------------------------------\n");

    // print table with energy levels
    // only energy levels with eigenvalue Ei <= max_energy+1e-5 will be printed
    ilevel = 1;
    for (i = 0; i < n_eigenvalues; i++) {
        if (creal(eigenvalues[i].eigval) > max_energy + 1.0e-5) {
            break;
        }

        // calculate degeneracy
        int deg = 0;
        for (j = i; j < n_eigenvalues; j++) {
            if (cabs(eigenvalues[i].eigval - eigenvalues[j].eigval) > degen_thresh) {
                break;
            }
            else {
                deg++;
            }
        }

        // print info about (maybe degenerate) level
        ei = eigenvalues[i].eigval;
        abs_energy = eref + creal(ei);
        rel_energy = creal(ei) - creal(e0);
        rel_energy_cm = rel_energy * CODATA_AU_TO_CM;
        rel_energy_ev = rel_energy * CODATA_AU_TO_EV;
        printf("@%5d%16.10f%10.2e%25.17f%16.10f%16.10f%16.6f  %2d  ",
               ilevel, creal(ei), cimag(ei), abs_energy, rel_energy,
               rel_energy_ev, rel_energy_cm, deg);
        for (irep = 0; irep < n_nz_reps; irep++) {
            rep_deg = 0;
            for (j = 0; j < deg; j++) {
                if (eigenvalues[i + j].repno == nz_reps[irep]) {
                    rep_deg++;
                }
            }
            if (rep_deg == 0) { continue; }
            printf(" %s", rep_names[nz_reps[irep]]);
            if (rep_deg > 1) {
                printf("(%d)", rep_deg);
            }
        }
        printf("\n");

        // go to the next (maybe degenerate) energy level
        i += deg - 1;
        ilevel++;
    }
}


/**
 * Prints model vectors.
 */
void print_model_vector(
        FILE *f_out,
        int sect_h, int sect_p,
        int rep_no, char *rep_name, int state_no,
        double complex energy,
        size_t len, double complex *coeffs, slater_det_t *det_list
)
{
    double complex coef;
    double const coef_thresh = 1e-4;   // threshold for printing model vec-s coeff-s

    fprintf(f_out, "\n");
    fprintf(f_out, " Irrep %d (%4s) State %d Energy %23.15f Eigenvalue %14.8f%14.4E\n",
            rep_no, rep_name, state_no + 1, cc_opts->eref + creal(energy), creal(energy), cimag(energy));
    for (size_t j = 0; j < len; j++) {
        coef = coeffs[j];
        if (cabs(coef) < coef_thresh) {
            continue;
        }
        fprintf(f_out, " %10.5f%10.5f ", creal(coef), cimag(coef));
        print_slater_det(f_out, sect_h, sect_p, det_list + j);
    }
}


int get_nroots_for_irrep(char *irrep_name)
{
    int nroots_irep = 0;

    for (int ii = 0; ii < nsym; ii++) {
        if (strcmp(irrep_name, cc_opts->nroots_specs.rep_names[ii]) == 0) {
            nroots_irep = cc_opts->nroots_specs.dim[ii];
            break;
        }
    }

    return nroots_irep;
}


/*
 * For transition moments between different Fock space sectors.
 * When using this code, please, cite the following paper:
 *
 * A. Zaitsevskii, A. V. Oleynichenko, E. Eliav,
 * Finite-field calculations of transition properties by the fock space
 * relativistic coupled cluster method: transitions between different
 * Fock space sectors.
 * Symmetry, 2020 (submitted).
 */


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

    // < 0h0p | Omega | 0h0p > = 1
    omega[0] = 1.0 + 0.0 * I;

    // < 0h0p | Omega | 1h1p > = S1[a,i]
    for (size_t iket = 0; iket < dim; iket++) {
        slater_det_t *det_ia = det_list_1h1p + iket;
        int i = det_ia->indices[0];
        int a = det_ia->indices[1];
        //omega[iket+1] = diagram_get_2(dg_deexc, a, i);
        omega[(dim + 1) * (iket + 1)] = diagram_get_2(dg_deexc, a, i);
    }

    // < 1h1p | Omega | 0h0p > = T1[i,a]
    for (size_t ibra = 0; ibra < dim; ibra++) {
        slater_det_t *det_ia = det_list_1h1p + ibra;
        int i = det_ia->indices[0];
        int a = det_ia->indices[1];
        //omega[(dim+1)*(ibra+1)] = diagram_get_2(dg_exc, i, a);
        omega[ibra + 1] = diagram_get_2(dg_exc, i, a);
    }

    // < 1h1p i,a | Omega | 1h1p j,b > = \delta_ij \delta_ab + T1[i,a] * S1[j,b]
    for (size_t ibra = 0; ibra < dim; ibra++) {
        for (size_t iket = 0; iket < dim; iket++) {
            slater_det_t *det_ia = det_list_1h1p + ibra;
            slater_det_t *det_jb = det_list_1h1p + iket;
            int i = det_ia->indices[0];
            int a = det_ia->indices[1];
            int j = det_jb->indices[0];
            int b = det_jb->indices[1];
            //omega[(dim+1) * (ibra+1) + iket + 1] = (i==j)*(a==b) + diagram_get_2(dg_exc, i, a) * diagram_get_2(dg_deexc, b, j);
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

    // allocate working arrays, init with zeros
    memset(omega, 0, sizeof(double complex) * dimx * dimx);
    double complex *omega_inv = zzeros(dimx, dimx);
    double complex *heff_0h0p_1h1p = zzeros(dimx, dimx);
    double complex *buf = zzeros(dimx, dimx);

    // construct 0h0p+1h1p (extended) block-diagonal Heff matrix
    // Heff 1h1p block will be placed into the right lower angle.
    // Ecorr is subtracted from the whole diagonal.
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            heff_0h0p_1h1p[dimx * (i + 1) + j + 1] = heff[dim * i + j];
        }
    }

    // construct P \Omega P matrix (P=0h0p+1h1p)
    omega_0h0p_0h1p(dim, det_list, omega, "t1c", "e1c");

    // construct inverse matrix (P \Omega P)^{-1}
    inv(dimx, omega, omega_inv);

    // Heff' = { (P\OmegaP) * { Heff * (P\OmegaP)^-1 } }
    // two matrix multiplications
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                dimx, dimx, dimx, &alpha, heff_0h0p_1h1p, dimx, omega_inv, dimx, &beta, buf, dimx);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                dimx, dimx, dimx, &alpha, omega, dimx, buf, dimx, &beta, heff_prime, dimx);

    // deallocate working arrays
    cc_free(omega_inv);
    cc_free(heff_0h0p_1h1p);
    cc_free(buf);
}

