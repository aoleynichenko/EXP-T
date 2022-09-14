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
 * Calculation of density matrices and natural orbitals.
 * Associated header file: heff.h
 *
 * 2019-2022 Alexander Oleynichenko
 */

#include <complex.h>
#include <string.h>

#include "slater_rules.h"
#include "spinors.h"
#include "symmetry.h"
#include "options.h"


// functions used in this module
int cmp_natorb_occupations(double complex occ1, double complex occ2);

void write_natural_spinors_formatted(char *natorb_file_name, int sect_h, int sect_p,
                                     int dim, double *occ_numbers, double complex *nat_orbs,
                                     double occ_thresh);

int density_matrix_element(int sect_h, int sect_p, int p, int q, slater_det_t *bra, slater_det_t *ket);

void sort_vectors(int n, double complex *ev, double complex *vl, double complex *vr,
                  int (*cmp)(double complex, double complex));


/*
 * @param dm (output) density matrix, size n_active x n_active
 * @param dim_dm (output) n_active
 */
void construct_model_space_density_matrix(int sect_h, int sect_p,
                                          int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
                                          int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
                                          double complex *dm, int *dim_dm)
{
    // only 'valence' parts of DMs are be constructed
    // dim DMs = n_active x n_active
    size_t n_active = 0;
    int active_spinors[CC_MAX_SPINORS]; // local -> global spinor index mapping
    int nspinors = get_num_spinors();

    // TODO: refactor, separate utility function (move to spinors.c)
    for (int i = 0; i < nspinors; i++) {
        if ((is_act_hole(i) && sect_h > 0) || (is_act_particle(i) && sect_p > 0)) {
            active_spinors[n_active] = i;
            n_active++;
        }
    }

    for (int p = 0; p < n_active; p++) {
        for (int q = 0; q < n_active; q++) {
            double complex d_pq = 0.0 + 0.0 * I;
            // loops over model vectors
            for (size_t i = 0; i < dim_bra; i++) {
                for (size_t j = 0; j < dim_ket; j++) {
                    slater_det_t *bra = &dets_bra[i];
                    slater_det_t *ket = &dets_ket[j];
                    int bra_d_ket = density_matrix_element(sect_h, sect_p, active_spinors[p], active_spinors[q], bra,
                                                           ket);
                    d_pq += conj(coef_bra[i]) * coef_ket[j] * bra_d_ket;
                }
            }
            dm[p * n_active + q] = d_pq;
        }
    }
}


/**
 * density_matrix_element*
 *
 * matrix elements of the a_p^+ a_q operator (in the determinant basis)
 * these matrix elements are simply integers (since all formulas contain only
 * Kronecker delta symbols)
 *
 * encapsulates code for different sectors
 **/
int density_matrix_element_0h1p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_1h0p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_1h1p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_0h2p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_2h0p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_3h0p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_0h0p_1h1p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_1h1p_0h0p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_0h3p(int p, int q, slater_det_t *bra, slater_det_t *ket);

int density_matrix_element_1h2p(int p, int q, slater_det_t *bra, slater_det_t *ket);


int density_matrix_element(int sect_h, int sect_p, int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    if (sect_h == 0 && sect_p == 1) {
        return density_matrix_element_0h1p(p, q, bra, ket);
    }
    else if (sect_h == 1 && sect_p == 0) {
        return density_matrix_element_1h0p(p, q, bra, ket);
    }
    else if (sect_h == 1 && sect_p == 1) {
        if (is_vacuum_det(bra) && is_vacuum_det(ket)) {
            // < vac | ap^+ ap | vac > = occ(p)
            if (p != q) {
                return 0;
            }
            if (is_act_hole(p)) {
                return 1;
            }
            if (is_act_particle(p)) {
                return 0;
            }
        }
        else if (is_vacuum_det(bra)) {
            return density_matrix_element_0h0p_1h1p(p, q, bra, ket);
        }
        else if (is_vacuum_det(ket)) {
            return density_matrix_element_1h1p_0h0p(p, q, bra, ket);
        }
        else {
            return density_matrix_element_1h1p(p, q, bra, ket);
        }
    }
    else if (sect_h == 0 && sect_p == 2) {
        return density_matrix_element_0h2p(p, q, bra, ket);
    }
    else if (sect_h == 0 && sect_p == 3) {
        return density_matrix_element_0h3p(p, q, bra, ket);
    }
    else if (sect_h == 2 && sect_p == 0) {
        return density_matrix_element_2h0p(p, q, bra, ket);
    }
    else if (sect_h == 3 && sect_p == 0) {
        return density_matrix_element_3h0p(p, q, bra, ket);
    }
    else if (sect_h == 1 && sect_p == 2) {
        return density_matrix_element_1h2p(p, q, bra, ket);
    }
    else {
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


// bra = vac
int density_matrix_element_0h0p_1h1p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int j = ket->indices[0];
    int b = ket->indices[1];

    // {p+ q} {b+ j}

    return (p == j) * (q == b);
}


// ket = vac
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

    return -(j == q) * (i == k) * (l == p)
           + (j == q) * (i == l) * (k == p)
           + (j == k) * (i == q) * (l == p)
           - (j == l) * (i == q) * (k == p);
}


int density_matrix_element_3h0p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int i = bra->indices[0];
    int j = bra->indices[1];
    int k = ket->indices[2];
    int l = ket->indices[0];
    int m = ket->indices[1];
    int n = ket->indices[2];

    /*
    + d_kq d_jl d_im d_np
    - d_kq d_jl d_in d_mp
    - d_kq d_jm d_il d_np
    + d_kq d_jm d_in d_lp
    + d_kq d_jn d_il d_mp
    - d_kq d_jn d_im d_lp
    - d_kl d_jq d_im d_np
    + d_kl d_jq d_in d_mp
    + d_kl d_jm d_iq d_np
    - d_kl d_jn d_iq d_mp
    + d_km d_jq d_il d_np
    - d_km d_jq d_in d_lp
    - d_km d_jl d_iq d_np
    + d_km d_jn d_iq d_lp
    - d_kn d_jq d_il d_mp
    + d_kn d_jq d_im d_lp
    + d_kn d_jl d_iq d_mp
    - d_kn d_jm d_iq d_lp
    */

    return +(k == q) * (j == l) * (i == m) * (n == p)
           - (k == q) * (j == l) * (i == n) * (m == p)
           - (k == q) * (j == m) * (i == l) * (n == p)
           + (k == q) * (j == m) * (i == n) * (l == p)
           + (k == q) * (j == n) * (i == l) * (m == p)
           - (k == q) * (j == n) * (i == m) * (l == p)
           - (k == l) * (j == q) * (i == m) * (n == p)
           + (k == l) * (j == q) * (i == n) * (m == p)
           + (k == l) * (j == m) * (i == q) * (n == p)
           - (k == l) * (j == n) * (i == q) * (m == p)
           + (k == m) * (j == q) * (i == l) * (n == p)
           - (k == m) * (j == q) * (i == n) * (l == p)
           - (k == m) * (j == l) * (i == q) * (n == p)
           + (k == m) * (j == n) * (i == q) * (l == p)
           - (k == n) * (j == q) * (i == l) * (m == p)
           + (k == n) * (j == q) * (i == m) * (l == p)
           + (k == n) * (j == l) * (i == q) * (m == p)
           - (k == n) * (j == m) * (i == q) * (l == p);
}


int density_matrix_element_0h3p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int a = bra->indices[0];
    int b = bra->indices[1];
    int c = bra->indices[2];
    int d = ket->indices[0];
    int e = ket->indices[1];
    int f = ket->indices[2];

    /*
    - d_cp d_bd d_ae d_fq
    + d_cp d_bd d_af d_eq
    + d_cp d_be d_ad d_fq
    - d_cp d_be d_af d_dq
    - d_cp d_bf d_ad d_eq
    + d_cp d_bf d_ae d_dq
    + d_cd d_bp d_ae d_fq
    - d_cd d_bp d_af d_eq
    - d_cd d_be d_ap d_fq
    + d_cd d_bf d_ap d_eq
    - d_ce d_bp d_ad d_fq
    + d_ce d_bp d_af d_dq
    + d_ce d_bd d_ap d_fq
    - d_ce d_bf d_ap d_dq
    + d_cf d_bp d_ad d_eq
    - d_cf d_bp d_ae d_dq
    - d_cf d_bd d_ap d_eq
    + d_cf d_be d_ap d_dq
    */

    return -(c == p) * (b == d) * (a == e) * (f == q)
           + (c == p) * (b == d) * (a == f) * (e == q)
           + (c == p) * (b == e) * (a == d) * (f == q)
           - (c == p) * (b == e) * (a == f) * (d == q)
           - (c == p) * (b == f) * (a == d) * (e == q)
           + (c == p) * (b == f) * (a == e) * (d == q)
           + (c == d) * (b == p) * (a == e) * (f == q)
           - (c == d) * (b == p) * (a == f) * (e == q)
           - (c == d) * (b == e) * (a == p) * (f == q)
           + (c == d) * (b == f) * (a == p) * (e == q)
           - (c == e) * (b == p) * (a == d) * (f == q)
           + (c == e) * (b == p) * (a == f) * (d == q)
           + (c == e) * (b == d) * (a == p) * (f == q)
           - (c == e) * (b == f) * (a == p) * (d == q)
           + (c == f) * (b == p) * (a == d) * (e == q)
           - (c == f) * (b == p) * (a == e) * (d == q)
           - (c == f) * (b == d) * (a == p) * (e == q)
           + (c == f) * (b == e) * (a == p) * (d == q);
}


int density_matrix_element_1h2p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int i = bra->indices[0];
    int a = bra->indices[1];
    int b = bra->indices[2];
    int j = ket->indices[0];
    int c = ket->indices[1];
    int d = ket->indices[2];

    /*
    + d_iq d_bc d_ad d_jp
    - d_iq d_bd d_ac d_jp
    + d_ij d_bp d_ac d_dq
    - d_ij d_bp d_ad d_cq
    - d_ij d_bc d_ap d_dq
    + d_ij d_bd d_ap d_cq
    */

    return +(i == q) * (b == c) * (a == d) * (j == p)
           - (i == q) * (b == d) * (a == c) * (j == p)
           + (i == j) * (b == p) * (a == c) * (d == q)
           - (i == j) * (b == p) * (a == d) * (c == q)
           - (i == j) * (b == c) * (a == p) * (d == q)
           + (i == j) * (b == d) * (a == p) * (c == q);
}


int density_matrix_element_0h1p_1h2p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int a = bra->indices[0];
    int i = ket->indices[0];
    int b = ket->indices[1];
    int c = ket->indices[2];

    /*
    + d_ab d_ip d_cq
    - d_ac d_ip d_bq
    */

    return + (a == b) * (i == p) * (c == q)
           - (a == c) * (i == p) * (b == q);
}


int density_matrix_element_1h2p_0h1p(int p, int q, slater_det_t *bra, slater_det_t *ket)
{
    int i = bra->indices[0];
    int a = bra->indices[1];
    int b = bra->indices[2];
    int c = ket->indices[0];

    /*
    + d_iq d_bp d_ac
    - d_iq d_bc d_ap
    */

    return + (i == q) * (b == p) * (a == c)
           - (i == q) * (b == c) * (a == p);
}
