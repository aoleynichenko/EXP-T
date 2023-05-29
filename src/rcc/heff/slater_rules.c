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
 * Slater rules: evaluation of matrix elements in the basis of Slater
 * determinants.
 *
 * 2019-2021 Alexander Oleynichenko
 */

#include <stdlib.h>
#include <stdio.h>

#include "slater_rules.h"

void *source_matrix;

double complex (*get_element)(void *source, void *indices);

double complex (*slater_rule)(slater_det_t *bra_det, slater_det_t *ket_det) = NULL;

double complex slater_01_1_01(slater_det_t *bra, slater_det_t *ket);

double complex slater_10_1_10(slater_det_t *bra, slater_det_t *ket);

double complex slater_02_1_02(slater_det_t *bra, slater_det_t *ket);

double complex slater_02_2_02(slater_det_t *bra, slater_det_t *ket);

double complex slater_20_1_20(slater_det_t *bra, slater_det_t *ket);

double complex slater_20_2_20(slater_det_t *bra, slater_det_t *ket);

double complex slater_00_1_11(slater_det_t *bra, slater_det_t *ket);

double complex slater_11_1_00(slater_det_t *bra, slater_det_t *ket);

double complex slater_11_1_11(slater_det_t *bra, slater_det_t *ket);

double complex slater_11_2_11(slater_det_t *bra, slater_det_t *ket);

double complex slater_03_1_03(slater_det_t *bra, slater_det_t *ket);

double complex slater_03_2_03(slater_det_t *bra, slater_det_t *ket);

double complex slater_03_3_03(slater_det_t *bra, slater_det_t *ket);

double complex slater_30_1_30(slater_det_t *bra, slater_det_t *ket);

double complex slater_30_2_30(slater_det_t *bra, slater_det_t *ket);

double complex slater_30_3_30(slater_det_t *bra, slater_det_t *ket);

double complex slater_12_1_12(slater_det_t *bra, slater_det_t *ket);

double complex slater_12_2_12(slater_det_t *bra, slater_det_t *ket);

double complex slater_12_3_12(slater_det_t *bra, slater_det_t *ket);

double complex slater_01_1_12(slater_det_t *bra, slater_det_t *ket);

double complex slater_12_1_01(slater_det_t *bra, slater_det_t *ket);


/**
 * Setups module for evaluation of matrix elements in the basis of Slater
 * determinants. To use this routine, you must define the 'getter' interface to
 * the 'source' object. For example, use can use diagram_t type object with the
 * 'diagram_get' routine or ordinary matrix with the appropriate accessor func-n
 * Must be invoked before matrix construction.
 * Arguments:
 *   source    matrix-like object containing matrix elements in spin-orbital
 *             (spinor) basis
 *   getter    function operating at the 'source' object, extracts matrix
 *             element at index 'indices'
 *   sect_h, sect_p number of Fock-space sector
 *   npart     "number of particles" of the operator consider (1 of one-particle
 *             operator, 2 for two-particle etc)
 */
void setup_slater(void *source, double complex (*getter)(void *source, void *indices),
                  int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p, int npart)
{
    source_matrix = source;
    get_element = getter;

    slater_rule = NULL;

    /*
     * 0h1p - 0h1p
     */
    if (bra_sect_h == 0 && bra_sect_p == 1 && ket_sect_h == 0 && ket_sect_p == 1) {
        if (npart == 1) {
            slater_rule = slater_01_1_01;
        }
    }

    /*
     * 1h0p - 1h0p
     */
    if (bra_sect_h == 1 && bra_sect_p == 0 && ket_sect_h == 1 && ket_sect_p == 0) {
        if (npart == 1) {
            slater_rule = slater_10_1_10;
        }
    }

    /*
     * 0h0p - 1h1p
     */
    if (bra_sect_h == 0 && bra_sect_p == 0 && ket_sect_h == 1 && ket_sect_p == 1) {
        if (npart == 1) {
            slater_rule = slater_00_1_11;
        }
    }

    /*
     * 1h1p - 0h0p
     */
    if (bra_sect_h == 1 && bra_sect_p == 1 && ket_sect_h == 0 && ket_sect_p == 0) {
        if (npart == 1) {
            slater_rule = slater_11_1_00;
        }
    }

    /*
     * 1h1p - 1h1p
     */
    if (bra_sect_h == 1 && bra_sect_p == 1 && ket_sect_h == 1 && ket_sect_p == 1) {
        if (npart == 1) {
            slater_rule = slater_11_1_11;
        }
        else if (npart == 2) {
            slater_rule = slater_11_2_11;
        }
    }

    /*
     * 0h2p - 0h2p
     */
    if (bra_sect_h == 0 && bra_sect_p == 2 && ket_sect_h == 0 && ket_sect_p == 2) {
        if (npart == 1) {
            slater_rule = slater_02_1_02;
        }
        else if (npart == 2) {
            slater_rule = slater_02_2_02;
        }
    }

    /*
     * 2h0p - 2h0p
     */
    if (bra_sect_h == 2 && bra_sect_p == 0 && ket_sect_h == 2 && ket_sect_p == 0) {
        if (npart == 1) {
            slater_rule = slater_20_1_20;
        }
        else if (npart == 2) {
            slater_rule = slater_20_2_20;
        }
    }

    /*
     * 0h3p - 0h3p
     */
    if (bra_sect_h == 0 && bra_sect_p == 3 && ket_sect_h == 0 && ket_sect_p == 3) {
        if (npart == 1) {
            slater_rule = slater_03_1_03;
        }
        else if (npart == 2) {
            slater_rule = slater_03_2_03;
        }
        else if (npart == 3) {
            slater_rule = slater_03_3_03;
        }
    }

    /*
     * 3h0p - 3h0p
     */
    if (bra_sect_h == 3 && bra_sect_p == 0 && ket_sect_h == 3 && ket_sect_p == 0) {
        if (npart == 1) {
            slater_rule = slater_30_1_30;
        }
        else if (npart == 2) {
            slater_rule = slater_30_2_30;
        }
        else if (npart == 3) {
            slater_rule = slater_30_3_30;
        }
    }

    /*
     * 1h2p - 1h2p
     */
    if (bra_sect_h == 1 && bra_sect_p == 2 && ket_sect_h == 1 && ket_sect_p == 2) {
        if (npart == 1) {
            slater_rule = slater_12_1_12;
        }
        else if (npart == 2) {
            slater_rule = slater_12_2_12;
        }
        else if (npart == 3) {
            slater_rule = slater_12_3_12;
        }
    }

    /*
     * 0h1p - 1h2p
     */
    if (bra_sect_h == 0 && bra_sect_p == 1 && ket_sect_h == 1 && ket_sect_p == 2) {
        if (npart == 1) {
            slater_rule = slater_01_1_12;
        }
    }

    /*
     * 1h2p - 0h1p
     */
    if (bra_sect_h == 1 && bra_sect_p == 2 && ket_sect_h == 0 && ket_sect_p == 1) {
        if (npart == 1) {
            slater_rule = slater_12_1_01;
        }
    }

    if (slater_rule == NULL) {
        printf("Cannot evaluate matrix element of a %d-particle operator "
               "between Slater determinants of type |%dh%dp> and |%dh%dp>\n",
               npart, bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p);
        printf("see %s for details, line %d\n", __FILE__, __LINE__);
        exit(0);
    }
}


/*******************************************************************************
 * Slater rules for different sectors and n-particle operators
 ******************************************************************************/

/*
 * Sector 0h1p
 * 1-particle operator
 */
double complex slater_01_1_01(slater_det_t *d1, slater_det_t *d2)
{
    int idx[2];
    idx[0] = d1->indices[0];
    idx[1] = d2->indices[0];
    return get_element(source_matrix, idx);
}


/*
 * Sector 1h0p
 * 1-particle operator
 */
double complex slater_10_1_10(slater_det_t *d1, slater_det_t *d2)
{
    int idx[2];
    idx[1] = d1->indices[0];
    idx[0] = d2->indices[0];
    return -get_element(source_matrix, idx);
}


/*
 * Sectors 0h0p - 1h1p
 * 1-particle operator
 */
double complex slater_00_1_11(slater_det_t *d1, slater_det_t *d2)
{
    int idx[2];
    idx[0] = d2->indices[1];
    idx[1] = d2->indices[0];

    return get_element(source_matrix, idx);
}


/*
 * Sectors 1h1p - 0h0p
 * 1-particle operator
 */
double complex slater_11_1_00(slater_det_t *d1, slater_det_t *d2)
{
    int idx[2];
    idx[0] = d1->indices[0];
    idx[1] = d1->indices[1];

    return get_element(source_matrix, idx);
}


/*
 * Sector 1h1p
 * 1-particle operator
 */
double complex slater_11_1_11(slater_det_t *bra, slater_det_t *ket)
{
    int idx2[2];
    int i = bra->indices[0];
    int a = bra->indices[1];
    int j = ket->indices[0];
    int b = ket->indices[1];
    double complex matr_elem = 0.0 + 0.0 * I;

    // (1) - 1.0 heff1 [ j i ] d_ab
    if (a == b) {
        idx2[0] = j;
        idx2[1] = i;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (2) + 1.0 heff1 [ a b ] d_ij
    if (i == j) {
        idx2[0] = a;
        idx2[1] = b;
        matr_elem += get_element(source_matrix, idx2);
    }

    return matr_elem;
}


/*
 * Sector 1h1p
 * 2-particle operator
 */
double complex slater_11_2_11(slater_det_t *bra, slater_det_t *ket)
{
    int idx4[4];
    int i = bra->indices[0];
    int a = bra->indices[1];
    int j = ket->indices[0];
    int b = ket->indices[1];
    double complex matr_elem = 0.0 + 0.0 * I;

    // (1) - heff2 [ j a i b ]
    // = - a j b i
    idx4[0] = a;
    idx4[1] = j;
    idx4[2] = b;
    idx4[3] = i;
    matr_elem -= get_element(source_matrix, idx4);
    //printf("%d %d %d %d = %f\n", j, a, b, i, creal(matr_elem));

    return matr_elem;
}


/*
 * Sector 0h2p
 * 1-particle operator
 */
double complex slater_02_1_02(slater_det_t *bra, slater_det_t *ket)
{
    int idx2[2];
    int a = bra->indices[0];
    int b = bra->indices[1];
    int c = ket->indices[0];
    int d = ket->indices[1];
    double complex matr_elem = 0.0 + 0.0 * I;
    /*
    (0) + 1.0 heff1 [ b d ] d_ac
    (1) - 1.0 heff1 [ b c ] d_ad
    (2) - 1.0 heff1 [ a d ] d_bc
    (3) + 1.0 heff1 [ a c ] d_bd
    */
    if (a == c) {
        idx2[0] = b;
        idx2[1] = d;
        matr_elem += get_element(source_matrix, idx2);
    }
    if (a == d) {
        idx2[0] = b;
        idx2[1] = c;
        matr_elem -= get_element(source_matrix, idx2);
    }
    if (b == c) {
        idx2[0] = a;
        idx2[1] = d;
        matr_elem -= get_element(source_matrix, idx2);
    }
    if (b == d) {
        idx2[0] = a;
        idx2[1] = c;
        matr_elem += get_element(source_matrix, idx2);
    }
    return matr_elem;
}


/*
 * Sector 0h2p
 * 2-particle operator
 */
double complex slater_02_2_02(slater_det_t *bra, slater_det_t *ket)
{
    int idx4[4];
    int a = bra->indices[0];
    int b = bra->indices[1];
    int c = ket->indices[0];
    int d = ket->indices[1];
    /*
    (2) - 1.0 heff2 [ a b d c ]
    (3) + 1.0 heff2 [ a b c d ]   | => <ab||cd>
    */
    idx4[0] = a;
    idx4[1] = b;
    idx4[2] = c;
    idx4[3] = d;
    return get_element(source_matrix, idx4);
}


/*
 * Sector 2h0p
 * 1-particle operator
 */
double complex slater_20_1_20(slater_det_t *bra, slater_det_t *ket)
{
    int idx2[2];
    int i = bra->indices[0];
    int j = bra->indices[1];
    int k = ket->indices[0];
    int l = ket->indices[1];
    double complex matr_elem = 0.0 + 0.0 * I;
    /*
    (0) - 1.0 heff1 [ l j ] d_ik
    (1) + 1.0 heff1 [ k j ] d_il
    (2) + 1.0 heff1 [ l i ] d_jk
    (3) - 1.0 heff1 [ k i ] d_jl
    */
    if (i == k) {
        idx2[0] = l;
        idx2[1] = j;
        matr_elem -= get_element(source_matrix, idx2);
    }
    if (i == l) {
        idx2[0] = k;
        idx2[1] = j;
        matr_elem += get_element(source_matrix, idx2);
    }
    if (j == k) {
        idx2[0] = l;
        idx2[1] = i;
        matr_elem += get_element(source_matrix, idx2);
    }
    if (j == l) {
        idx2[0] = k;
        idx2[1] = i;
        matr_elem -= get_element(source_matrix, idx2);
    }
    return matr_elem;
}


/*
 * Sector 2h0p
 * 2-particle operator
 */
double complex slater_20_2_20(slater_det_t *bra, slater_det_t *ket)
{
    int idx4[4];
    int i = bra->indices[0];
    int j = bra->indices[1];
    int k = ket->indices[0];
    int l = ket->indices[1];
    /*
    (0) + 1.0 heff2 [ k l i j ]
    (1) - 1.0 heff2 [ l k i j ]   => <kl||ij>
    */
    idx4[0] = k;
    idx4[1] = l;
    idx4[2] = i;
    idx4[3] = j;
    return get_element(source_matrix, idx4);
}


/*
 * Sector 0h3p
 * 1-particle operator
 */
double complex slater_03_1_03(slater_det_t *bra, slater_det_t *ket)
{
    int idx2[2];
    int a = bra->indices[0];
    int b = bra->indices[1];
    int c = bra->indices[2];
    int d = ket->indices[0];
    int e = ket->indices[1];
    int f = ket->indices[2];
    double complex matr_elem = 0.0 + 0.0 * I;

    // (0) - 1.0 heff1 [ c f ] d_bd d_ae
    if ((b == d) && (a == e)) {
        idx2[0] = c;
        idx2[1] = f;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (1) + 1.0 heff1 [ c e ] d_bd d_af
    if ((b == d) && (a == f)) {
        idx2[0] = c;
        idx2[1] = e;
        matr_elem += get_element(source_matrix, idx2);
    }
    // (2) + 1.0 heff1 [ c f ] d_be d_ad
    if ((b == e) && (a == d)) {
        idx2[0] = c;
        idx2[1] = f;
        matr_elem += get_element(source_matrix, idx2);
    }
    // (3) - 1.0 heff1 [ c d ] d_be d_af
    if ((b == e) && (a == f)) {
        idx2[0] = c;
        idx2[1] = d;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (4) - 1.0 heff1 [ c e ] d_bf d_ad
    if ((b == f) && (a == d)) {
        idx2[0] = c;
        idx2[1] = e;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (5) + 1.0 heff1 [ c d ] d_bf d_ae
    if ((b == f) && (a == e)) {
        idx2[0] = c;
        idx2[1] = d;
        matr_elem += get_element(source_matrix, idx2);
    }
    // (6) + 1.0 heff1 [ b f ] d_cd d_ae
    if ((c == d) && (a == e)) {
        idx2[0] = b;
        idx2[1] = f;
        matr_elem += get_element(source_matrix, idx2);
    }
    // (7) - 1.0 heff1 [ b e ] d_cd d_af
    if ((c == d) && (a == f)) {
        idx2[0] = b;
        idx2[1] = e;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (8) - 1.0 heff1 [ a f ] d_cd d_be
    if ((c == d) && (b == e)) {
        idx2[0] = a;
        idx2[1] = f;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (9) + 1.0 heff1 [ a e ] d_cd d_bf
    if ((c == d) && (b == f)) {
        idx2[0] = a;
        idx2[1] = e;
        matr_elem += get_element(source_matrix, idx2);
    }
    // (10) - 1.0 heff1 [ b f ] d_ce d_ad
    if ((c == e) && (a == d)) {
        idx2[0] = b;
        idx2[1] = f;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (11) + 1.0 heff1 [ b d ] d_ce d_af
    if ((c == e) && (a == f)) {
        idx2[0] = b;
        idx2[1] = d;
        matr_elem += get_element(source_matrix, idx2);
    }
    // (12) + 1.0 heff1 [ a f ] d_ce d_bd
    if ((b == d) && (c == e)) {
        idx2[0] = a;
        idx2[1] = f;
        matr_elem += get_element(source_matrix, idx2);
    }
    // (13) - 1.0 heff1 [ a d ] d_ce d_bf
    if ((b == f) && (c == e)) {
        idx2[0] = a;
        idx2[1] = d;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (14) + 1.0 heff1 [ b e ] d_cf d_ad
    if ((c == f) && (a == d)) {
        idx2[0] = b;
        idx2[1] = e;
        matr_elem += get_element(source_matrix, idx2);
    }
    // (15) - 1.0 heff1 [ b d ] d_cf d_ae
    if ((c == f) && (a == e)) {
        idx2[0] = b;
        idx2[1] = d;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (16) - 1.0 heff1 [ a e ] d_cf d_bd
    if ((b == d) && (c == f)) {
        idx2[0] = a;
        idx2[1] = e;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (17) + 1.0 heff1 [ a d ] d_cf d_be
    if ((c == f) && (b == e)) {
        idx2[0] = a;
        idx2[1] = d;
        matr_elem += get_element(source_matrix, idx2);
    }

    return matr_elem;
}


/*
 * Sector 0h3p
 * 2-particle operator
 */
double complex slater_03_2_03(slater_det_t *bra, slater_det_t *ket)
{
    int idx4[4];
    int a = bra->indices[0];
    int b = bra->indices[1];
    int c = bra->indices[2];
    int d = ket->indices[0];
    int e = ket->indices[1];
    int f = ket->indices[2];
    double complex matr_elem = 0.0 + 0.0 * I;

    // sum of nine matrix elements (prefactor = 1/2)

    // d_ad <bc||ef>
    if (a == d) {
        idx4[0] = b;
        idx4[1] = c;
        idx4[2] = e;
        idx4[3] = f;
        matr_elem += get_element(source_matrix, idx4);
    }
    // d_ae <bc||fd>
    if (a == e) {
        idx4[0] = b;
        idx4[1] = c;
        idx4[2] = f;
        idx4[3] = d;
        matr_elem += get_element(source_matrix, idx4);
    }
    // d_af <bc||de>
    if (a == f) {
        idx4[0] = b;
        idx4[1] = c;
        idx4[2] = d;
        idx4[3] = e;
        matr_elem += get_element(source_matrix, idx4);
    }
    // d_bd <ac||fe>
    if (b == d) {
        idx4[0] = a;
        idx4[1] = c;
        idx4[2] = f;
        idx4[3] = e;
        matr_elem += get_element(source_matrix, idx4);
    }
    // d_be <ac||df>
    if (b == e) {
        idx4[0] = a;
        idx4[1] = c;
        idx4[2] = d;
        idx4[3] = f;
        matr_elem += get_element(source_matrix, idx4);
    }
    // d_bf <ac||ed>
    if (b == f) {
        idx4[0] = a;
        idx4[1] = c;
        idx4[2] = e;
        idx4[3] = d;
        matr_elem += get_element(source_matrix, idx4);
    }
    // d_cd <ab||ef>
    if (c == d) {
        idx4[0] = a;
        idx4[1] = b;
        idx4[2] = e;
        idx4[3] = f;
        matr_elem += get_element(source_matrix, idx4);
    }
    // d_ce <ab||fd>
    if (c == e) {
        idx4[0] = a;
        idx4[1] = b;
        idx4[2] = f;
        idx4[3] = d;
        matr_elem += get_element(source_matrix, idx4);
    }
    // d_cf <ab||de>
    if (c == f) {
        idx4[0] = a;
        idx4[1] = b;
        idx4[2] = d;
        idx4[3] = e;
        matr_elem += get_element(source_matrix, idx4);
    }

    return matr_elem;
}


/*
 * Sector 0h3p
 * 3-particle operator
 */
double complex slater_03_3_03(slater_det_t *bra, slater_det_t *ket)
{
    int idx6[6];
    int a = bra->indices[0];
    int b = bra->indices[1];
    int c = bra->indices[2];
    int d = ket->indices[0];
    int e = ket->indices[1];
    int f = ket->indices[2];
    /*
    prefactor = 1/3! = 1/6
    - 1.0 heff3 [ a b c f e d ]
    + 1.0 heff3 [ a b c e f d ]
    + 1.0 heff3 [ a b c f d e ]
    - 1.0 heff3 [ a b c d f e ]
    - 1.0 heff3 [ a b c e d f ]
    + 1.0 heff3 [ a b c d e f ]   | => +<abc||def>
    */
    idx6[0] = a;
    idx6[1] = b;
    idx6[2] = c;
    idx6[3] = d;
    idx6[4] = e;
    idx6[5] = f;

    return get_element(source_matrix, idx6);
}


/*
 * Sector 0h3p
 * 1-particle operator
 */
double complex slater_30_1_30(slater_det_t *bra, slater_det_t *ket)
{
    int idx[2];
    int i = bra->indices[0];
    int j = bra->indices[1];
    int k = bra->indices[2];
    int l = ket->indices[0];
    int m = ket->indices[1];
    int n = ket->indices[2];
    double complex matr_elem = 0.0 + 0.0 * I;

    // (0) + 1.0 heff1 [ n k ] d_jl d_im
    if ((j == l) && (i == m)) {
        idx[0] = n;
        idx[1] = k;
        matr_elem += get_element(source_matrix, idx);
    }

    // (1) - 1.0 heff1 [ m k ] d_jl d_in
    if ((j == l) && (i == n)) {
        idx[0] = m;
        idx[1] = k;
        matr_elem -= get_element(source_matrix, idx);
    }

    // (2) - 1.0 heff1 [ n k ] d_jm d_il
    if ((j == m) && (i == l)) {
        idx[0] = n;
        idx[1] = k;
        matr_elem -= get_element(source_matrix, idx);
    }

    // (3) + 1.0 heff1 [ l k ] d_jm d_in
    if ((j == m) && (i == n)) {
        idx[0] = l;
        idx[1] = k;
        matr_elem += get_element(source_matrix, idx);
    }

    // (4) + 1.0 heff1 [ m k ] d_jn d_il
    if ((j == n) && (i == l)) {
        idx[0] = m;
        idx[1] = k;
        matr_elem += get_element(source_matrix, idx);
    }

    // (5) - 1.0 heff1 [ l k ] d_jn d_im
    if ((j == n) && (i == m)) {
        idx[0] = l;
        idx[1] = k;
        matr_elem -= get_element(source_matrix, idx);
    }

    // (6) - 1.0 heff1 [ n j ] d_kl d_im
    if ((k == l) && (i == m)) {
        idx[0] = n;
        idx[1] = j;
        matr_elem -= get_element(source_matrix, idx);
    }

    // (7) + 1.0 heff1 [ m j ] d_kl d_in
    if ((k == l) && (i == n)) {
        idx[0] = m;
        idx[1] = j;
        matr_elem += get_element(source_matrix, idx);
    }

    // (8) + 1.0 heff1 [ n i ] d_kl d_jm
    if ((k == l) && (j == m)) {
        idx[0] = n;
        idx[1] = i;
        matr_elem += get_element(source_matrix, idx);
    }

    // (9) - 1.0 heff1 [ m i ] d_kl d_jn
    if ((k == l) && (j == n)) {
        idx[0] = m;
        idx[1] = i;
        matr_elem -= get_element(source_matrix, idx);
    }

    // (10) + 1.0 heff1 [ n j ] d_km d_il
    if ((k == m) && (i == l)) {
        idx[0] = n;
        idx[1] = j;
        matr_elem += get_element(source_matrix, idx);
    }

    // (11) - 1.0 heff1 [ l j ] d_km d_in
    if ((k == m) && (i == n)) {
        idx[0] = l;
        idx[1] = j;
        matr_elem -= get_element(source_matrix, idx);
    }

    // (12) - 1.0 heff1 [ n i ] d_km d_jl
    if ((k == m) && (j == l)) {
        idx[0] = n;
        idx[1] = i;
        matr_elem -= get_element(source_matrix, idx);
    }

    // (13) + 1.0 heff1 [ l i ] d_km d_jn
    if ((k == m) && (j == n)) {
        idx[0] = l;
        idx[1] = i;
        matr_elem += get_element(source_matrix, idx);
    }

    // (14) - 1.0 heff1 [ m j ] d_kn d_il
    if ((k == n) && (i == l)) {
        idx[0] = m;
        idx[1] = j;
        matr_elem -= get_element(source_matrix, idx);
    }

    // (15) + 1.0 heff1 [ l j ] d_kn d_im
    if ((k == n) && (i == m)) {
        idx[0] = l;
        idx[1] = j;
        matr_elem += get_element(source_matrix, idx);
    }

    // (16) + 1.0 heff1 [ m i ] d_kn d_jl
    if ((k == n) && (j == l)) {
        idx[0] = m;
        idx[1] = i;
        matr_elem += get_element(source_matrix, idx);
    }

    // (17) - 1.0 heff1 [ l i ] d_kn d_jm
    if ((k == n) && (j == m)) {
        idx[0] = l;
        idx[1] = i;
        matr_elem -= get_element(source_matrix, idx);
    }

    return matr_elem;
}


/*
 * Sector 0h3p
 * 2-particle operator
 */
double complex slater_30_2_30(slater_det_t *bra, slater_det_t *ket)
{
    int idx[4];
    int i = bra->indices[0];
    int j = bra->indices[1];
    int k = bra->indices[2];
    int l = ket->indices[0];
    int m = ket->indices[1];
    int n = ket->indices[2];
    double complex matr_elem = 0.0 + 0.0 * I;

    // (0) + 1.0 heff2 [ m n j k ] d_il
    if (i == l) {
        idx[0] = m;
        idx[1] = n;
        idx[2] = j;
        idx[3] = k;
        matr_elem += get_element(source_matrix, idx);
    }

    // (1) + 1.0 heff2 [ n l j k ] d_im
    if (i == m) {
        idx[0] = n;
        idx[1] = l;
        idx[2] = j;
        idx[3] = k;
        matr_elem += get_element(source_matrix, idx);
    }

    // (2) + 1.0 heff2 [ l m j k ] d_in
    if (i == n) {
        idx[0] = l;
        idx[1] = m;
        idx[2] = j;
        idx[3] = k;
        matr_elem += get_element(source_matrix, idx);
    }

    // (3) + 1.0 heff2 [ n m i k ] d_jl
    if (j == l) {
        idx[0] = n;
        idx[1] = m;
        idx[2] = i;
        idx[3] = k;
        matr_elem += get_element(source_matrix, idx);
    }

    // (4) + 1.0 heff2 [ l n i k ] d_jm
    if (j == m) {
        idx[0] = l;
        idx[1] = n;
        idx[2] = i;
        idx[3] = k;
        matr_elem += get_element(source_matrix, idx);
    }

    // (5) + 1.0 heff2 [ m l i k ] d_jn
    if (j == n) {
        idx[0] = m;
        idx[1] = l;
        idx[2] = i;
        idx[3] = k;
        matr_elem += get_element(source_matrix, idx);
    }

    // (6) + 1.0 heff2 [ m n i j ] d_kl
    if (k == l) {
        idx[0] = m;
        idx[1] = n;
        idx[2] = i;
        idx[3] = j;
        matr_elem += get_element(source_matrix, idx);
    }

    // (7) + 1.0 heff2 [ n l i j ] d_km
    if (k == m) {
        idx[0] = n;
        idx[1] = l;
        idx[2] = i;
        idx[3] = j;
        matr_elem += get_element(source_matrix, idx);
    }

    // (8) + 1.0 heff2 [ l m i j ] d_kn
    if (k == n) {
        idx[0] = l;
        idx[1] = m;
        idx[2] = i;
        idx[3] = j;
        matr_elem += get_element(source_matrix, idx);
    }

    return matr_elem;
}


/*
 * Sector 0h3p
 * 3-particle operator
 */
double complex slater_30_3_30(slater_det_t *bra, slater_det_t *ket)
{
    int idx[6];
    int i = bra->indices[0];
    int j = bra->indices[1];
    int k = bra->indices[2];
    int l = ket->indices[0];
    int m = ket->indices[1];
    int n = ket->indices[2];

    /*
     prefactor = 1/3! = 1/6
     - 6.0 heff3 [ l m n i j k ]
     + 6.0 heff3 [ l n m i j k ]
     + 6.0 heff3 [ m l n i j k ]
     - 6.0 heff3 [ m n l i j k ]
     - 6.0 heff3 [ n l m i j k ]
     + 6.0 heff3 [ n m l i j k ]
     => + <nml||ijk>
     */
    idx[0] = n;
    idx[1] = m;
    idx[2] = l;
    idx[3] = i;
    idx[4] = j;
    idx[5] = k;

    return get_element(source_matrix, idx);
}


/*
 * Sector 1h2p
 * 1-particle operator
 */
double complex slater_12_1_12(slater_det_t *bra, slater_det_t *ket)
{
    // (0) + 1.0 heff1 [ j i ] d_bc d_ad
    // (1) - 1.0 heff1 [ j i ] d_bd d_ac
    // (2) + 1.0 heff1 [ b d ] d_ij d_ac
    // (3) - 1.0 heff1 [ b c ] d_ij d_ad
    // (4) - 1.0 heff1 [ a d ] d_ij d_bc
    // (5) + 1.0 heff1 [ a c ] d_ij d_bd

    int idx2[2];
    int i = bra->indices[0];
    int a = bra->indices[1];
    int b = bra->indices[2];
    int j = ket->indices[0];
    int c = ket->indices[1];
    int d = ket->indices[2];
    double complex matr_elem = 0.0 + 0.0 * I;

    // (0) + 1.0 heff1 [ j i ] d_bc d_ad
    if ((b == c) && (a == d)) {
        idx2[0] = j;
        idx2[1] = i;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (1) - 1.0 heff1 [ j i ] d_bd d_ac
    if ((b == d) && (a == c)) {
        idx2[0] = j;
        idx2[1] = i;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (2) + 1.0 heff1 [ b d ] d_ij d_ac
    if ((i == j) && (a == c)) {
        idx2[0] = b;
        idx2[1] = d;
        matr_elem += get_element(source_matrix, idx2);
    }
    // (3) - 1.0 heff1 [ b c ] d_ij d_ad
    if ((i == j) && (a == d)) {
        idx2[0] = b;
        idx2[1] = c;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (4) - 1.0 heff1 [ a d ] d_ij d_bc
    if ((i == j) && (b == c)) {
        idx2[0] = a;
        idx2[1] = d;
        matr_elem -= get_element(source_matrix, idx2);
    }
    // (5) + 1.0 heff1 [ a c ] d_ij d_bd
    if ((i == j) && (b == d)) {
        idx2[0] = a;
        idx2[1] = c;
        matr_elem += get_element(source_matrix, idx2);
    }

    return matr_elem;
}


/*
 * Sector 1h2p
 * 2-particle operator
 */
double complex slater_12_2_12(slater_det_t *bra, slater_det_t *ket)
{
    int idx4[4];
    int i = bra->indices[0];
    int a = bra->indices[1];
    int b = bra->indices[2];
    int j = ket->indices[0];
    int c = ket->indices[1];
    int d = ket->indices[2];
    double complex matr_elem = 0.0 + 0.0 * I;

    // (0) - heff2 [ b j d i ] d_ac
    if (a == c) {
        idx4[0] = b;
        idx4[1] = j;
        idx4[2] = d;
        idx4[3] = i;
        matr_elem -= get_element(source_matrix, idx4);
    }
    // (1) + heff2 [ b j c i ] d_ad
    if (a == d) {
        idx4[0] = b;
        idx4[1] = j;
        idx4[2] = c;
        idx4[3] = i;
        matr_elem += get_element(source_matrix, idx4);
    }
    // (4) + heff2 [ a j d i ] d_bc
    if (b == c) {
        idx4[0] = a;
        idx4[1] = j;
        idx4[2] = d;
        idx4[3] = i;
        matr_elem += get_element(source_matrix, idx4);
    }
    // (6) - heff2 [ a j c i ] d_bd
    if (b == d) {
        idx4[0] = a;
        idx4[1] = j;
        idx4[2] = c;
        idx4[3] = i;
        matr_elem -= get_element(source_matrix, idx4);
    }
    // (19) + heff2 [ a b c d ] d_ij
    if (i == j) {
        idx4[0] = a;
        idx4[1] = b;
        idx4[2] = c;
        idx4[3] = d;
        matr_elem += get_element(source_matrix, idx4);
    }

    return matr_elem;
}


/*
 * Sector 1h2p
 * 3-particle operator
 */
double complex slater_12_3_12(slater_det_t *bra, slater_det_t *ket)
{
    int idx6[6];
    int i = bra->indices[0];
    int a = bra->indices[1];
    int b = bra->indices[2];
    int j = ket->indices[0];
    int c = ket->indices[1];
    int d = ket->indices[2];

    idx6[0] = a;
    idx6[1] = b;
    idx6[2] = j;
    idx6[3] = c;
    idx6[4] = d;
    idx6[5] = i;

    // - 1.0 heff3 [ a b j c d i ]
    // prefactor = 1/6 is accounted for
    return -get_element(source_matrix, idx6);
}


/*
 * Sectors 0h1p - 1h2p
 * 1-particle operator
 */
double complex slater_01_1_12(slater_det_t *bra, slater_det_t *ket)
{
    int idx[2];
    double complex matr_elem = 0.0 + 0.0 * I;

    int a = bra->indices[0];
    int i = ket->indices[0];
    int b = ket->indices[1];
    int c = ket->indices[2];

    // (0) + 1.0 heff1 [ i c ] d_ab
    if (a == b) {
        idx[0] = i;
        idx[1] = c;
        matr_elem += get_element(source_matrix, idx);
    }

    // (1) - 1.0 heff1 [ i b ] d_ac
    if (a == c) {
        idx[0] = i;
        idx[1] = b;
        matr_elem -= get_element(source_matrix, idx);
    }

    return matr_elem;
}


/*
 * Sectors 1h2p - 0h1p
 * 1-particle operator
 */
double complex slater_12_1_01(slater_det_t *bra, slater_det_t *ket)
{
    int idx[2];
    double complex matr_elem = 0.0 + 0.0 * I;

    int i = bra->indices[0];
    int a = bra->indices[1];
    int b = bra->indices[2];
    int c = ket->indices[0];

    // (0) + 1.0 heff1 [ b i ] d_ac
    if (a == c) {
        idx[0] = b;
        idx[1] = i;
        matr_elem += get_element(source_matrix, idx);
    }

    // (1) - 1.0 heff1 [ a i ] d_bc
    if (b == c) {
        idx[0] = a;
        idx[1] = i;
        matr_elem -= get_element(source_matrix, idx);
    }

    return matr_elem;
}


