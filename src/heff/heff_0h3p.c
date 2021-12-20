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
#include "slater_rules.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"

// for slater rules for the 0h3p sector
extern void *source_matrix;
extern double complex (*get_element)(void *source, void *indices);


// Sector 0h3p
// 1-particle operator
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


// Sector 0h3p
// 2-particle operator
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


// Sector 0h3p
// 3-particle operator
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


// Sector 1h2p
// 1-particle operator
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


// Sector 1h2p
// 2-particle operator
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


// Sector 1h2p
// 3-particle operator
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
    return - get_element(source_matrix, idx6);
}

