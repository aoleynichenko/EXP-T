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
 * slater_rules.c
 * ==============
 *
 * Slater rules: evaluation of matrix elements in the basis of Slater
 * determinants.
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "slater_rules.h"

void *source_matrix;

double complex (*get_element)(void *source, void *indices);

double complex (*slater)(slater_det_t *d1, slater_det_t *d2) = NULL;

double complex slater_01_1_01(slater_det_t *bra, slater_det_t *ket);

double complex slater_10_1_10(slater_det_t *bra, slater_det_t *ket);

double complex slater_02_1_02(slater_det_t *bra, slater_det_t *ket);
double complex slater_02_2_02(slater_det_t *bra, slater_det_t *ket);

double complex slater_20_1_20(slater_det_t *bra, slater_det_t *ket);
double complex slater_20_2_20(slater_det_t *bra, slater_det_t *ket);

double complex slater_11_1_11(slater_det_t *bra, slater_det_t *ket);
double complex slater_11_2_11(slater_det_t *bra, slater_det_t *ket);

double complex slater_03_1_03(slater_det_t *bra, slater_det_t *ket);
double complex slater_03_2_03(slater_det_t *bra, slater_det_t *ket);
double complex slater_03_3_03(slater_det_t *bra, slater_det_t *ket);

double complex slater_12_1_12(slater_det_t *bra, slater_det_t *ket);
double complex slater_12_2_12(slater_det_t *bra, slater_det_t *ket);
double complex slater_12_3_12(slater_det_t *bra, slater_det_t *ket);


/*******************************************************************************
 * setup_slater
 *
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
 ******************************************************************************/
void setup_slater(void *source, double complex (*getter)(void *source, void *indices),
                  int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p, int npart)
{
    source_matrix = source;
    get_element = getter;

    if (bra_sect_h == 0 && bra_sect_p == 1 && ket_sect_h == 0 && ket_sect_p == 1) {
        slater = slater_01_1_01;
    }
    else if (bra_sect_h == 1 && bra_sect_p == 0 && ket_sect_h == 1 && ket_sect_p == 0) {
        slater = slater_10_1_10;
    }
    else if (bra_sect_h == 1 && bra_sect_p == 1 && ket_sect_h == 1 && ket_sect_p == 1) {
        if (npart == 1) {
            slater = slater_11_1_11;
        }
        else if (npart == 2) {
            slater = slater_11_2_11;
        }
    }
    else if (bra_sect_h == 0 && bra_sect_p == 2 && ket_sect_h == 0 && ket_sect_p == 2) {
        if (npart == 1) {
            slater = slater_02_1_02;
        }
        else if (npart == 2) {
            slater = slater_02_2_02;
        }
    }
    else if (bra_sect_h == 2 && bra_sect_p == 0 && ket_sect_h == 2 && ket_sect_p == 0) {
        if (npart == 1) {
            slater = slater_20_1_20;
        }
        else if (npart == 2) {
            slater = slater_20_2_20;
        }
    }
#ifdef VERSION_DEVEL
    else if (bra_sect_h == 0 && bra_sect_p == 3 && ket_sect_h == 0 && ket_sect_p == 3) {
        if (npart == 1) {
            slater = slater_03_1_03;
        }
        else if (npart == 2) {
            slater = slater_03_2_03;
        }
        else if (npart == 3) {
            slater = slater_03_3_03;
        }
    }
    else if (bra_sect_h == 1 && bra_sect_p == 2 && ket_sect_h == 1 && ket_sect_p == 2) {
        if (npart == 1) {
            slater = slater_12_1_12;
        }
        else if (npart == 2) {
            slater = slater_12_2_12;
        }
        else if (npart == 3) {
            slater = slater_12_3_12;
        }
    }
#endif
    else {
        printf("Cannot evaluate matrix element between Slater determinants of type |%dh%dp> and |%dh%dp>\n",
               bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p);
        printf("see %s for details, line %d\n", __FILE__, __LINE__);
        exit(0);
    }
}


/*******************************************************************************
 * Slater rules for different sector and n-particle operators
 ******************************************************************************/

// Sector 0h1p
// 1-particle operator
double complex slater_01_1_01(slater_det_t *d1, slater_det_t *d2)
{
    int idx[2];
    idx[0] = d1->indices[0];
    idx[1] = d2->indices[0];
    return get_element(source_matrix, idx);
}


// Sector 1h0p
// 1-particle operator
double complex slater_10_1_10(slater_det_t *d1, slater_det_t *d2)
{
    int idx[2];
    idx[1] = d1->indices[0];
    idx[0] = d2->indices[0];
    return -get_element(source_matrix, idx);
}


// Sector 1h1p
// 1-particle operator
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


// Sector 1h1p
// 2-particle operator
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


// Sector 0h2p
// 1-particle operator
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


// Sector 0h2p
// 2-particle operator
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


// Sector 2h0p
// 1-particle operator
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


// Sector 2h0p
// 2-particle operator
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
