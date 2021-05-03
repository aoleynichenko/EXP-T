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
 * slater_det.c
 * ============
 *
 * Tools for operating with Slater determinants.
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#include "slater_det.h"

#include <stdio.h>

#include "spinors.h"
#include "symmetry.h"

/**
 * checks if 'det' represents the vacuum determinant (no holes, no particles)
 */
int is_vacuum_det(slater_det_t *det)
{
    if (det->indices[0] == 0 && det->indices[1] == 0) {
        return 1;
    }
    else {
        return 0;
    }
}


/**
 * creates the vacuum determinant (no holes, no particles)
 */
int set_vacuum_det(slater_det_t *det)
{
    det->indices[0] = 0;
    det->indices[1] = 0;
}


/**
 * compares Slater determinants.
 * the lower the number of irrep => the "older" the determinant
 */
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


/**
 * Prints Slater determinant 'det' to the stream 'f'.
 */
void print_slater_det(FILE *f, int sect_h, int sect_p, slater_det_t *det)
{
    // здесь нужно сделать исключение -- печать детерминанта из сектора 00, который у нас считается вакуумным
    // детектировать его очень просто -- det->indices[0] == 0 && det->indices[1] == 0
    // нужна отдельная функция, говорящая нам, это вакуумный детерминант или нет
    // и еще отдельная функция, создающая вакуумный детерминант
    if (is_vacuum_det(det)) {
        fprintf(f, "| vac > (%4s)\n", get_irrep_name(det->sym));
        return;
    }

    fprintf(f, "|");
    for (int ih = 0; ih < sect_h; ih++) {
        int idx = det->indices[ih];
        int rep = spinor_info[idx].repno;
        fprintf(f, " %4s #%4d (%12.6f)", get_irrep_name(rep), idx + 1, spinor_info[idx].eps);
    }
    if (sect_h != 0 && sect_p != 0) {
        fprintf(f, " -> ");
    }
    for (int ip = 0; ip < sect_p; ip++) {
        int idx = det->indices[sect_h + ip];
        int rep = spinor_info[idx].repno;
        fprintf(f, " %4s #%4d (%12.6f)", get_irrep_name(rep), idx + 1, spinor_info[idx].eps);
    }
    fprintf(f, " > (%4s)\n", get_irrep_name(det->sym));
}
