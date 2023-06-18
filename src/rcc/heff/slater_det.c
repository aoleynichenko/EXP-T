/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2023 The EXP-T developers.
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
 * Tools for operating with Slater determinants.
 */

#include "slater_det.h"

#include <stdio.h>
#include <string.h>

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
 * calculate zero-order energy of the model space determinant
 */
double model_det_energy(int sect_h, int sect_p, slater_det_t *det)
{
    double energy = 0.0;

    // holes
    for (int h = 0; h < sect_h; h++) {
        int idx = det->indices[h];
        double eps = get_eps(idx);
        energy -= eps;
    }

    // particles
    for (int p = 0; p < sect_p; p++) {
        int idx = det->indices[sect_h + p];
        double eps = get_eps(idx);
        energy += eps;
    }

    return energy;
}


/**
 * Prints Slater determinant 'det' to the buffer 'str'
 */
void str_print_slater_det(char *str, int sect_h, int sect_p, slater_det_t *det)
{
    char buf[1024];
    str[0] = '\0';

    // здесь нужно сделать исключение -- печать детерминанта из сектора 00, который у нас считается вакуумным
    // детектировать его очень просто -- det->indices[0] == 0 && det->indices[1] == 0
    // нужна отдельная функция, говорящая нам, это вакуумный детерминант или нет
    // и еще отдельная функция, создающая вакуумный детерминант
    if (is_vacuum_det(det)) {
        sprintf(str, "| vac > (%4s)\n", get_irrep_name(det->sym));
        return;
    }

    sprintf(buf, "|");
    strcat(str, buf);

    /*
     * hole spinors
     */
    for (int ih = 0; ih < sect_h; ih++) {
        int idx = det->indices[ih];
        int rep = spinor_info[idx].repno;

        if (cc_opts->spinor_labels[idx] != NULL) {
            sprintf(buf, " %-6s #%4d %-14s", get_irrep_name(rep), idx + 1, cc_opts->spinor_labels[idx]);
        }
        else {
            sprintf(buf, " %-6s #%4d (%12.6f)", get_irrep_name(rep), idx + 1, spinor_info[idx].eps);
        }

        strcat(str, buf);
    }

    if (sect_h != 0 && sect_p != 0) {
        sprintf(buf, " -> ");
        strcat(str, buf);
    }

    /*
     * particle spinors
     */
    for (int ip = 0; ip < sect_p; ip++) {
        int idx = det->indices[sect_h + ip];
        int rep = spinor_info[idx].repno;

        if (cc_opts->spinor_labels[idx] != NULL) {
            sprintf(buf, " %-6s #%4d %-14s", get_irrep_name(rep), idx + 1, cc_opts->spinor_labels[idx]);
        }
        else {
            sprintf(buf, " %-6s #%4d (%12.6f)", get_irrep_name(rep), idx + 1, spinor_info[idx].eps);
        }

        strcat(str, buf);
    }

    /*
     * symmetry of the determinant
     */
    sprintf(buf, " > (%4s)", get_irrep_name(det->sym));
    strcat(str, buf);
}


/**
 * Prints Slater determinant 'det' to the stream 'f'.
 */
void print_slater_det(FILE *f, int sect_h, int sect_p, slater_det_t *det)
{
    int intham_imms_is_main_space_det(int sect_h, int sect_p, slater_det_t *det);
    char buf[1024];

    /*
     * print determinant
     */
    str_print_slater_det(buf, sect_h, sect_p, det);
    fprintf(f, "%s", buf);

    /*
     * main or intermediate determinant?
     */
    if (intham_imms_is_main_space_det(sect_h, sect_p, det)) {
        fprintf(f, "  main");
    }

    fprintf(f, "\n");
}


/**
 * Evaluates overlap integral between two Slater determinants,
 * < det_1 | det_2 >
 */
double kronecker_delta(int sect_h, int sect_p, slater_det_t *det_1, slater_det_t *det_2)
{
    for (int i = 0; i < sect_h + sect_p; i++) {
        if (det_1->indices[i] != det_2->indices[i]) {
            return 0.0;
        }
    }

    return 1.0;
}
