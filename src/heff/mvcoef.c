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
 * mvcoef.c
 *
 * Model vectors are stored on disk in binary files MVCOEF**.
 * This file provides interface for operating with the MVCOEF file format.
 *
 * 2020-2021 Alexander Oleynichenko
 ******************************************************************************/

#include "mvcoef.h"

#include <string.h>
#include <stdio.h>

#include "io.h"
#include "linalg.h"
#include "slater.h"
#include "symmetry.h"
#include "spinors.h"
#include "options.h"


/**
 * Creates MVCOEF* unformatted file containing model vectors coefficients
 * and related meta information.
 * @return file descriptor.
 */
int mvcoef_open(int sect_h, int sect_p)
{
    char mv_file_name[256];

    sprintf(mv_file_name, "MVCOEF%d%d", sect_h, sect_p);

    int f_mvcoef = io_open(mv_file_name, "w");

    return f_mvcoef;
}


void mvcoef_close(int file_descr, double eigval_0)
{
    char *eof_str = "EOF";
    size_t eof_len = strlen(eof_str) + 1;
    io_write(file_descr, &eof_len, sizeof(eof_len));
    io_write(file_descr, eof_str, eof_len);
    io_write(file_descr, &eigval_0, sizeof(eigval_0));
    io_close(file_descr);
}


void mvcoef_write_vectors_unformatted(int file_descr, char *rep_name,
                                      size_t nroots, size_t dim, slater_det_t *det_list,
                                      double complex *ev, double complex *vl, double complex *vr)
{
    size_t rep_name_len = strlen(rep_name) + 1;
    size_t nb = sizeof(double complex) * nroots * dim;

    io_write(file_descr, &rep_name_len, sizeof(rep_name_len));
    io_write(file_descr, rep_name, rep_name_len);

    io_write(file_descr, &dim, sizeof(dim));
    io_write(file_descr, &nroots, sizeof(nroots));
    io_write(file_descr, det_list, sizeof(slater_det_t) * dim);

    io_write(file_descr, ev, sizeof(double complex) * nroots);
    io_write(file_descr, vr, nb);
    io_write(file_descr, vl, nb);
}


/**
 * Reads model vectors for the given (sect_h,sect_p) sector from the MVCOEF-**
 * unformatted file. All data will be collected as the array of struct_mvcoef type
 * entries.
 *
 * You can set the name of the file (char *file_name). If file_name == NULL,
 * the name of the file for the (h,p) FS sector will be MVCOEF<h><p>, example: MVCOEF01
 *
 * @note: this function allocates arrays! they must be deallocated!
 */
void read_model_vectors_unformatted(int sect_h, int sect_p, char *file_name, int *nrep, struct mv_block *mv_blocks)
{
    char mvcoef_file_name[64];
    double eigval_0;
    double const AU2CM = 219474.6313702;

    *nrep = 0;

    if (file_name == NULL) {
        sprintf(mvcoef_file_name, "MVCOEF%d%d", sect_h, sect_p);
    }
    else {
        sprintf(mvcoef_file_name, "%s", file_name);
    }

    // read model vectors from unformatted file
    int f_mvcoef = io_open(mvcoef_file_name, "r");
    while (1) {
        struct mv_block *curr_block = &mv_blocks[*nrep];
        size_t rep_name_len;
        size_t ms_size, nroots;

        io_read(f_mvcoef, &rep_name_len, sizeof(rep_name_len));
        io_read(f_mvcoef, curr_block->rep_name, rep_name_len);
        if (strcmp(curr_block->rep_name, "EOF") == 0) {
            break;
        }

        io_read(f_mvcoef, &ms_size, sizeof(ms_size));
        io_read(f_mvcoef, &nroots, sizeof(nroots));
        curr_block->ms_size = ms_size;
        curr_block->nroots = nroots;

        curr_block->dets = (slater_det_t *) cc_malloc(sizeof(slater_det_t) * ms_size);
        curr_block->eigval = (double complex *) cc_malloc(sizeof(double complex) * nroots);
        curr_block->energy_cm = (double *) cc_malloc(sizeof(double) * nroots);
        curr_block->vl = (double complex *) cc_malloc(sizeof(double complex) * ms_size * nroots);
        curr_block->vr = (double complex *) cc_malloc(sizeof(double complex) * ms_size * nroots);

        io_read(f_mvcoef, curr_block->dets, sizeof(slater_det_t) * ms_size);
        io_read(f_mvcoef, curr_block->eigval, sizeof(double complex) * nroots);
        size_t nb = sizeof(double complex) * nroots * ms_size;
        io_read(f_mvcoef, curr_block->vr, nb);
        io_read(f_mvcoef, curr_block->vl, nb);

        *nrep = *nrep + 1;
    }
    io_read(f_mvcoef, &eigval_0, sizeof(eigval_0));
    io_close(f_mvcoef);

    // calculate energy of each state wrt ground state (eigval_0)
    for (size_t irep = 0; irep < *nrep; irep++) {
        struct mv_block *b = &mv_blocks[irep];
        for (size_t i = 0; i < b->nroots; i++) {
            b->energy_cm[i] = (creal(b->eigval[i]) - eigval_0) * AU2CM;
        }
    }
}
