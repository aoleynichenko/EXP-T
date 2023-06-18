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
 * Model vectors are stored on disk in binary files MVCOEF**.
 * This file provides interface for operating with the MVCOEF file format.
 */

#include "mvcoef.h"

#include <string.h>
#include <stdio.h>

#include "io.h"
#include "linalg.h"
#include "slater_rules.h"
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


void mvcoef_close(int file_descr, double lowest_root)
{
    char *eof_str = "EOF";
    size_t eof_len = strlen(eof_str) + 1;
    io_write(file_descr, &eof_len, sizeof(eof_len));
    io_write(file_descr, eof_str, eof_len);
    io_write(file_descr, &lowest_root, sizeof(lowest_root));
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
void mvcoef_read_vectors_unformatted(int sect_h, int sect_p, char *file_name, int *nrep, struct mv_block *mv_blocks)
{
    char mvcoef_file_name[CC_MAX_NUM_IRREPS];
    double eigval_0;
    double const AU2CM = 219474.6313702;

    /*
     * special case of the 0h0p sector
     */
    if (sect_h == 0 && sect_p == 0) {
        *nrep = 1;

        mv_blocks[0].eigval = z_zeros(1, 1);
        mv_blocks[0].energy_cm = d_zeros(1, 1);
        mv_blocks[0].nroots = 1;
        mv_blocks[0].ms_size = 1;

        int vac_irrep = get_vacuum_irrep();
        char *vac_irrep_name = get_irrep_name(vac_irrep);
        strcpy(mv_blocks[0].rep_name, vac_irrep_name);

        mv_blocks[0].dets = (slater_det_t *) cc_calloc(1, sizeof(slater_det_t));
        mv_blocks[0].dets[0].sym = vac_irrep;

        mv_blocks[0].vl = z_zeros(1, 1);
        mv_blocks[0].vr = z_zeros(1, 1);
        mv_blocks[0].vl[0] = 1.0;
        mv_blocks[0].vr[0] = 1.0;

        return;
    }

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


void mvcoef_read_vectors_unformatted_for_state(
        int sect_h, int sect_p, char *file_name, char *irrep_name, int state,
        int *ms_dim, slater_det_t **det_list, double *eigenvalue, double *exc_energy_cm,
        double complex **coef_left, double complex **coef_right
)
{
    /*
     * extract model vectors and eigenvalues from the MVCOEF* unformatted file
     */
    int nrep;
    struct mv_block mv_blocks[CC_MAX_NUM_IRREPS];
    mvcoef_read_vectors_unformatted(sect_h, sect_p, file_name, &nrep, mv_blocks);

    /*
     * get model vectors for the required irrep
     */
    struct mv_block *mvb1 = NULL;
    for (size_t iblock = 0; iblock < nrep; iblock++) {
        if (strcmp(mv_blocks[iblock].rep_name, irrep_name) == 0) {
            mvb1 = mv_blocks + iblock;
        }
    }

    /*
     * return data via working arrays
     */
    size_t ms_size = mvb1->ms_size;
    *ms_dim = (int) ms_size;
    *det_list = (slater_det_t *) cc_memdup(mvb1->dets, sizeof(slater_det_t) * ms_size);
    *eigenvalue = creal(mvb1->eigval[state]);
    *exc_energy_cm = mvb1->energy_cm[state];
    *coef_left = (double complex *) cc_memdup(mvb1->vl + ms_size * state, sizeof(double complex) * ms_size);
    *coef_right = (double complex *) cc_memdup(mvb1->vr + ms_size * state, sizeof(double complex) * ms_size);

    /*
     * clean up
     */
    for (size_t irep = 0; irep < nrep; irep++) {
        mvblock_free(mv_blocks + irep);
    }
}


/**
 * Deallocates memory used by the mvblock_t structure.
 */
void mvblock_free(struct mv_block *block)
{
    cc_free(block->dets);
    cc_free(block->eigval);
    cc_free(block->energy_cm);
    cc_free(block->vl);
    cc_free(block->vr);
}
