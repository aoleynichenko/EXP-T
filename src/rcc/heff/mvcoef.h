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

#ifndef CC_MVCOEF_H_INCLUDED
#define CC_MVCOEF_H_INCLUDED

#include "slater_rules.h"

/*
  data structure which is used to store data about model vectors read from
  unformatted files (see also code in heff.c)
 */
typedef struct mv_block {
    char rep_name[64];
    size_t ms_size;
    size_t nroots;
    double complex *eigval;
    double *energy_cm;
    slater_det_t *dets;
    double complex *vl;
    double complex *vr;
} mv_block_t;

void mvblock_free(struct mv_block *block);

int mvcoef_open(int sect_h, int sect_p);

void mvcoef_close(int file_descr, double lowest_root);

void mvcoef_write_vectors_unformatted(int file_descr, char *rep_name,
                                      size_t nroots, size_t dim, slater_det_t *det_list,
                                      double complex *ev, double complex *vl, double complex *vr);

void mvcoef_read_vectors_unformatted(int sect_h, int sect_p, char *file_name, int *nrep, struct mv_block *mv_blocks);

void mvcoef_read_vectors_unformatted_for_state(
        int sect_h, int sect_p, char *file_name, char *irrep_name, int state,
        int *ms_dim, slater_det_t **det_list, double *eigenvalue, double *exc_energy_cm,
        double complex **coef_left, double complex **coef_right
);

#endif /* CC_MVCOEF_H_INCLUDED */
