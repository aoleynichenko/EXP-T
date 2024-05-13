/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2024 The EXP-T developers.
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

#ifndef CC_CROP_H_INCLUDED
#define CC_CROP_H_INCLUDED

#include "engine.h"

#define CROP_MAX 1000

/*
 * all the amplitudes to be used in extrapolation with corresponding error
 * vectors are stored in the stack-like structure crop_queue_t
 *
 * t1, t2, t3 -- names of diagrams with amplitudes
 * e1, e2, e3 -- error vectors
 * n -- queue length
 * do_t1 -- extrapolate singles amplitudes (T1)
 * do_t2 -- extrapolate doubles amplitudes (T2)
 * do_t3 -- extrapolate triples amplitudes (T3)
 */
typedef struct {
    char t1[CROP_MAX][CC_DIAGRAM_MAX_NAME];
    char e1[CROP_MAX][CC_DIAGRAM_MAX_NAME];
    char t2[CROP_MAX][CC_DIAGRAM_MAX_NAME];
    char e2[CROP_MAX][CC_DIAGRAM_MAX_NAME];
    char t3[CROP_MAX][CC_DIAGRAM_MAX_NAME];
    char e3[CROP_MAX][CC_DIAGRAM_MAX_NAME];
    int n;
    int do_t1;
    int do_t2;
    int do_t3;
} crop_queue_t;

crop_queue_t *new_crop_queue(int do_t1, int do_t2, int do_t3);


void delete_crop_queue(crop_queue_t *q);

void crop_push(crop_queue_t *q,
               char *diag_t1new, char *diag_r1new,
               char *diag_t2new, char *diag_r2new,
               char *diag_t3new, char *diag_r3new);

void crop_put(crop_queue_t *q,
              char *diag_t1new, char *diag_t1old,
              char *diag_t2new, char *diag_t2old,
              char *diag_t3new, char *diag_t3old,
              int iter);

void crop_truncate(crop_queue_t *q, int len);
void crop_pop(crop_queue_t *q);

void crop_extrapolate(crop_queue_t *q, char *extrap_t1, char *extrap_t2, char *extrap_t3);

#endif // CC_CROP_H_INCLUDED
