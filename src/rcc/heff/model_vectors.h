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

#ifndef CC_MODEL_VECTORS_H_INCLUDED
#define CC_MODEL_VECTORS_H_INCLUDED

#include "slater_det.h"

void print_model_vectors_stdout(
        int sector_h,
        int sector_p,
        slater_det_t **det_list,
        size_t *block_dims,
        double complex **eigvalues,
        double complex **coef_left,
        double complex **coef_right
);

void write_model_vectors_unformatted(
        int sector_h,
        int sector_p,
        slater_det_t **det_list,
        size_t *block_dims,
        double complex **eigvalues,
        double complex **coef_left,
        double complex **coef_right
);

void write_formatted_heff(
        int sect_h,
        int sect_p,
        size_t *block_dims,
        double complex **heff
);

void print_model_vector(
        FILE *f_out,
        int sect_h,
        int sect_p,
        int rep_no,
        char *rep_name,
        int state_no,
        double complex energy,
        size_t len,
        double complex *coeffs,
        slater_det_t *det_list,
        double coef_thresh
);

#endif // CC_MODEL_VECTORS_H_INCLUDED
