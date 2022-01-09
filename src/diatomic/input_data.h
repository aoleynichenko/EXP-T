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

#ifndef EXPT_INPUT_DATA_H
#define EXPT_INPUT_DATA_H

#include "cubic_spline.h"

typedef struct {
    double reduced_mass;
    int v_min;
    int v_max;
    int J_min;
    int J_max;
    int n_points;
    double *r;
    cubic_spline_t *pot1;
    cubic_spline_t *pot2;
    cubic_spline_t *prop;
    int grid_size;
    double min_energy;
} input_data_t;

void print_input_data(input_data_t *data);

void delete_input_data(input_data_t *data);

#endif //EXPT_INPUT_DATA_H
