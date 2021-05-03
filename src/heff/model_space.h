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
 * model_space.h
 * =============
 *
 * Model space in FS-CC consists of Slater determinants.
 * Current implementation deals with only the quasi-complete model spaces.
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#ifndef CC_MODEL_SPACE_H_INCLUDED
#define CC_MODEL_SPACE_H_INCLUDED

#include "comdef.h"
#include "slater_det.h"
#include "spinors.h"

/**
 * Returns total number of model-space determinants for the given FS sector
 */
size_t get_model_space_size(int sect_h, int sect_p, int nspinors, spinor_attr_t *spinor_info);

/**
 * Returns list of model-space determinants.
 */
slater_det_t **create_model_dets(int sect_h, int sect_p, size_t *ms_rep_sizes);

#endif /* CC_MODEL_SPACE_H_INCLUDED */
