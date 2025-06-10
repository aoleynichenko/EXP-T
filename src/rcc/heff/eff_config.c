/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2025 The EXP-T developers.
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

#include "eff_config.h"

#include <complex.h>
#include <string.h>

#include "model_space.h"
#include "slater_det.h"


/**
 * calculates effective occupation numbers of active spinors
 * for the given electronic state (represented by a model vector)
 */
void get_eff_configuration(
        int sect_h,
        int sect_p,
        int ms_size,
        slater_det_t *det_list,
        double complex *model_vector_coefs,
        double *config
)
{
    /*
     * construct list of indices of active spinors
     */
    int n_active = 0;
    int active_spinors[CC_MAX_SPINORS]; // local -> global spinor index mapping
    get_active_space(sect_h, sect_p, &n_active, active_spinors);

    memset(config, 0, n_active * sizeof(double));

    /*
     * loop over model determinants
     */
    for (int i = 0; i < ms_size; i++) {
        double det_coefft = model_vector_coefs[i];
        double det_weight = cabs(det_coefft) * cabs(det_coefft);
        slater_det_t *det = det_list + i;

        /*
         * for each determinant: determine which active spinors are occupied
         */
        for (int j = 0; j < sect_h + sect_p; j++) {
            int act_spinor_index = det->indices[j];

            for (int k = 0; k < n_active; k++) {
                if (active_spinors[k] == act_spinor_index) {
                    config[k] += det_weight;
                    break;
                }
            }
        }
    }
}
