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

#ifndef EXPT_MAPPING_H
#define EXPT_MAPPING_H

enum {
    MAPPING_IDENTITY,
    MAPPING_MESHKOV_08
};


/*
 * For details on the mapping procedure, see:
 * V. V. Meshkov, A. V. Stolyarov, R. J. Le Roy
 * Adaptive analytical mapping procedure for efficiently solving the radial Schrodinger equation
 * Phys. Rev. A 78, 052510 (2008)
 * doi: 10.1103/PhysRevA.78.052510
 */
typedef struct mapping {

    /*
     * mapping formula
     */
    int type;

    /*
     * array containing parameters of the mapping
     */
    double *params;

    /*
     * direct mapping function
     * y = y(r)
     */
    double (*y2r)(struct mapping *map, double y);

    /*
     * inverse mapping funciton
     * r = r(y)
     */
    double (*r2y)(struct mapping *map, double r);

    /*
     * factor for the transformed wavefunction (squared)
     * g(y) = (dr/dy)(y)
     */
    double (*drdy)(struct mapping *map, double y);

    /*
     * additive term in the transformed Schrodinger equation
     * F(y) = g''/(2g) - 3/4 (g'/g)^2
     */
    double (*F)(struct mapping *map, double y);

} mapping_t;


mapping_t *new_mapping(int type, double *params);

void delete_mapping(mapping_t *map);

double *construct_radial_grid(int n_points, double rmin, double rmax, mapping_t *map);

void restore_wavefunction(int n_points, double *radial_grid, mapping_t *mapping, double *phi, double *psi);

#endif //EXPT_MAPPING_H
