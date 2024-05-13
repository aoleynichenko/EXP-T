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

#include "mapping.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "units.h"

mapping_t *new_mapping_identity();

mapping_t *new_mapping_meshkov_08(double re, double beta);


/**
 * constructs new mapping with given parameters 'params' and return
 * pointer to it.
 *
 * possible values of the 'type' argument:
 *
 * MAPPING_IDENTITY     identity mapping (= nothing will be done)
 * MAPPING_MESHKOV_08   reduced variable mapping
 */
mapping_t *new_mapping(int type, double *params)
{
    assert(type == MAPPING_IDENTITY || type == MAPPING_MESHKOV_08);

    if (type == MAPPING_IDENTITY) {
        return new_mapping_identity();
    }
    else if (type == MAPPING_MESHKOV_08) {
        return new_mapping_meshkov_08(params[0], params[1]);
    }
}


/**
 * deallocates memory associates with the 'map' object
 */
void delete_mapping(mapping_t *map)
{
    if (map->params != NULL) {
        free(map->params);
    }
    free(map);
}


/**
 * constructs radial grid for the given mapping
 */
double *construct_radial_grid(int n_points, double rmin, double rmax, mapping_t *map)
{
    double ymin = map->r2y(map, rmin);
    double ymax = map->r2y(map, rmax);

    double h = (ymax - ymin) / (n_points + 1);

    double *r = (double *) calloc(n_points, sizeof(double));

    for (int i = 0; i < n_points; i++) {
        double yi = ymin + (i + 1) * h;
        r[i] = map->y2r(map, yi);
    }

    return r;
}


/**
 * restores wavefunction psi(r) = sqrt( dr/dy(y) ) * phi(y).
 * psi[] and phi[] arrays must be pre-allocated.
 */
void restore_wavefunction(int n_points, double *radial_grid, mapping_t *mapping, double *phi, double *psi)
{
    for (int i = 0; i < n_points; i++) {
        double ri = radial_grid[i];
        double yi = mapping->r2y(mapping, ri);
        double gi = mapping->drdy(mapping, yi);

        psi[i] = sqrt(gi) * phi[i];
    }
}


/**
 * "Identity" mapping:
 * y = r
 */
mapping_t *new_mapping_identity()
{
    double y2r_identity(mapping_t *map, double y);
    double r2y_identity(mapping_t *map, double r);
    double drdy_identity(mapping_t *map, double y);
    double F_identity(mapping_t *map, double y);

    mapping_t *mapping = (mapping_t *) malloc(sizeof(mapping_t));

    mapping->type = MAPPING_IDENTITY;
    mapping->params = NULL;
    mapping->y2r = y2r_identity;
    mapping->r2y = r2y_identity;
    mapping->drdy = drdy_identity;
    mapping->F = F_identity;

    return mapping;
}


double y2r_identity(mapping_t *map, double y)
{
    return y;
}


double r2y_identity(mapping_t *map, double r)
{
    return r;
}


double drdy_identity(mapping_t *map, double y)
{
    return 1.0;
}


double F_identity(mapping_t *map, double y)
{
    return 0.0;
}


/**
 * mapping procedure proposed by V. V. Meshkov and co-authors
 * for calculations of highly excited rovibrational states near the dissociation limit
 *
 * V. V. Meshkov, A. V. Stolyarov, R. J. Le Roy
 * Adaptive analytical mapping procedure for efficiently solving the radial Schrodinger equation
 * Phys. Rev. A 78, 052510 (2008)
 * doi: 10.1103/PhysRevA.78.052510
 */
mapping_t *new_mapping_meshkov_08(double re, double beta)
{
    double y2r_meshkov_08(mapping_t *map, double y);
    double r2y_meshkov_08(mapping_t *map, double r);
    double drdy_meshkov_08(mapping_t *map, double y);
    double F_meshkov_08(mapping_t *map, double y);

    mapping_t *mapping = (mapping_t *) malloc(sizeof(mapping_t));

    mapping->type = MAPPING_MESHKOV_08;
    mapping->params = (double *) calloc(2, sizeof(double));
    mapping->params[0] = re;
    mapping->params[1] = beta;

    mapping->y2r = y2r_meshkov_08;
    mapping->r2y = r2y_meshkov_08;
    mapping->drdy = drdy_meshkov_08;
    mapping->F = F_meshkov_08;

    return mapping;
}


double y2r_meshkov_08(mapping_t *map, double y)
{
    double re = map->params[0];
    double beta = map->params[1];

    return re * pow((1.0 + y)/(1.0 - y), 1.0/beta);
}


double r2y_meshkov_08(mapping_t *map, double r)
{
    double re = map->params[0];
    double beta = map->params[1];

    double rbeta = pow(r, beta);

    return 2.0 * rbeta / (rbeta + pow(re,beta)) - 1.0;
}


double drdy_meshkov_08(mapping_t *map, double y)
{
    double re = map->params[0];
    double beta = map->params[1];

    double numer = pow(1.0 + y, 1.0/beta - 1.0);
    double denom = pow(1.0 - y, 1.0/beta + 1.0);

    return (2 * re / beta) * numer / denom;
}


double F_meshkov_08(mapping_t *map, double y)
{
    double re = map->params[0];
    double beta = map->params[1];

    return (1.0 - 1.0 / (beta * beta)) / pow(1.0 - y * y, 2);
}

