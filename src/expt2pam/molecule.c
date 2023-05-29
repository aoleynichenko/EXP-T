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

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "molecule.h"


/**
 * constructor
 */
molecule_t *molecule_new()
{
    molecule_t *mol = (molecule_t *) malloc(sizeof(molecule_t));

    mol->n_atoms = 0;
    mol->sym_group.group = SYMMETRY_AUTO;
    mol->sym_group.xyz[0] = 0;
    mol->sym_group.xyz[1] = 0;
    mol->sym_group.xyz[2] = 0;

    mol->cage.total_charge = 0.0;
    mol->cage.n_points = 0;
    mol->cage.points = NULL;

    mol->n_point_charges;

    return mol;
}


/**
 * destructor
 */
void molecule_delete(molecule_t *mol)
{
    free(mol);
}


/**
 * adds one atom to the 'molecule_t' object
 */
void molecule_add_atom(molecule_t *mol, int nuc_charge, double x, double y, double z)
{
    int n_at = mol->n_atoms;

    assert(nuc_charge >= 0 && nuc_charge < N_CHEM_ELEMENTS);

    mol->charges[n_at] = nuc_charge;
    mol->x[n_at] = x;
    mol->y[n_at] = y;
    mol->z[n_at] = z;
    mol->n_atoms++;
}


/**
 * adds one point charge to the 'molecule_t' object
 */
void molecule_add_point_charge(molecule_t *mol, double point_charge, double x, double y, double z)
{
    int n_q = mol->n_point_charges;

    mol->point_charges[n_q] = point_charge;
    mol->qx[n_q] = x;
    mol->qy[n_q] = y;
    mol->qz[n_q] = z;

    mol->n_point_charges++;
}


/**
 * number of types of elements
 */
int molecule_n_atom_types(molecule_t *mol)
{
    int n_atoms[N_CHEM_ELEMENTS];
    int n_atom_types;

    memset(n_atoms, 0, sizeof(n_atoms));

    for (int i = 0; i < mol->n_atoms; i++) {
        n_atoms[mol->charges[i]]++;
    }

    n_atom_types = 0;
    for (int i = 0; i < N_CHEM_ELEMENTS; i++) {
        if (n_atoms[i] > 0) {
            n_atom_types++;
        }
    }

    if (mol->cage.n_points > 0) {
        n_atom_types++;
    }

    return n_atom_types;
}


/**
 * number of types of point charges
 */
int molecule_n_point_charge_types(molecule_t *mol)
{
    double point_charge_types[1000];
    int n_point_charge_types = 0;

    memset(point_charge_types, 0, sizeof(point_charge_types));

    for (int i = 0; i < mol->n_point_charges; i++) {
        double q = mol->point_charges[i];

        int charge_found = 0;
        for (int j = 0; j < n_point_charge_types; j++) {
            if (fabs(point_charge_types[j] - q) < 1e-6) {
                charge_found = 1;
                break;
            }
        }

        // new type of point charge
        if (!charge_found) {
            point_charge_types[n_point_charge_types] = q;
            n_point_charge_types++;
        }
    }

    return n_point_charge_types;
}


/**
 * counts number of atoms of the element given
 */
int molecule_n_atoms_of(molecule_t *mol, int element)
{
    int n_atoms = 0;

    for (int i = 0; i < mol->n_atoms; i++) {
        if (mol->charges[i] == element) {
            n_atoms++;
        }
    }

    return n_atoms;
}


void molecule_print(molecule_t *mol)
{
    printf("\ngeometry:\n");
    printf("---------\n");
    for (int i = 0; i < mol->n_atoms; i++) {
        printf("  %3d%12.8f%12.8f%12.8f\n", mol->charges[i], mol->x[i], mol->y[i], mol->z[i]);
    }
    for (int i = 0; i < mol->n_point_charges; i++) {
        printf("  q=%12.8f   %12.8f%12.8f%12.8f\n", mol->point_charges[i], mol->qx[i], mol->qy[i], mol->qz[i]);
    }
    printf("Symmetry: %d\n (orientation %d %d %d)\n", mol->sym_group.group,
           mol->sym_group.xyz[0], mol->sym_group.xyz[1], mol->sym_group.xyz[2]);
    if (mol->cage.n_points == 0) {
        printf("Cage: no cage\n");
    }
    else {
        printf("Cage:\n");
        printf("  total charge = %f\n", mol->cage.total_charge);
        printf("  npoints = %d\n", mol->cage.n_points);
        for (int i = 0; i < mol->cage.n_points; i++) {
            printf("  %12.8f%12.8f%12.8f\n", mol->cage.points[3 * i], mol->cage.points[3 * i + 1],
                   mol->cage.points[3 * i + 2]);
        }
    }
    printf("\n");
}
