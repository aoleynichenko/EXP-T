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

/*
 * Tools for conversion of coordinates: zmatrix -> cartesian.
 *
 * The algorithm is based on the work of Robert Shaw, see for details:
 * https://github.com/robashaw/geomConvert/blob/master/gcutil.py
 *
 * 2023 Alexander Oleynichenko
 */

#include "zmatrix.h"

#include <math.h>
#include "molecule.h"


static double vector_norm_3d(double *r);

static void cross_product_3d(double *r1, double *r2, double *prod);

static void normalize_3d(double *r);

static void subtract_vectors_3d(double *a, double *b, double *ba);


void convert_zmatrix_to_xyz(int n_atoms, int *charges, int *rconnect, double *rlist,
                            int *aconnect, double *alist, int *dconnect, double *dlist, molecule_t *mol)
{
    if (n_atoms > 0) {
        molecule_add_atom(mol, charges[0], 0.0, 0.0, 0.0);
    }

    if (n_atoms > 1) {
        molecule_add_atom(mol, charges[1], rlist[1], 0.0, 0.0);
    }

    if (n_atoms > 2) {

        // third atom in the xy-plane
        // such that the angle a012 is correct
        int iatom = 2;
        int atom_i = rconnect[iatom] - 1;
        int atom_j = aconnect[iatom] - 1;

        double r = rlist[iatom];
        double theta = alist[iatom] * M_PI / 180.0;
        double x = r * cos(theta);
        double y = r * sin(theta);

        double r1[3];
        r1[0] = mol->x[atom_i];
        r1[1] = mol->y[atom_i];
        r1[2] = mol->z[atom_i];

        double r2[3];
        r2[0] = mol->x[atom_j];
        r2[1] = mol->y[atom_j];
        r2[2] = mol->z[atom_j];

        double r21[3];
        subtract_vectors_3d(r1, r2, r21);

        if (r21[0] < 0) {
            x = r1[0] - x;
            y = r1[1] - y;
        }
        else {
            x = r1[0] + x;
            y = r1[1] + y;
        }

        molecule_add_atom(mol, charges[iatom], x, y, 0.0);
    }

    for (int iatom = 3; iatom < n_atoms; iatom++) {

        double ab[3], bc[3], nv[3], ncbc[3];

        // distance_3d, angle, dihedral angle for this atom
        double r = rlist[iatom];
        double theta = alist[iatom] * M_PI / 180.0;
        double phi = dlist[iatom] * M_PI / 180.0;

        double x = r * cos(theta);
        double y = r * cos(phi) * sin(theta);
        double z = r * sin(phi) * sin(theta);

        // XYZ coordinates of three atoms connected with this one
        int atom_i = rconnect[iatom] - 1;
        int atom_j = aconnect[iatom] - 1;
        int atom_k = dconnect[iatom] - 1;

        double a[3];
        a[0] = mol->x[atom_k];
        a[1] = mol->y[atom_k];
        a[2] = mol->z[atom_k];

        double b[3];
        b[0] = mol->x[atom_j];
        b[1] = mol->y[atom_j];
        b[2] = mol->z[atom_j];

        double c[3];
        c[0] = mol->x[atom_i];
        c[1] = mol->y[atom_i];
        c[2] = mol->z[atom_i];

        // ab = b - a
        subtract_vectors_3d(a, b, ab);

        // bc = (c - b) / ||c - b||
        subtract_vectors_3d(b, c, bc);
        normalize_3d(bc);

        // nv = ab x bc / || ab x bc ||
        cross_product_3d(ab, bc, nv);
        normalize_3d(nv);

        // ncbc = nv x bc
        cross_product_3d(nv, bc, ncbc);

        double new_x = c[0] - bc[0] * x + ncbc[0] * y + nv[0] * z;
        double new_y = c[1] - bc[1] * x + ncbc[1] * y + nv[1] * z;
        double new_z = c[2] - bc[2] * x + ncbc[2] * y + nv[2] * z;

        molecule_add_atom(mol, charges[iatom], new_x, new_y, new_z);
    }
}


/*
 * some necessary subroutines for operators with 3D vectors
 */


/*
 * ba = b - a
 */
static void subtract_vectors_3d(double *a, double *b, double *ba)
{
    ba[0] = b[0] - a[0];
    ba[1] = b[1] - a[1];
    ba[2] = b[2] - a[2];
}


static void normalize_3d(double *r)
{
    double norm = vector_norm_3d(r);
    r[0] /= norm;
    r[1] /= norm;
    r[2] /= norm;
}


static void cross_product_3d(double *r1, double *r2, double *prod)
{
    prod[0] = + r1[1] * r2[2] - r1[2] * r2[1];
    prod[1] = - r1[0] * r2[2] + r1[2] * r2[0];
    prod[2] = + r1[0] * r2[1] - r1[1] * r2[0];
}


static double vector_norm_3d(double *r)
{
    return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
}
