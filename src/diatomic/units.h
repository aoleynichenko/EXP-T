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

#ifndef EXPT_UNITS_H
#define EXPT_UNITS_H

enum {
    UNITS_ANGSTROM,
    UNITS_ATOMIC,
    UNITS_CM
};

/*
 * identity
 */
#define ATOMIC_TO_ATOMIC 1.0

/*
 * units of length
 */
#define ATOMIC_TO_ANGSTROM 0.529177210903
#define ANGSTROM_TO_ATOMIC (1.0 / (ATOMIC_TO_ANGSTROM))

/*
 * units of energy
 */
#define ATOMIC_TO_CM 219474.63136320
#define CM_TO_ATOMIC (1.0 / (ATOMIC_TO_CM))

/*
 * units of mass
 */
#define AMU_TO_ELECTRON_MASS 1822.888486209
#define ELECTRON_MASS_TO_AMU (1.0 / (AMU_TO_ELECTRON_MASS))

/*
 * speed of light in atomic units
 */
#define SPEED_OF_LIGHT_AU 137.0359998

#endif //EXPT_UNITS_H
