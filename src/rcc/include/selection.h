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

/*******************************************************************************
* selection.h
* ===========
*
* Tools for selection of cluster amplitudes which satisfy some condition.
* All other amplitudes will be set to zero.
* These subroutines not affect performance; they are designed just for fast
* construction and testing of new approximate schemes.
*
* 2020-2021 Alexander Oleynichenko
******************************************************************************/

#ifndef CC_SELECTION_H_INCLUDED
#define CC_SELECTION_H_INCLUDED

enum {
    CC_SELECTION_SET_ZERO,
    CC_SELECTION_SET_ZERO_EXCEPT,
    CC_SELECTION_ALL,
    CC_SELECTION_SPECTATOR,
    CC_SELECTION_ACT_TO_ACT,
    CC_SELECTION_MAX_2_INACT,
    CC_SELECTION_EXC_WINDOW,
    CC_SELECTION_EPS_WINDOW
};

typedef struct ampl_selection {
    int sect_h;
    int sect_p;
    int rank;
    int task;  // ZERO or ZERO_EXCEPT
    int rule;  // ALL, SPECTATOR, ACT_TO_ACT, EXC_WINDOW or EPS_WINDOW
    double e1;
    double e2;
} ampl_selection_t;

void apply_selections(int sect_h, int sect_p, char *diag_name);

#endif /* CC_SELECTION_H_INCLUDED */
