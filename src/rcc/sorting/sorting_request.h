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

#ifndef CC_SORTING_REQUEST_H_INCLUDED
#define CC_SORTING_REQUEST_H_INCLUDED

#include "engine.h"

#define CC_MAX_SORTING_REQUESTS 64

typedef struct {
    diagram_t *dg;
    char dg_name[256];
    char hp[CC_DIAGRAM_MAX_RANK];
    char valence[CC_DIAGRAM_MAX_RANK];
    char order[CC_DIAGRAM_MAX_RANK];
    int done;
} sorting_request_t;

extern sorting_request_t sorting_requests[CC_MAX_SORTING_REQUESTS];
extern int n_requests;

int num_twoelec_requests(sorting_request_t *requests, int num_requests);

sorting_request_t *append_sorting_request(sorting_request_t *requests, int *num_requests,
                                          char *name, char *qparts, char *valence, char *order);

#endif /* CC_SORTING_REQUEST_H_INCLUDED */
