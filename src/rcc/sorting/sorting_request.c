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

#include <string.h>

#include "sorting_request.h"

sorting_request_t sorting_requests[CC_MAX_SORTING_REQUESTS];
int n_requests = 0;


sorting_request_t *append_sorting_request(sorting_request_t *requests, int *num_requests,
                                          char *name, char *qparts, char *valence, char *order)
{
    sorting_request_t *req;

    req = sorting_requests + (*num_requests);
    req->dg = diagram_stack_find(name);
    strcpy(req->dg_name, name);
    strcpy(req->hp, qparts);
    strcpy(req->valence, valence);
    strcpy(req->order, order);
    req->done = 0;

    *num_requests = *num_requests + 1;
    return req;
}

