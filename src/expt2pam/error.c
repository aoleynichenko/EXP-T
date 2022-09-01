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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "error.h"

#define MAX_ERR_LEN 1024


// writes error message & aborts execution
void errquit(char *fmt, ...)
{
    va_list args;
    char errbuf[MAX_ERR_LEN];

    va_start(args, fmt);
    vsprintf(errbuf, fmt, args);
    errbuf[MAX_ERR_LEN - 1] = '\0';
    va_end(args);
    printf("\nERROR:\n");
    printf("%s\n", errbuf);

    // abort execution
    exit(1);
}
