/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2021 The EXP-T developers.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basis.h"
#include "ecp.h"
#include "molecule.h"

#define MAX_FILE_NAME 1024

enum {
    TARGET_CFOUR,
    TARGET_DIRAC
};

void print_usage(char *binary_name);
void parse_argv(int argc, char **argv, char *file_name, int maxlen, int *target);
void expt_parse(char *file_name, molecule_t *mol, basis_lib_t *basis, ecp_lib_t *ecp_lib);
void gen_pam(FILE *out, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib);
void gen_cfour(FILE *out, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib);


int main(int argc, char **argv)
{
    char inp_file[MAX_FILE_NAME];
    int target = TARGET_DIRAC;
    basis_lib_t *bas;
    ecp_lib_t *ecp;
    molecule_t *mol;

    parse_argv(argc, argv, inp_file, MAX_FILE_NAME, &target);

    bas = basis_lib_new();
    ecp = ecp_lib_new();
    mol = molecule_new();

    expt_parse(inp_file, mol, bas, ecp);

    /*basis_lib_print(bas);
    ecp_lib_print(ecp);
    molecule_print(mol);*/

    if (target == TARGET_DIRAC) {
        gen_pam(stdout, mol, bas, ecp);
    }
    else if (target == TARGET_CFOUR) {
        gen_cfour(stdout, mol, bas, ecp);
    }

    molecule_delete(mol);
    ecp_lib_delete(ecp);
    basis_lib_delete(bas);
}


void parse_argv(int argc, char **argv, char *file_name, int maxlen, int *target)
{
    if (argc == 2) {
        strncpy(file_name, argv[1], maxlen);
        file_name[maxlen - 1] = '\0';
        *target = TARGET_DIRAC;
    }
    else if (argc == 3) {
        strncpy(file_name, argv[2], maxlen);
        file_name[maxlen - 1] = '\0';

        if (strcmp(argv[1], "--pam") == 0 || strcmp(argv[1], "--dirac") == 0) {
            *target = TARGET_DIRAC;
        }
        else if (strcmp(argv[1], "--cfour") == 0) {
            *target = TARGET_CFOUR;
        }
        else {
            print_usage(argv[0]);
            exit(1);
        }
    }
    else {
        print_usage(argv[0]);
        exit(1);
    }
}


void print_usage(char *binary_name)
{
    printf("Usage: %s [target] <input-file>\n", binary_name);
    printf(" target: --dirac or --cfour\n");
}