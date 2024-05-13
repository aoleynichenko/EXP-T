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

#include "parse_argv.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Will cause problems on Windows, works for Linux and Mac
#include <getopt.h>

#include "version.h"


/*
 * declarations of functions used in this file
 */
void print_blas_info();

void display_version();

void display_usage();

void display_usage_short();

void display_help();


/**
 * Parses command-line arguments:
 * -s, --scratch=PATH         Path to scratch directory (default: ./scratch)
 * -n, --no-clean             Do not clean scratch directory on exit
 * -?, --help                 Print this help list and exit
 *     --usage                Print a short usage message and exit
 * -V, --version              Print program version and exit
 * Usage: expt.x [-n?V] [-s PATH] [--no-clean-scratch] [--scratch-dir=PATH]
 *               [--help] [--usage] [--version] <input-file>
 */
void parse_argv(int argc, char **argv, char *input_name, char **scratch_dir_path, int *do_clean_scratch)
{
    int args_clean;              /* =0 if do not clean scratch directory on exit */
    char *args_scratch_dir;      /* Path to scratch directory (default: ./scratch) */
    char **args_inputFiles;      /* input files */
    int args_numInputFiles;      /* number of input files */
    int args_print_usage;        /* --usage option */

    int opt = 0;
    int longIndex = 0;

    // tells getopt_long() which options we support, and which options have arguments
    static const char *optString = "nIl:o:s:vh?";

    // tells getopt_long() which long options we support, and which long options have arguments
    static const struct option longOpts[] = {
            {"no-clean", no_argument,       NULL, 'n'},
            {"version",  no_argument,       NULL, 'V'},
            {"scratch",  required_argument, NULL, 's'},
            {"usage",    no_argument,       NULL, 0},
            {"help",     no_argument,       NULL, 'h'},
            {NULL,       no_argument,       NULL, 0}
    };

    /* Initialize args before we get to work. */
    args_inputFiles = NULL;
    args_scratch_dir = NULL;
    args_numInputFiles = 0;
    args_print_usage = 0;
    args_clean = 1;

    /* Process the arguments with getopt_long(), then
     * populate args_*.
     */
    opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    while (opt != -1) {
        switch (opt) {
            case 'n':
                args_clean = 0;    /* true */
                break;
            case 'V':
                display_version();
                exit(0);
                break;
            case 's':
                args_scratch_dir = optarg;
                break;
            case 'h':    /* fall-through is intentional */
            case '?':
                display_help();
                exit(0);
                break;
            case 0:        /* long option without a short arg */
                if (strcmp("usage", longOpts[longIndex].name) == 0) {
                    args_print_usage = 1;
                    display_usage();
                    exit(0);
                }
                break;
            default:
                /* You won't actually get here. */
                break;
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }

    args_inputFiles = argv + optind;
    args_numInputFiles = argc - optind;

    if (args_numInputFiles != 1) {
        display_usage_short();
        exit(0);
    }

    if (args_print_usage) {
        display_usage();
        exit(0);
    }

    strcpy(input_name, args_inputFiles[0]);
    *scratch_dir_path = args_scratch_dir;
    *do_clean_scratch = args_clean;
}


void display_version()
{
    printf("expt v. %d.%d.%d (%d %s %d) public\n",
           CC_VERSION_MAJOR, CC_VERSION_MINOR, CC_VERSION_REVISION,
           CC_VERSION_DAY, CC_VERSION_MONTH, CC_VERSION_YEAR);

#if defined __ICC
    printf("Intel C Compiler %d (%s)\n", __ICC, __VERSION__);
#elif defined __GNUC__
    printf("GNU C Compiler %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#else
    printf("undetected compiler\n");
#endif
    print_blas_info();
    printf("Build date: %s %s\n", __DATE__, __TIME__);
    exit(0);
}


void display_usage()
{
    printf("Usage: expt.x [-n?V] [-s PATH] [--no-clean] [--scratch=PATH]\n");
    printf("       [--help] [--usage] [--version] <input-file>\n");
    exit(0);
}


void display_usage_short()
{
    printf("Usage: expt.x [OPTION...] <input-file>\n");
    printf("Try `expt.x --help' or `expt.x --usage' for more information.\n");
    exit(0);
}


void display_help()
{
    printf("Usage: expt.x [OPTION...] <input-file>\n");
    printf("expt -- A Fock-Space Multireference Relativistic Coupled-Cluster Program\n");
    printf("\n");
    printf("  -n, --no-clean             Do not clean scratch directory on exit\n");
    printf("                             (use this option to preserve cluster amplitudes etc)\n");
    printf("  -s, --scratch=PATH         Path to scratch directory (default: ./scratch)\n");
    printf("  -?, --help                 Print this help list and exit\n");
    printf("      --usage                Print a short usage message and exit\n");
    printf("  -V, --version              Print program version and exit\n");
    printf("\n");
    printf("Mandatory or optional arguments to long options are also mandatory or optional\n");
    printf("for any corresponding short options.\n");
    printf("\n");
    printf("Please report bugs to <exp-t-program@googlegroups.com>.\n");
    exit(0);
}




