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

/*
 * Fock-space multireference coupled cluster program.
 *
 * Authors:
 *   Alexander Oleynichenko (NRC "Kurchatov Institute" - PNPI, Gatchina, Russia)
 * Mailto:
 *   alexvoleynichenko@gmail.com
 *   oleynichenko_av@pnpi.nrcki.ru
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "parse_argv.h"
#include "platform.h"
#include "engine.h"
#include "io.h"
#include "interfaces.h"
#include "methods.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "symmetry.h"
#include "timer.h"
#include "utils.h"
#include "version.h"


/*
 * declarations of functions used in this file
 */
int readinp(char *file_name, cc_options_t *opts);

void setup_scratch();

int cuda_device_query();

void print_header_banner();

void print_footer_banner(int exit_code);


/**
 * entry point
 */
int main(int argc, char **argv)
{
    cc_options_t *opts;
    char input_name[256];
    int exit_code = EXIT_SUCCESS;
    int do_clean_scratch = 0;
    char *scratch_dir_path = NULL;

    // turn off buffering for stdout (for seeing output immediately)
    setvbuf(stdout, NULL, _IONBF, 0);

    // parse command line arguments
    parse_argv(argc, argv, input_name, &scratch_dir_path, &do_clean_scratch);

    print_header_banner();

    // ???
    execute_external_program("rm", "-rf HEFF", NULL);

    // setup timers
    timer_new_entry("tot", "Total time of execution");
    timer_start("tot");
    timer_new_entry("mult_pppp", "Diagram contraction (mult) for <PP||PP>");

    // read input file
    opts = new_options();
    cc_opts = opts;
    readinp(input_name, opts);

    // apply options from command line
    if (scratch_dir_path != NULL) {
        strcpy(opts->scratch_dir, scratch_dir_path);
    }
    opts->clean_scratch = do_clean_scratch;

    print_options(opts);

    cc_init_allocator(cc_opts->max_memory_size);

    // setup scratch directory: create if needed and cd to it
    setup_scratch();

    // load transformed MO integrals
    // prepare raw integrals for sorting
    if (opts->int_source == CC_INTEGRALS_DIRAC) {
        dirac_interface(opts);
    }
    else {
        printf("unknown integral interface\n");
        return 1;
    }

    // solve FS-CC equations in the required sector
    if (opts->sector_h == 0 && opts->sector_p == 0) {
        exit_code = sector_0h0p(opts);
    }
    else if (opts->sector_h == 0 && opts->sector_p == 1) {
        exit_code = sector_0h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_0h1p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
    }
    else if (opts->sector_h == 1 && opts->sector_p == 0) {
        exit_code = sector_0h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_1h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
    }
    else if (opts->sector_h == 1 && opts->sector_p == 1) {
        exit_code = sector_0h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_1h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_0h1p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_1h1p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
    }
    else if (opts->sector_h == 1 && opts->sector_p == 2) {
        exit_code = sector_0h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_1h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_0h1p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_1h1p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_0h2p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_1h2p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
    }
    else if (opts->sector_h == 0 && opts->sector_p == 2) {
        exit_code = sector_0h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_0h1p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_0h2p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
    }
    else if (opts->sector_h == 2 && opts->sector_p == 0) {
        exit_code = sector_0h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_1h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_2h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
    }
    else if (opts->sector_h == 0 && opts->sector_p == 3) {
        exit_code = sector_0h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_0h1p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_0h2p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_0h3p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
    }
    else if (opts->sector_h == 3 && opts->sector_p == 0) {
        exit_code = sector_0h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_1h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_2h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
        exit_code = sector_3h0p(opts);
        if (exit_code == EXIT_FAILURE) {
            goto finalize;
        }
    }
    else {
        errquit("Fock space sector %dh%dp is not implemented yet", opts->sector_h, opts->sector_p);
    }

    finalize:

    timer_stop("tot");
    if (opts->print_level >= CC_PRINT_MEDIUM) {
        timer_stats();
    }

    if (opts->compress != CC_COMPRESS_NONE && opts->print_level >= CC_PRINT_HIGH) {
        print_compression_stats();
    }

    // final clean-up and exit
    delete_options(opts);

    symmetry_cleanup();
    spinors_cleanup();
    cc_finalize_allocator();

    // clean working directory
    if (do_clean_scratch) {
        execute_external_program("rm", "-rf VINT* *sb *dg HINT", NULL);
    }

    // print total time of execution: days, hours, min, sec, millisec
    print_footer_banner(exit_code);
}


/**
 * Creates scratch directory for temporary files (if needed)
 */
void setup_scratch()
{
    int status;
    char *scratch_dir = cc_opts->scratch_dir;
    char cwd_path[CC_MAX_PATH_LENGTH];
    char tmp[CC_MAX_PATH_LENGTH];

    if (!io_directory_exists(scratch_dir)) {
        printf(" Creating scratch directory '%s' ...\n", scratch_dir);
        status = io_mkdir(scratch_dir);
        if (status != 0) {
            errquit("Unable to create scratch directory "
                    "(errno = %d, strerror = %s)", errno, strerror(errno));
        }
    }
    else {
        printf(" Scratch directory '%s' exists\n", scratch_dir);
    }

    io_getcwd(cwd_path, CC_MAX_PATH_LENGTH);
    printf(" Current working directory = %s\n", cwd_path);

    /*
     * determine absolute path to files with integrals:
     */

    // "MRCONEE" - 1el integrals
    strcpy(tmp, cwd_path);
    strcat(tmp, "/");
    strcat(tmp, cc_opts->integral_file_1);
    strcpy(cc_opts->integral_file_1, tmp);
    printf(" Full path to one-electron integrals = %s\n", cc_opts->integral_file_1);

    // "MDCINT" - 2el (Coulomb) integrals
    strcpy(tmp, cwd_path);
    strcat(tmp, "/");
    strcat(tmp, cc_opts->integral_file_2);
    strcpy(cc_opts->integral_file_2, tmp);
    printf(" Full path to Coulomb integrals = %s\n", cc_opts->integral_file_2);

    // "MDPROP" - properties integrals
    strcpy(tmp, cwd_path);
    strcat(tmp, "/");
    strcat(tmp, cc_opts->integral_file_prop);
    strcpy(cc_opts->integral_file_prop, tmp);
    printf(" Full path to properties integrals = %s\n", cc_opts->integral_file_prop);

    // "MGINT" - Gaunt integrals
    if (cc_opts->gaunt_defined) {
        strcpy(tmp, cwd_path);
        strcat(tmp, "/");
        strcat(tmp, cc_opts->integral_file_gaunt);
        strcpy(cc_opts->integral_file_gaunt, tmp);
        printf(" Full path to Gaunt integrals = %s\n", cc_opts->integral_file_gaunt);
    }

    // OneProp files (1el integrals, real and imaginary parts)
    if (cc_opts->oneprop_on) {
        for (int i = 0; i < cc_opts->n_oneprop; i++) {
            strcpy(tmp, cwd_path);
            strcat(tmp, "/");
            strcat(tmp, cc_opts->oneprop_file_re[i]);
            strcpy(cc_opts->oneprop_file_re[i], tmp);
            printf(" Full path to OneProp 1-el integrals (real part) = %s\n", cc_opts->oneprop_file_re[i]);
            strcpy(tmp, cwd_path);
            strcat(tmp, "/");
            strcat(tmp, cc_opts->oneprop_file_im[i]);
            strcpy(cc_opts->oneprop_file_im[i], tmp);
            printf(" Full path to OneProp 1-el integrals (imag part) = %s\n", cc_opts->oneprop_file_im[i]);
        }
    }

    if (cc_opts->n_model_space_props > 0) {
        for (int i = 0; i < cc_opts->n_model_space_props; i++) {
            if (cc_opts->prop_queries[i].source == CC_PROP_FROM_TXTPROP) {
                strcpy(tmp, cwd_path);
                strcat(tmp, "/");
                strcat(tmp, cc_opts->prop_queries[i].file_real);
                strcpy(cc_opts->prop_queries[i].file_real, tmp);
                printf(" Full path to OneProp 1-el integrals (real part) = %s\n", cc_opts->prop_queries[i].file_real);
                strcpy(tmp, cwd_path);
                strcat(tmp, "/");
                strcat(tmp, cc_opts->prop_queries[i].file_imag);
                strcpy(cc_opts->prop_queries[i].file_imag, tmp);
                printf(" Full path to OneProp 1-el integrals (imag part) = %s\n", cc_opts->prop_queries[i].file_imag);
            }
        }
    }

    // TwoProp files (2el property integrals)
    if (cc_opts->twoprop_on) {
        for (int i = 0; i < cc_opts->n_twoprop; i++) {
            strcpy(tmp, cwd_path);
            strcat(tmp, "/");
            strcat(tmp, cc_opts->twoprop_file[i]);
            strcpy(cc_opts->twoprop_file[i], tmp);
            printf(" Full path to 2-el property integrals = %s\n", cc_opts->twoprop_file[i]);
        }
    }

    // change working directory
    printf(" Changing working directory to %s ...\n", scratch_dir);
    status = io_chdir(scratch_dir);
    if (status != 0) {
        errquit("Unable to use scratch directory '%s' as working directory "
                "(errno = %d, strerror = %s)", scratch_dir, errno, strerror(errno));
    }
    io_getcwd(cwd_path, CC_MAX_PATH_LENGTH);
    printf(" Current working directory = %s\n", cwd_path);
}


/**
 * Prints header banner containing title, authors and basic compilation info.
 */
void print_header_banner()
{
    time_t t = time(0);
    const size_t MAX_HOST_NAME_LEN = 1024;
    char host_name[MAX_HOST_NAME_LEN];

    printf("\n");
    printf("\t\t**********************************************************************************\n");
    printf("\t\t**                                                                              **\n");
    printf("\t\t**                                   E X P - T                                  **\n");
    printf("\t\t**        Relativistic Fock-Space Multireference Coupled Cluster Program        **\n");
    printf("\t\t**                                                                              **\n");
    printf("\t\t**                         version %d.%d.%d (%d %s %d)                         **\n",
           CC_VERSION_MAJOR, CC_VERSION_MINOR, CC_VERSION_REVISION,
           CC_VERSION_DAY, CC_VERSION_MONTH, CC_VERSION_YEAR);
    printf("\t\t**                                                                              **\n");
    printf("\t\t**********************************************************************************\n");
    printf("\t\t**                                                                              **\n");
    printf("\t\t** EXP-T is free software: you can redistribute it and/or modify                **\n");
    printf("\t\t** it under the terms of the GNU Lesser General Public License as published by  **\n");
    printf("\t\t** the Free Software Foundation, either version 3 of the License, or            **\n");
    printf("\t\t** (at your option) any later version.                                          **\n");
    printf("\t\t**                                                                              **\n");
    printf("\t\t** EXP-T is distributed in the hope that it will be useful,                     **\n");
    printf("\t\t** but WITHOUT ANY WARRANTY; without even the implied warranty of               **\n");
    printf("\t\t** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                **\n");
    printf("\t\t** GNU Lesser General Public License for more details.                          **\n");
    printf("\t\t**                                                                              **\n");
    printf("\t\t** You should have received a copy of the GNU Lesser General Public License     **\n");
    printf("\t\t** along with EXP-T.  If not, see <http://www.gnu.org/licenses/>.               **\n");
    printf("\t\t**                                                                              **\n");
    printf("\t\t**********************************************************************************\n");
    printf("\t\t**                                                                              **\n");
    printf("\t\t** For more information about the EXP-T program system see                      **\n");
    printf("\t\t** http://www.qchem.pnpi.spb.ru/expt                                            **\n");
    printf("\t\t**                                                                              **\n");
    printf("\t\t** This is an experimental code. The authors accept no responsibility           **\n");
    printf("\t\t** for the performance of the code or for the correctness of the results.       **\n");
    printf("\t\t**                                                                              **\n");
    printf("\t\t** If results obtained with this code are published, an                         **\n");
    printf("\t\t** appropriate citation would be:                                               **\n");
    printf("\t\t**                                                                              **\n");
    printf("\t\t** [1] A. Oleynichenko, A. Zaitsevskii, E. Eliav.                               **\n");
    printf("\t\t**     EXP-T, An Extensible Code for Fock Space Relativistic Coupled Cluster    **\n");
    printf("\t\t**     Calculations. See http://www.qchem.pnpi.spb.ru/expt.                     **\n");
    printf("\t\t** [2] A. V. Oleynichenko, A. Zaitsevskii, E. Eliav.                            **\n");
    printf("\t\t**     Towards high performance relativistic electronic structure modelling:    **\n");
    printf("\t\t**     the EXP-T program package. Commun. Comput. Inf. Sci. 1331, 375 (2020).   **\n");
    printf("\t\t**     doi: 10.1007/978-3-030-64616-5_33                                        **\n");
    printf("\t\t**                                                                              **\n");
    printf("\t\t**********************************************************************************\n");
    printf("\n");

    printf("Authors:\n");
    printf("  Alexander Oleynichenko [alexvoleynichenko@gmail.com]\n");
    printf("  Andrei Zaitsevskii\n");
    printf("  Ephraim Eliav\n");
#if defined __ICC
    printf("Compiler: Intel C Compiler %d (%s)\n", __ICC, __VERSION__);
#elif defined __GNUC__
    printf("Compiler: GNU C Compiler %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#else
    printf("Compiler: undetected\n");
#endif
    print_blas_info();
    printf("Build date: %s %s\n", __DATE__, __TIME__);
    printf("Run date: %s", asctime(localtime(&t)));
    get_host_name(host_name, MAX_HOST_NAME_LEN);
    printf("Hostname: %s", host_name);
    printf("\n");

    // CUDA settings
#ifdef CUDA_FOUND
    cuda_device_query();
#else
    printf("CUDA disabled\n");
#endif
}


/**
 * prints final information on exit:
 * termination code, wall time, amount of memory used
 */
void print_footer_banner(int exit_code)
{
    int d, h, m, s, ms;
    double total_time;
    io_stat_t io_st;
    time_t t = time(0);

    // disk usage report
    io_statistics(&io_st);
    printf(" Disk I/O:\n");
    printf("   files created: %ld   files removed: %ld\n", io_st.n_created, io_st.n_removed);
    printf("   read   %15ld bytes = %.3f Gb\n", io_st.n_read, io_st.n_read / (1024.0 * 1024.0 * 1024.0));
    printf("   write  %15ld bytes = %.3f Gb\n", io_st.n_written, io_st.n_written / (1024.0 * 1024.0 * 1024.0));
    printf("\n");

    // calculate number of days, hours, minutes, seconds, milliseconds
    total_time = timer_get("tot");
    decompose_time(total_time, &d, &h, &m, &s, &ms);

    if (exit_code == EXIT_SUCCESS) {
        printf(" EXP-T terminated normally at %s", asctime(localtime(&t)));
    }
    else {
        printf(" EXP-T terminated ABNORMALLY at %s", asctime(localtime(&t)));
    }
    printf(" Total run time: %d days %d hours %d minutes %d seconds %d milliseconds\n", d, h, m, s, ms);
}
