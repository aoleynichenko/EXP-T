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

/**
 * Parser for the 'ih_imms' block.
 *
 * The 'ih_imms' block specifies parameters for the Intermediate Hamiltonian
 * Fock space coupled cluster method for incomplete main model spaces.
 * For details, see:
 * A. V. Zaitsevskii, N. S. Mosyagin, A. V. Oleynichenko, E. Eliav,
 * Int. J. Quantum Chem, submitted (2022)
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "intham_imms.h"
#include "lexer.h"
#include "options.h"

void yyerror(char *s);

int next_token();

void put_back(int token_type);

int match(int required_type);

void str_tolower(char *s);

int parse_state_spec(char *buf, char *rep_name, int *state_no);

void read_space_specification(cc_space_t *space);

void set_input_file_name(char *file_name);

void parse_sector_label(char *s, int *h, int *p);

void directive_ih_imms_target_sectors(ih_imms_options_t *ih1_opts);

void directive_ih_imms_shift_type(ih_imms_options_t *ih1_opts);

void directive_ih_imms_npower(ih_imms_options_t *ih1_opts);

void directive_ih_imms_det_energy(ih_imms_options_t *ih1_opts);

void directive_ih_imms_subspace(ih_imms_options_t *ih1_opts);

void directive_ih_imms_main_occ(ih_imms_options_t *ih1_opts);

void directive_ih_imms_scale_shift(ih_imms_options_t *ih1_opts);


/**
 * Syntax:
 * ih_imms
 *   sectors <list-of-target-sectors>
 *   shift_type [real|realimag|imag|...]
 *   npower <int attenuation_param>
 *   subspace energy <float lower_bound> <float upper_bound>
 *   main_occ <list-of-integers>
 *   frontier_energy <real energy>
 * end
 */
void directive_ih_incomplete_main_model_spaces(cc_options_t *opts)
{
    static char *err_end_of_line = "end of line is expected";
    static char *err_keyword = "unknown keyword";

    int token_type;
    int main_specified = 0;
    int sectors_specified = 0;
    int shift_to_specified = 0;
    int subspaces_specified = 0;

    // defaults
    opts->intham_imms_opts.shift_type = CC_SHIFT_REALIMAG;
    opts->intham_imms_opts.npower = 5;
    opts->intham_imms_opts.subspaces_definition = IH_IMMS_SUBSPACES_DEF_ENERGY;

    // end of line after the opening keyword is required
    token_type = next_token();
    if (token_type != END_OF_LINE) {
        yyerror(err_end_of_line);
    }

    // loop over sub-directives
    token_type = next_token();
    while (token_type != END_OF_FILE) {

        str_tolower(yytext);

        if (token_type == TT_WORD && (strcmp(yytext, "sectors") == 0)) {
            directive_ih_imms_target_sectors(&opts->intham_imms_opts);
            sectors_specified = 1;
        }
        else if (token_type == TT_WORD && (strcmp(yytext, "main_occ") == 0)) {
            directive_ih_imms_main_occ(&opts->intham_imms_opts);
            main_specified = 1;
        }
        else if (token_type == KEYWORD_SHIFT_TYPE) {
            directive_ih_imms_shift_type(&opts->intham_imms_opts);
        }
        else if (token_type == TT_WORD && (strcmp(yytext, "npower") == 0)) {
            directive_ih_imms_npower(&opts->intham_imms_opts);
        }
        else if (token_type == TT_WORD && (strcmp(yytext, "scale") == 0)) {
            directive_ih_imms_scale_shift(&opts->intham_imms_opts);
        }
        else if (token_type == TT_WORD && (strcmp(yytext, "subspace") == 0)) {
            directive_ih_imms_subspace(&opts->intham_imms_opts);
            subspaces_specified = 1;
        }
        else if (token_type == TT_WORD &&
                 (strcmp(yytext, "det_energy") == 0 ||
                  strcmp(yytext, "frontier_energy") == 0)) {
            directive_ih_imms_det_energy(&opts->intham_imms_opts);
            shift_to_specified = 1;
        }
        else if (token_type == KEYWORD_END) {
            break;
        }
        else if (token_type == END_OF_LINE) {
            // nothing to do
        }
        else {
            yyerror(err_keyword);
        }

        /* each "non-newline" directive must end with End-Of-Line */
        if (token_type != END_OF_LINE) {
            token_type = next_token();
            if (token_type != END_OF_LINE && token_type != END_OF_FILE) {
                yyerror(err_end_of_line);
            }
        }

        /* go to the next directive */
        token_type = next_token();
    }

    if (!sectors_specified) {
        yyerror("at least one target sector for IH must be specified");
    }
    if (!main_specified &&
        !shift_to_specified/* && opts->intham_imms_opts.sum_formula != IH1_SHIFT_FORMULA_DET_ENERGY*/) {
        yyerror("main subspace of active spinors must be specified");
    }
    if (!subspaces_specified &&
        !shift_to_specified/* && opts->intham_imms_opts.sum_formula != IH1_SHIFT_FORMULA_DET_ENERGY*/) {
        yyerror("at least one spinor subspace must be specified");
    }

    /*
     * post-processing of IH1 options and re-check
     */
    if (!(opts->intham_imms_opts.shift_type == CC_SHIFT_REALIMAG ||
          opts->intham_imms_opts.shift_type == CC_SHIFT_IMAG) &&
        fabs(opts->intham_imms_opts.scale_shift - 1.0) > 1e-10) {
        yyerror("scaling factors can be used with imaginary shifts only");
    }

    opts->do_intham_imms = 1;
    opts->shift_type = opts->intham_imms_opts.shift_type;
}


/**
 * Syntax:
 * sectors <list-of-sector-labels>
 *
 * Example:
 * sectors 0h1p 0h2p
 */
void directive_ih_imms_target_sectors(ih_imms_options_t *ih1_opts)
{
    int token_type;
    int sect_h, sect_p;
    int count = 0;

    while ((token_type = next_token()) == TT_SECTOR) {
        parse_sector_label(yytext, &sect_h, &sect_p);
        if (sect_h >= MAX_SECTOR_RANK || sect_p >= MAX_SECTOR_RANK) {
            yyerror("sector is too high");
        }
        ih1_opts->sectors[sect_h][sect_p] = 1;
        count++;
    }

    if (count == 0) {
        yyerror("at least one sector label is required");
    }

    put_back(token_type);
}


/**
 * Syntax:
 * shift_type ( none || real || realimag || imag || taylor )
 */
void directive_ih_imms_shift_type(ih_imms_options_t *ih1_opts)
{
    static char *msg = "wrong specification of the denominators shift formula!\n"
                       "Possible values: none, taylor, real, realimag, imag";
    cc_shifttype_t type;

    if (!match(TT_WORD)) {
        yyerror(msg);
    }
    str_tolower(yytext);
    if (strcmp(yytext, "none") == 0) {
        type = CC_SHIFT_NONE;
    }
    else if (strcmp(yytext, "taylor") == 0) {
        type = CC_SHIFT_TAYLOR;
    }
    else if (strcmp(yytext, "real") == 0) {
        type = CC_SHIFT_REAL;
    }
    else if (strcmp(yytext, "realimag") == 0) {
        type = CC_SHIFT_REALIMAG;
    }
    else if (strcmp(yytext, "imag") == 0) {
        type = CC_SHIFT_IMAG;
    }
    else {
        yyerror(msg);
    }

    ih1_opts->shift_type = type;
}


/**
 * attenuation parameter
 */
void directive_ih_imms_npower(ih_imms_options_t *ih1_opts)
{
    static char *msg = "wrong specification of the attenuation parameter";
    int token_type;

    token_type = next_token();
    if (token_type == TT_INTEGER) {
        ih1_opts->npower = atoi(yytext);
    }
    else {
        yyerror(msg);
    }
}


/**
 * specification of frontier energy in a given FS sector
 */
void directive_ih_imms_det_energy(ih_imms_options_t *ih1_opts)
{
    static char *msg1 = "wrong specification of the shift parameter\n"
                        "allowed values: lower_bound, upper_bound, real number";
    int sect_h = 0, sect_p = 0;

    // first: sector label
    int token_type = next_token();
    if (token_type == TT_SECTOR) {
        parse_sector_label(yytext, &sect_h, &sect_p);
        if (sect_h >= MAX_SECTOR_RANK || sect_p >= MAX_SECTOR_RANK) {
            yyerror("sector is too high");
        }
    }
    else {
        yyerror("sector label is required");
    }

    // second: frontier energy for intermediate space determinants
    token_type = next_token();
    if (token_type == TT_WORD) {
        str_tolower(yytext);
        if (strcmp(yytext, "lower_bound") == 0) {
            ih1_opts->det_shift_auto[sect_h][sect_p] = IH_IMMS_FRONTIER_ENERGY_LOWER_BOUND;
        }
        else if (strcmp(yytext, "upper_bound") == 0) {
            ih1_opts->det_shift_auto[sect_h][sect_p] = IH_IMMS_FRONTIER_ENERGY_UPPER_BOUND;
        }
        else {
            yyerror(msg1);
        }
    }
    else if (token_type == TT_INTEGER || token_type == TT_FLOAT) {
        ih1_opts->det_shift_to[sect_h][sect_p] = atof(yytext);
        ih1_opts->det_shift_auto[sect_h][sect_p] = IH_IMMS_FRONTIER_ENERGY_MANUAL;
    }
    else {
        yyerror(msg1);
    }
}


/**
 * defines subspaces of one-electron functions (spinors)
 */
void directive_ih_imms_subspace(ih_imms_options_t *ih1_opts)
{
    static char *msg1 = "wrong specification of active space!\n"
                        "Two real numbers are expected";
    static char *msg2 = "wrong specification of active space!\n"
                        "expected the 'energy' keyword";

    int token_type = next_token();

    if (token_type == TT_WORD) {
        str_tolower(yytext);
        if (strcmp(yytext, "energy") != 0) {
            yyerror(msg2);
        }

        // lower limit
        token_type = next_token();
        if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
            yyerror(msg1);
        }
        double emin = atof(yytext);

        // upper limit
        token_type = next_token();
        if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
            yyerror(msg1);
        }
        double emax = atof(yytext);

        ih1_opts->subspaces_definition = IH_IMMS_SUBSPACES_DEF_ENERGY;
        int idx = ih1_opts->n_spinor_subspaces;
        ih1_opts->subspace_energy_ranges[idx][0] = emin;
        ih1_opts->subspace_energy_ranges[idx][1] = emax;
        ih1_opts->n_spinor_subspaces++;
    }
    else if (token_type == TT_INTEGER) {
        printf("number of spinors in the subspace = %s\n", yytext);

        ih1_opts->subspaces_definition = IH_IMMS_SUBSPACES_DEF_NTOTAL;
        int idx = ih1_opts->n_spinor_subspaces;
        ih1_opts->subspace_total_nspinors[idx] = atoi(yytext);
        ih1_opts->n_spinor_subspaces++;
    }
    else {
        yyerror(msg2);
    }
}


/**
 * defines main model space determinants
 */
void directive_ih_imms_main_occ(ih_imms_options_t *ih1_opts)
{
    static char *msg1 = "wrong specification of active space!\n"
                        "Two real numbers are expected";
    static char *msg2 = "wrong specification of active space!\n"
                        "expected the 'energy' keyword";

    int n_ints_read = 0;
    int buf[IH_IMMS_MAX_SPINOR_SUBSPACES];

    // read list of integers
    int token_type = next_token();
    while (token_type == TT_INTEGER) {
        buf[n_ints_read++] = atoi(yytext);
        token_type = next_token();
    }
    put_back(token_type);

    // save the acquired list
    int imain = ih1_opts->n_main_subspaces;
    for (int i = 0; i < n_ints_read; i++) {
        ih1_opts->main_occ[imain][i] = buf[i];
    }
    ih1_opts->n_main_subspaces++;
}


/**
 * defines scaling factor for the determinantal shift
 */
void directive_ih_imms_scale_shift(ih_imms_options_t *ih1_opts)
{
    static char *msg1 = "wrong specification of the scaling parameter\n"
                        "Real number is expected";

    int token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg1);
    }

    ih1_opts->scale_shift = atof(yytext);
}
