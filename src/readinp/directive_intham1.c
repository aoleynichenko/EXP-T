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

/*******************************************************************************
 * directive_intham1.c
 *
 * Parser for the 'intham1' block.
 *
 * 2021 Alexander Oleynichenko
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "intham1.h"
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
void parse_sector(char *s, int *h, int *p);

void directive_intham1_sectors(ih1_options_t *ih1_opts);
void directive_intham1_shift_type(ih1_options_t *ih1_opts);
void directive_intham1_npower(ih1_options_t *ih1_opts);
void directive_intham1_main(ih1_options_t *ih1_opts);
void directive_intham1_shift_to(ih1_options_t *ih1_opts);
void directive_intham1_box_line(ih1_options_t *ih1_opts);


/**
 * Syntax:
 * intham1
 *   sectors <list-of-target-sectors>
 *   shift_type [real|realimag|imag|...]
 *   npower <int attenuation_param>
 *   main <int n_main_holes> <int n_main_particles>
 *   shift_to <float upper_for_holes> <float lower_for_particles>
 *   formula [box|line]
 * end
 */
void directive_intham1(cc_options_t *opts)
{
    static char *err_no_sector   = "sector label is expected";
    static char *err_end_of_line = "end of line is expected";
    static char *err_keyword     = "unknown keyword";

    int token_type;
    int main_specified = 0;
    int sectors_specified = 0;
    int shift_to_specified = 0;

    // defaults
    opts->ih1_opts.shift_type = CC_SHIFT_REAL;
    opts->ih1_opts.sum_formula = IH1_SUM_FORMULA_BOX;
    opts->ih1_opts.npower = 3;
    opts->ih1_opts.shift_to_formula = IH1_SHIFT_TO_MAIN_BOUNDS;

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
            directive_intham1_sectors(&opts->ih1_opts);
            sectors_specified = 1;
        }
        else if (token_type == KEYWORD_MAIN) {
            directive_intham1_main(&opts->ih1_opts);
            main_specified = 1;
        }
        else if (token_type == TT_WORD && (strcmp(yytext, "shift_type") == 0)) {
            directive_intham1_shift_type(&opts->ih1_opts);
        }
        else if (token_type == TT_WORD && (strcmp(yytext, "npower") == 0)) {
            directive_intham1_npower(&opts->ih1_opts);
        }
        else if (token_type == TT_WORD && (strcmp(yytext, "shift_to") == 0)) {
            directive_intham1_shift_to(&opts->ih1_opts);
            shift_to_specified = 1;
        }
        else if (token_type == TT_WORD && (strcmp(yytext, "formula") == 0)) {
            directive_intham1_box_line(&opts->ih1_opts);
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
    if (!main_specified) {
        yyerror("main subspace of active spinors must be specified");
    }
    if (!shift_to_specified) {
        yyerror("bounds to shift intermediate spinor energies must be specified");
    }

    opts->do_intham1 = 1;
    opts->shift_type = opts->ih1_opts.shift_type;
}


/**
 * Syntax:
 * sectors <list-of-sector-labels>
 *
 * Example:
 * sectors 0h1p 0h2p
 */
void directive_intham1_sectors(ih1_options_t *ih1_opts)
{
    int token_type;
    int sect_h, sect_p;
    int count = 0;

    while ((token_type = next_token()) == TT_SECTOR) {
        parse_sector(yytext, &sect_h, &sect_p);
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
void directive_intham1_shift_type(ih1_options_t *ih1_opts)
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


/*
 * attenuation parameter
 */
void directive_intham1_npower(ih1_options_t *ih1_opts)
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


void directive_intham1_main(ih1_options_t *ih1_opts)
{
    static char *msg = "wrong specification of the main active space";
    int token_type;

    // number of active "main" hole (occupied) spinors
    token_type = next_token();
    if (token_type == TT_INTEGER) {
        ih1_opts->main_holes = atoi(yytext);
    }
    else {
        yyerror(msg);
    }

    // number of active "main" particle (virtual) spinors
    token_type = next_token();
    if (token_type == TT_INTEGER) {
        ih1_opts->main_particles = atoi(yytext);
    }
    else {
        yyerror(msg);
    }
}


void directive_intham1_shift_to(ih1_options_t *ih1_opts)
{
    static char *err_number  = "floating-point number is expected";
    static char *err_unknown = "unknown option";
    int token_type;

    // "shift_to main"
    // denominators will be shifted to the lowest main hole energy for holes
    // and to the highest main particle energy for particles
    token_type = next_token();
    if (token_type == TT_WORD || token_type == KEYWORD_MAIN) {
        str_tolower(yytext);
        if (strcmp(yytext, "main") == 0) {
            ih1_opts->shift_to_formula = IH1_SHIFT_TO_MAIN_BOUNDS;
            return;
        }
        else if (strcmp(yytext, "barycenter") == 0) {
            ih1_opts->shift_to_formula = IH1_SHIFT_TO_BARYCENTER;
            return;
        }
        else {
            yyerror(err_unknown);
        }
    }
    else {
        put_back(token_type);
    }

    // bound for active holes
    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(err_number);
    }
    ih1_opts->holes_shift_to = atof(yytext);

    // bound for active particles
    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(err_number);
    }
    ih1_opts->shift_to_formula = IH1_SHIFT_TO_GIVEN_BOUNDS;
    ih1_opts->particles_shift_to = atof(yytext);
}


void directive_intham1_box_line(ih1_options_t *ih1_opts)
{
    next_token();

    str_tolower(yytext);

    if (strcmp(yytext, "box") == 0) {
        ih1_opts->sum_formula = IH1_SUM_FORMULA_BOX;
    }
    else if (strcmp(yytext, "line") == 0) {
        ih1_opts->sum_formula = IH1_SUM_FORMULA_LINE;
    }
    else {
        yyerror("unknown formula");
    }
}
