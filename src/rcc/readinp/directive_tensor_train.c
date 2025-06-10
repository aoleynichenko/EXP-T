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
 * Parser for the 'tensor_trains' block.
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

double match_float_number();

double match_positive_float_number();

int match_positive_integer();

int match_non_negative_integer();

void match_end_of_line();


/**
 * Syntax:
 *
 * tensor_trains
 *   accuracy ( low | medium | high )
 *   svd_thresh <real thresh>
 *   cholesky_thresh <real thresh>
 *   no_tt
 * end
 */
void directive_tensor_train(cc_options_t *opts)
{
    static char *err_end_of_line = "end of line is expected";
    static char *err_keyword = "unknown keyword";

    // defaults
    opts->tt_options.tt_module_enabled = 1;
    opts->tt_options.use_tensor_trains = 1;
    opts->tt_options.svd_thresh = 1e-9;
    opts->tt_options.cholesky_thresh = 1e-6;
    opts->tt_options.use_pyscf_integrals = 0;
    opts->tt_options.pyscf_integrals_path[0] = '\0';

    // end of line after the opening keyword is required
    match_end_of_line();

    // loop over sub-directives
    int token_type = next_token();
    while (token_type != END_OF_FILE) {

        str_tolower(yytext);

        if (token_type == TT_WORD && (strcmp(yytext, "no_tt") == 0)) {
            opts->tt_options.use_tensor_trains = 0;
        }
        else if (token_type == TT_WORD && (strcmp(yytext, "pyscf") == 0)) {
            opts->tt_options.use_pyscf_integrals = 1;
            token_type = next_token();
            if (token_type == END_OF_LINE || token_type == END_OF_FILE) {
                yyerror("path to the file with pyscf integrals is expected");
            }
            strcpy(opts->tt_options.pyscf_integrals_path, yytext);
        }
        else if (token_type == TT_WORD && (strcmp(yytext, "svd_thresh") == 0)) {
            opts->tt_options.svd_thresh = match_positive_float_number();
        }
        else if (token_type == TT_WORD && (strcmp(yytext, "cholesky_thresh") == 0)) {
            opts->tt_options.cholesky_thresh = match_positive_float_number();
        }
        else if (token_type == TT_WORD && (strcmp(yytext, "accuracy") == 0)) {
            token_type = next_token();

            if (token_type == TT_WORD) {
                str_tolower(yytext);
                if (strcmp(yytext, "low") == 0) {
                    opts->tt_options.svd_thresh = 1e-6;
                    opts->tt_options.cholesky_thresh = 1e-3;
                }
                else if (strcmp(yytext, "medium") == 0) {
                    opts->tt_options.svd_thresh = 1e-8;
                    opts->tt_options.cholesky_thresh = 1e-4;
                }
                else if (strcmp(yytext, "high") == 0) {
                    opts->tt_options.svd_thresh = 1e-9;
                    opts->tt_options.cholesky_thresh = 1e-6;
                }
                else {
                    yyerror("wrong accuracy! Possible values are: low, medium, high");
                }
            }
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
}
