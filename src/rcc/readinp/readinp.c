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
 * Reads EXP-T input file.
 * Lexical analysis is performed with the automatically generated code
 * (by lex, see expt.l)
 */

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "codata.h"
#include "lexer.h"
#include "platform.h"
#include "error.h"
#include "options.h"
#include "spinors.h"
#include "memory.h"
#include "cc_properies.h"

#define MAX_LINE_LEN 1024
#define MAXIMUM2(x, y) (((x) > (y)) ? (x) : (y))

/*
 * functions defined in this file
 */
int echo(char *file_name, FILE *dst);

int readinp(char *file_name, cc_options_t *opts);

void directive_title(cc_options_t *opts);

void directive_print(cc_options_t *opts);

void directive_scratch_dir(cc_options_t *opts);

void directive_sector(cc_options_t *opts);

void directive_model(cc_options_t *opts);

void directive_occ(cc_options_t *opts);

void directive_occ_irreps(cc_options_t *opts);

void directive_active(cc_options_t *opts);

void directive_nactp(cc_options_t *opts);

void directive_nacth(cc_options_t *opts);

void directive_degen_thresh(cc_options_t *opts);

void directive_natorb(cc_options_t *opts);

void directive_nroots(cc_options_t *opts);

void directive_roots_cutoff(cc_options_t *opts);

void directive_maxiter(cc_options_t *opts);

void directive_conv(cc_options_t *opts);

void directive_div_thresh(cc_options_t *opts);

void directive_damping(cc_options_t *opts);

void directive_diis(cc_options_t *opts);

void directive_crop(cc_options_t *opts);

void directive_shift_type(cc_options_t *opts);

void directive_shift(cc_options_t *opts);

void directive_reuse(cc_options_t *opts);

void directive_skip(cc_options_t *opts);

void directive_flush(cc_options_t *opts);

void directive_integrals(cc_options_t *opts);

void directive_gaunt(cc_options_t *opts);

void directive_oneprop(cc_options_t *opts);

void directive_twoprop(cc_options_t *opts);

void directive_memory(cc_options_t *opts);

void directive_tilesize(cc_options_t *opts);

void directive_disk_usage(cc_options_t *opts);

void directive_nthreads(cc_options_t *opts);

void directive_openmp_algorithm(cc_options_t *opts);

void directive_arith(cc_options_t *opts);

void directive_mdprop(cc_options_t *opts);

void directive_txtprop(cc_options_t *opts);

void directive_analyt_prop(cc_options_t *opts);

void directive_density(cc_options_t *opts);

void directive_lambda(cc_options_t *opts);

void directive_overlap(cc_options_t *opts);

void directive_interface(cc_options_t *opts);

void directive_ih_incomplete_main_model_spaces(cc_options_t *opts);

void directive_restrict_t3(cc_options_t *opts);

void directive_compress_triples(cc_options_t *opts);

void directive_spinor_labels(cc_options_t *opts);

void directive_tt_svd_tol(cc_options_t *opts);

void directive_tt_cholesky_tol(cc_options_t *opts);

void directive_tensor_train(cc_options_t *opts);

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
 * Read options from input file, store them in object 'opts'
 */
int readinp(char *file_name, cc_options_t *opts)
{
    int token_type;

    yyin = fopen(file_name, "r");
    if (yyin == NULL) {
        errquit("Input file \"%s\" not found", file_name);
    }

    printf("Reading input file: %s\n\n", file_name);
    printf("\n");
    printf("\t\t\t\t****************\n");
    printf("\t\t\t\t** Input file **\n");
    printf("\t\t\t\t****************\n\n");
    printf(" -----\n");
    echo(file_name, stdout);
    printf(" -----\n");
    printf("\n\n");

    set_input_file_name(file_name);

    token_type = next_token();
    while (token_type != END_OF_FILE) {
        switch (token_type) {
            case KEYWORD_TITLE:
                directive_title(opts);
                break;
            case KEYWORD_PRINT:
                directive_print(opts);
                break;
            case KEYWORD_SCRATCH_DIR:
                directive_scratch_dir(opts);
                break;
            case KEYWORD_SECTOR:
                directive_sector(opts);
                break;
            case KEYWORD_MODEL:
                directive_model(opts);
                break;
            case KEYWORD_GOLDSTONE:
                opts->use_goldstone = 1;
                break;
            case KEYWORD_OCC:
                directive_occ(opts);
                break;
            case KEYWORD_OCC_IRREPS:
                directive_occ_irreps(opts);
                break;
            case KEYWORD_ACTIVE:
                directive_active(opts);
                break;
            case KEYWORD_NACTP:
                directive_nactp(opts);
                break;
            case KEYWORD_NACTH:
                directive_nacth(opts);
                break;
            case KEYWORD_DEGEN_THRESH:
                directive_degen_thresh(opts);
                break;
            case KEYWORD_HERMIT:
                opts->do_hermit = 1;
                break;
            case KEYWORD_MSTDM:
                opts->calc_model_space_tdms = 1;
                break;
            case KEYWORD_NATORB:
                directive_natorb(opts);
                break;
            case KEYWORD_NROOTS:
                directive_nroots(opts);
                break;
            case KEYWORD_ROOTS_CUTOFF:
                directive_roots_cutoff(opts);
                break;
            case KEYWORD_MAXITER:
                directive_maxiter(opts);
                break;
            case KEYWORD_CONV:
                directive_conv(opts);
                break;
            case KEYWORD_DIV_THRESH:
                directive_div_thresh(opts);
                break;
            case KEYWORD_DAMPING:
                directive_damping(opts);
                break;
            case KEYWORD_DIIS:
                directive_diis(opts);
                break;
            case KEYWORD_CROP:
                directive_crop(opts);
                break;
            case KEYWORD_SHIFT_TYPE:
                directive_shift_type(opts);
                break;
            case KEYWORD_SHIFT:
                directive_shift(opts);
                break;
            case KEYWORD_REUSE:
                directive_reuse(opts);
                break;
            case KEYWORD_SKIP:
                directive_skip(opts);
                break;
            case KEYWORD_FLUSH:
                directive_flush(opts);
                break;
            case KEYWORD_INTEGRALS:
                directive_integrals(opts);
                break;
            case KEYWORD_X2CMMF:
                opts->x2cmmf = 1;
                break;
            case KEYWORD_NEW_SORTING:
                opts->new_sorting = 1;
                break;
            case KEYWORD_GAUNT:
                directive_gaunt(opts);
                break;
            case KEYWORD_ONEPROP:
                directive_oneprop(opts);
                break;
            case KEYWORD_TWOPROP:
                directive_twoprop(opts);
                break;
            case KEYWORD_MEMORY:
                directive_memory(opts);
                break;
            case KEYWORD_TILESIZE:
                directive_tilesize(opts);
                break;
            case KEYWORD_DISK_USAGE:
                directive_disk_usage(opts);
                break;
            case KEYWORD_COMPRESS:
                opts->compress = CC_COMPRESS_LZ4;
                break;
            case KEYWORD_COMPRESS_TRIPLES:
                directive_compress_triples(opts);
                break;
            case KEYWORD_CUDA:
                opts->cuda_enabled = 1;
                break;
            case KEYWORD_OPENMP:
            case KEYWORD_NTHREADS:
                directive_nthreads(opts);
                break;
            case KEYWORD_OPENMP_ALGORITHM:
                directive_openmp_algorithm(opts);
                break;
            case KEYWORD_ARITH:
                directive_arith(opts);
                break;
            case KEYWORD_MDPROP:
                directive_mdprop(opts);
                break;
            case KEYWORD_TXTPROP:
                directive_txtprop(opts);
                break;
            case KEYWORD_ANALYT_PROP:
                directive_analyt_prop(opts);
                break;
            case KEYWORD_DENSITY:
                directive_density(opts);
                break;
            case KEYWORD_LAMBDA:
                directive_lambda(opts);
                break;
            case KEYWORD_OVERLAP:
                directive_overlap(opts);
                break;
            case KEYWORD_INTERFACE:
                directive_interface(opts);
                break;
            case KEYWORD_IH_IMMS:  // simple IH-like technique by A.V.Zaitsevskii
                directive_ih_incomplete_main_model_spaces(opts);
                break;
            case KEYWORD_RESTRICT_T3:
                directive_restrict_t3(opts);
                break;
            case KEYWORD_FLUSH_AMPLITUDES_TXT:
                opts->do_flush_amplitudes_txt = 1;
                break;
            case KEYWORD_HUGHES_KALDOR_1H2P:
                opts->hughes_kaldor_1h2p = 1;
                break;
            case KEYWORD_HUGHES_KALDOR_2H1P:
                opts->hughes_kaldor_2h1p = 1;
                break;
            case KEYWORD_USE_ORB_ENERGIES:
                opts->use_oe = 1;
                break;
            case KEYWORD_RECALC_ORB_ENERGIES:
                opts->use_oe = 0;
                break;
            case KEYWORD_SPINOR_LABELS:
                directive_spinor_labels(opts);
                break;
            case KEYWORD_TENSOR_TRAINS:
                directive_tensor_train(opts);
                break;
            case END_OF_LINE:
                // nothing to do
                break;
            default:
                yyerror("unknown keyword");
                break;
        }

        /* each "non-newline" directive must end with EOL or EOF */
        if (token_type != END_OF_LINE && token_type != END_OF_FILE) {
            token_type = next_token();
            if (token_type != END_OF_LINE && token_type != END_OF_FILE) {
                yyerror("end of line is expected");
            }
        }

        /* go to the next directive */
        token_type = next_token();
    }

    fclose(yyin);

    // post-processing of options

    // by default, if denominator shifts are specified by the "shift" directive,
    // but the "shift_type" directive is absent, the real simulation of imaginary
    // shifts will be switched on
    if (opts->shift_type == CC_SHIFT_NONE) {
        for (int h = 0; h < MAX_SECTOR_RANK; h++) {
            for (int p = 0; p < MAX_SECTOR_RANK; p++) {
                if (opts->shifts[h][p].enabled) {
                    opts->shift_type = CC_SHIFT_REAL;
                    opts->shifts[h][p].type = CC_SHIFT_REAL;
                }
            }
        }
    }

    // if occupation number were setted up explicitly => ignore occ_irreps input
    if (opts->occ_defined) {
        opts->nelec_defined = 0;
    }
    // natural (transition) orbitals will be calculated only for the target sector
    // if the sector was not explicitly specified
    for (int i = 0; i < opts->n_denmat; i++) {
        cc_denmat_query_t *query = &opts->denmat_query[i];
        query->sect1[0] = opts->sector_h;
        query->sect1[1] = opts->sector_p;
        query->sect2[0] = opts->sector_h;
        query->sect2[1] = opts->sector_p;
    }
    // switch on complex arithmetic if imaginary shifts are enabled
    if (opts->shift_type == CC_SHIFT_IMAG) {
        opts->recommended_arith = CC_ARITH_COMPLEX;
        printf("Recommended arithmetic is COMPLEX due to IMAGINARY shifts\n");
    }
    // switch on complex arithmetic if perturbation parameters
    // associated with property operators are complex-valued
    for (int i = 0; i < opts->n_oneprop; i++) {
        if (fabs(cimag(opts->oneprop_lambda[i])) > 1e-14) {
            opts->recommended_arith = CC_ARITH_COMPLEX;
            printf("Recommended arithmetic is COMPLEX due to IMAGINARY perturbation parameter\n");
        }
    }
    for (int i = 0; i < opts->n_twoprop; i++) {
        if (fabs(cimag(opts->twoprop_lambda[i])) > 1e-14) {
            opts->recommended_arith = CC_ARITH_COMPLEX;
            printf("Recommended arithmetic is COMPLEX due to IMAGINARY perturbation parameter\n");
        }
    }

    // if calculation of properties is required, but no information on density matrix is given,
    // switch on default value (lambda equations)
    if (opts->n_analyt_prop > 0 && opts->calc_density[0][0] == CC_DENSITY_MATRIX_DISABLED) {
        opts->calc_density[0][0] = CC_DENSITY_MATRIX_LAMBDA;
    }

    // only for tensor trains: remove the tilesize limit
    if (opts->tt_options.tt_module_enabled) {
        opts->tile_size = CC_MAX_SPINORS;
    }
}


/**
 * echo some file to some another file; returns 0 if OK, 1 if error
 */
int echo(char *file_name, FILE *dst)
{
    FILE *f;
    char line[MAX_LINE_LEN];

    f = fopen(file_name, "r");
    if (f == NULL) {
        return 1;
    }

    while (fgets(line, MAX_LINE_LEN, f) != NULL) {
        fprintf(dst, " %s", line);
    }

    return 0;
}


/*******************************************************************************
                       PARSER OF INPUT FILES - DIRECTIVES
 ******************************************************************************/

/**
 * Syntax:
 * title <quoted-string>
 */
void directive_title(cc_options_t *opts)
{
    static char *msg = "wrong argument!\n"
                       "A quoted string is expected";
    int bufsize = sizeof(opts->title);

    if (!match(TT_QUOTE)) {
        yyerror(msg);
    }
    // copy string ignoring quotes
    size_t len = strlen(yytext) - 2;
    strncpy(opts->title, yytext + 1, bufsize - 1);
    opts->title[bufsize - 1] = '\0';
    if (len < bufsize - 1) {
        opts->title[len] = '\0';
    }
}


/**
 * Syntax:
 * print ( low || medium || high || debug )
 */
void directive_print(cc_options_t *opts)
{
    static char *msg = "wrong print mode!\n"
                       "Possible values are: low, medium, high, debug\n"
                       "Additional printing options:\n"
                       "  \"model space\"\n"
                       "  \"eff config\"\n"
                       "  \"model vectors\"\n";

    int token_type = next_token();

    if (token_type == TT_WORD) {
        str_tolower(yytext);
        if (strcmp(yytext, "low") == 0) {
            opts->print_level = CC_PRINT_LOW;
        }
        else if (strcmp(yytext, "medium") == 0) {
            opts->print_level = CC_PRINT_MEDIUM;
        }
        else if (strcmp(yytext, "high") == 0) {
            opts->print_level = CC_PRINT_HIGH;
        }
        else if (strcmp(yytext, "debug") == 0) {
            opts->print_level = CC_PRINT_DEBUG;
        }
        else {
            yyerror(msg);
        }
    }
    else if (token_type == TT_QUOTE) {

        // to lowercase and remove quotation
        str_tolower(yytext);
        yytext = yytext + 1;
        yytext[strlen(yytext) - 1] = '\0';

        if (strcmp(yytext, "model space") == 0) {
            opts->print_model_space = 1;
        }
        else if (strcmp(yytext, "model vectors") == 0) {
            opts->print_model_vectors = 1;
        }
        else if (strcmp(yytext, "eff config") == 0) {
            opts->print_eff_config = 1;
        }
        else {
            yyerror(msg);
        }
    }
    else {
        yyerror(msg);
    }
}


/**
 * Syntax:
 * scratch_dir <string directory>]
 */
void directive_scratch_dir(cc_options_t *opts)
{
    static char *msg = "path not specified";
    int token_type;

    token_type = next_token();
    if (token_type == END_OF_LINE || token_type == END_OF_FILE) {
        yyerror(msg);
    }
    strcpy(opts->scratch_dir, yytext);
}


/**
 * Syntax:
 * sector <integer H>h<integer P>p
 */
void directive_sector(cc_options_t *opts)
{
    static char *msg = "wrong Fock space sector specification!";

    if (!match(TT_SECTOR)) {
        yyerror(msg);
    }

    parse_sector_label(yytext, &opts->sector_h, &opts->sector_p);
}


/**
 * Syntax:
 * one of: ( ccs || ccd || ccsd || ccsd(t) || ccsd+t(3) || ccsdt-1 || ccsdt-1' || ccsdt-2 || ccsdt-3 )
 */
void parse_cc_model(cc_options_t *opts)
{
    static char *msg = "wrong coupled cluster model!\n"
                       "Possible values: ccs, ccd, ccsd, ccsd+t(3)=ccsd(t), ccsdt-1a, ccsdt-1b=ccsdt-1, ccsdt-2, ccsdt-3";

    str_tolower(yytext);
    if (strcmp(yytext, "ccs") == 0) {
        opts->cc_model = CC_MODEL_CCS;
    }
    else if (strcmp(yytext, "ccd") == 0) {
        opts->cc_model = CC_MODEL_CCD;
    }
    else if (strcmp(yytext, "ccsd") == 0) {
        opts->cc_model = CC_MODEL_CCSD;
    }
    else if (strcmp(yytext, "ccsd+t(3)") == 0) {
        opts->cc_model = CC_MODEL_CCSD_T3;
    }
    else if (strcmp(yytext, "ccsd+t*(3)") == 0) {
        opts->cc_model = CC_MODEL_CCSD_T3_STAR;
    }
    else if (strcmp(yytext, "ccsd+t(4)") == 0) {
        opts->cc_model = CC_MODEL_CCSD_T4;
    }
    else if (strcmp(yytext, "ccsd+t*(4)") == 0) {
        opts->cc_model = CC_MODEL_CCSD_T4_STAR;
    }
    else if (strcmp(yytext, "ccsd(t)") == 0) {
        opts->cc_model = CC_MODEL_CCSD_T3;
    }
    else if (strcmp(yytext, "ccsdt-1a") == 0) {
        opts->cc_model = CC_MODEL_CCSDT_1A;
    }
    else if (strcmp(yytext, "ccsdt-1b") == 0) {
        opts->cc_model = CC_MODEL_CCSDT_1B;
    }
    else if (strcmp(yytext, "ccsdt-1") == 0) {
        opts->cc_model = CC_MODEL_CCSDT_1B;
    }
    else if (strcmp(yytext, "ccsdt-1'") == 0) {
        opts->cc_model = CC_MODEL_CCSDT_1B_PRIME;
    }
    else if (strcmp(yytext, "ccsdt-2") == 0) {
        opts->cc_model = CC_MODEL_CCSDT_2;
    }
    else if (strcmp(yytext, "ccsdt-3") == 0) {
        opts->cc_model = CC_MODEL_CCSDT_3;
    }
    else if (strcmp(yytext, "ccsdt") == 0) {
        opts->cc_model = CC_MODEL_CCSDT;
    }
    else {
        yyerror(msg);
    }
}


/**
 * Syntax:
 * occ <list of the 0 and 1 digits separated by spaces>
 */
void directive_occ(cc_options_t *opts)
{
    static char *msg = "wrong specification of occupation numbers!\n"
                       "A list of 0 or 1 separated by spaces is expected";
    int token_type;
    int count = 0;

    opts->occ_defined = 1;

    token_type = next_token();
    while (token_type != END_OF_LINE && token_type != END_OF_FILE) {
        if (token_type != TT_INTEGER) {
            yyerror(msg);
        }

        if (strcmp(yytext, "0") == 0) {
            opts->occ[count++] = 0;
        }
        else if (strcmp(yytext, "1") == 0) {
            opts->occ[count++] = 1;
        }
        else {
            yyerror(msg);
        }

        token_type = next_token();
    }
    put_back(token_type);
}


/**
 * Syntax:
 * occ_irreps <list of integers>
 */
void directive_occ_irreps(cc_options_t *opts)
{
    opts->nelec_defined = 1;
    read_space_specification(&opts->irrep_occ_numbers);
}


/**
 * Syntax:
 * active <real eps_min> <real eps_max>
 */
void directive_active(cc_options_t *opts)
{
    static char *msg1 = "wrong specification of active space!\n"
                        "Two real numbers are expected";
    static char *msg2 = "wrong specification of active space!\n"
                        "expected the 'energy' keyword or list of 1 or 0";
    static char *msg3 = "wrong specification of active space!\n"
                        "only 0 and 1 digits as flags are allowed";
    int token_type;

    token_type = next_token();
    if (token_type == TT_WORD) {
        str_tolower(yytext);
        if (strcmp(yytext, "energy") == 0) {
            // lower limit
            token_type = next_token();
            if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                yyerror(msg1);
            }
            opts->act_energy_min = atof(yytext);

            // upper limit
            token_type = next_token();
            if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                yyerror(msg1);
            }
            opts->act_energy_max = atof(yytext);

            opts->actsp_defined = CC_ACT_SPACE_SPEC_ENERGY;
            return;
        }
        else {
            yyerror(msg2);
        }
    }

    int i = 0;
    while (token_type != END_OF_LINE && token_type != END_OF_FILE && i < CC_MAX_SPINORS) {
        if (token_type != TT_INTEGER) {
            yyerror(msg3);
        }
        int flag = atoi(yytext);
        if (flag != 0 && flag != 1) {
            yyerror(msg3);
        }
        opts->active_each[i++] = flag;
        token_type = next_token();
    }
    opts->actsp_defined = CC_ACT_SPACE_SPEC_BINARY;
}


/**
 * Syntax:
 * nactp <integer dim>
 */
void directive_nactp(cc_options_t *opts)
{
    static char *msg = "wrong specification of the active space!\n"
                       "A positive integer or a list of pairs [rep]:nact is expected";
    int token_type;

    token_type = next_token();
    if (token_type == TT_INTEGER) {
        opts->nactp = atoi(yytext);
        opts->actsp_defined = CC_ACT_SPACE_SPEC_TOTAL;
        return;
    }
    else if (token_type == TT_ELEC_STATE) {
        do {
            char rep_name[32];
            int nact;
            int written = 0;

            parse_state_spec(yytext, rep_name, &nact);
            // try to find entry
            for (int i = 0; i < CC_MAX_NUM_IRREPS; i++) {
                cc_active_spec_t *sp = &opts->active_specs[i];
                if (strcmp(sp->rep_name, rep_name) == 0) {
                    sp->nactp = nact;
                    written = 1;
                    break;
                }
            }
            // if entry not found -- write to the first empty slot
            if (!written) {
                for (int i = 0; i < CC_MAX_NUM_IRREPS; i++) {
                    cc_active_spec_t *sp = &opts->active_specs[i];
                    if (*(sp->rep_name) == '\0') {
                        sp->nactp = nact;
                        strcpy(sp->rep_name, rep_name);
                        written = 1;
                        break;
                    }
                }
            }
            token_type = next_token();
        } while (token_type == TT_ELEC_STATE);
        put_back(token_type);
        opts->actsp_defined = CC_ACT_SPACE_SPEC_IRREPS;
    }
    else {
        yyerror(msg);
    }
}


/**
 * Syntax:
 * nacth <integer dim>
 */
void directive_nacth(cc_options_t *opts)
{
    static char *msg = "wrong specification of the active space!\n"
                       "A positive integer or a list of pairs [rep]:nact is expected";
    int token_type;

    token_type = next_token();
    if (token_type == TT_INTEGER) {
        opts->nacth = atoi(yytext);
        opts->actsp_defined = CC_ACT_SPACE_SPEC_TOTAL;
        return;
    }
    else if (token_type == TT_ELEC_STATE) {
        do {
            char rep_name[32];
            int nact;
            int written = 0;

            parse_state_spec(yytext, rep_name, &nact);
            // try to find entry
            for (int i = 0; i < CC_MAX_NUM_IRREPS; i++) {
                cc_active_spec_t *sp = &opts->active_specs[i];
                if (strcmp(sp->rep_name, rep_name) == 0) {
                    sp->nacth = nact;
                    written = 1;
                }
            }
            // if entry not found -- write to the first empty slot
            if (!written) {
                for (int i = 0; i < CC_MAX_NUM_IRREPS; i++) {
                    cc_active_spec_t *sp = &opts->active_specs[i];
                    if (*(sp->rep_name) == '\0') {
                        sp->nacth = nact;
                        strcpy(sp->rep_name, rep_name);
                        written = 1;
                    }
                }
            }
            token_type = next_token();
        } while (token_type == TT_ELEC_STATE);
        put_back(token_type);
        opts->actsp_defined = CC_ACT_SPACE_SPEC_IRREPS;
    }
    else {
        yyerror(msg);
    }
}


/**
 * Syntax:
 * degen_thresh <real thresh> [au | ev | cm | mhz]
 * (default units - au)
 */
void directive_degen_thresh(cc_options_t *opts)
{
    static char *msg = "wrong specification of degeneracy threshold!\n"
                       "A positive real is expected";
    double thresh;
    int token_type;

    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg);
    }

    thresh = atof(yytext);
    if (thresh <= 0) {
        yyerror(msg);
    }

    // try to read units of energy
    double units_conversion_factor = 1.0;

    token_type = next_token();
    if (token_type == TT_WORD) {
        str_tolower(yytext);
        if (strcmp(yytext, "au") == 0) {
            units_conversion_factor = 1.0;
        }
        else if (strcmp(yytext, "cm") == 0) {
            units_conversion_factor = 1.0 / 219474.63136320;
        }
        else if (strcmp(yytext, "ev") == 0) {
            units_conversion_factor = 1.0 / 27.211386245988;
        }
        else if (strcmp(yytext, "mhz") == 0) {
            units_conversion_factor = 1.0 / 6579683920.502;
        }
        else {
            yyerror("unknown units of energy, allowed units are:\n"
                    " au   atomic units (hartree)\n"
                    " ev   electronvolts\n"
                    " cm   wavenumbers\n"
                    " mhz  megahertz");
        }
    }
    else {
        put_back(token_type); // default units: atomic
    }

    // final value of a threshold - in atomic units (hartree)
    opts->degen_thresh = thresh * units_conversion_factor;
}


/**
 * Syntax:
 * natorb [<rep>:<state> | <rep1:state1> - <rep2:state2>] ...
 */
void directive_natorb(cc_options_t *opts)
{
    static char *msg = "wrong specification of electronic states!\n"
                       "A pair [rep_name]:state_no or a hyphen are expected";
    int token_type;

    token_type = next_token();
    while (token_type != END_OF_LINE && token_type != END_OF_FILE) {
        if (token_type != TT_ELEC_STATE) {
            yyerror(msg);
        }

        cc_denmat_query_t *query = &opts->denmat_query[opts->n_denmat];
        query->sect1[0] = 0;
        query->sect1[1] = 0;
        parse_state_spec(yytext, query->rep1_name, &query->state1);
        query->state1 -= 1;

        token_type = next_token();
        if (token_type == TT_HYPHEN) {
            token_type = next_token();
            if (token_type != TT_ELEC_STATE) {
                yyerror(msg);
            }
            query->sect2[0] = 0;
            query->sect2[1] = 0;
            parse_state_spec(yytext, query->rep2_name, &query->state2);
            query->state2 -= 1;
        }
        else {
            put_back(token_type);
            strcpy(query->rep2_name, query->rep1_name);
            query->state2 = query->state1;
        }
        opts->n_denmat++;

        token_type = next_token();
    }
    put_back(token_type);
}


/**
 * Syntax:
 * nroots <list-of-integers>
 */
void directive_nroots(cc_options_t *opts)
{
    static char *msg = "wrong specification of number of roots!\n"
                       "A list of pairs [irrep]:nroots separated by spaces is expected";
    int token_type;
    int count = 0;

    opts->nroots_specified = 1;
    read_space_specification(&opts->nroots_specs);

    if (opts->nroots_specs.deftype != CC_SPACE_BY_IRREPS) {
        yyerror(msg);
    }
}


/**
 * Syntax:
 * roots_cutoff <real cutoff_energy (in cm-1)>
 */
void directive_roots_cutoff(cc_options_t *opts)
{
    static char *msg = "wrong specification of the energy cutoff!\n"
                       "A real positive number is expected";
    static char *msg2 = "wrong specification of the energy cutoff!\n"
                        "Only atomic (au), inverse centimeters (cm) and electron-volt (eV) units are allowed";
    double cutoff_energy;
    double factor = 1.0;
    int token_type;

    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg);
    }

    cutoff_energy = atof(yytext);
    if (cutoff_energy <= 0) {
        yyerror(msg);
    }

    // units: cm or ev
    next_token();
    str_tolower(yytext);
    if (strcmp(yytext, "au") == 0) {
        factor = 1.0;
    }
    else if (strcmp(yytext, "cm") == 0) {  // cm^-1 => a.u.
        factor = 1 / CODATA_AU_TO_CM;
    }
    else if (strcmp(yytext, "ev") == 0) {  // eV => a.u.
        factor = 1 / CODATA_AU_TO_EV;
    }
    else {
        yyerror(msg2);
    }

    opts->roots_cutoff = cutoff_energy * factor;
    opts->roots_cutoff_specified = 1;
}


/**
 * Syntax:
 * maxiter <integer max>
 */
void directive_maxiter(cc_options_t *opts)
{
    static char *msg = "wrong specification of maximum number of iterations!\n"
                       "A positive integer is expected";

    if (!match(TT_INTEGER)) {
        yyerror(msg);
    }
    opts->maxiter = atoi(yytext);
}


/**
 * Syntax:
 * conv <real thresh>
 */
void directive_conv(cc_options_t *opts)
{
    static char *msg = "wrong specification of convergence threshold!\n"
                       "A positive real is expected";
    double thresh;
    int token_type;

    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg);
    }

    thresh = atof(yytext);
    if (thresh <= 0) {
        yyerror(msg);
    }
    opts->conv_thresh = thresh;
}


/**
 * Syntax:
 * div_thresh <real thresh>
 */
void directive_div_thresh(cc_options_t *opts)
{
    static char *msg = "wrong specification of the divergence threshold!\n"
                       "A positive real is expected";
    double thresh;
    int token_type;

    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg);
    }

    thresh = atof(yytext);
    if (thresh <= 0) {
        yyerror(msg);
    }
    opts->div_thresh = thresh;
}


/**
 * Syntax:
 * damping [<H>h<P>p] <integer last_step> <real factor>
 */
void directive_damping(cc_options_t *opts)
{
    static char *msg1 = "wrong specification of damping parameters!\n"
                        "A positive integer is expected\n";
    static char *msg2 = "wrong specification of damping parameters!\n"
                        "A non-negative real value is expected\n";
    int token_type;
    int sector_specified = 0;
    int h, p;
    int last_step;
    double factor;

    token_type = next_token();
    if (token_type == TT_SECTOR) {
        sector_specified = 1;
        h = yytext[0] - '0';
        p = yytext[2] - '0';
    }
    else {
        put_back(token_type);
    }

    // last step of damping
    token_type = next_token();
    if (token_type != TT_INTEGER) {
        yyerror(msg1);
    }
    last_step = atoi(yytext);

    // damping factor
    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg2);
    }
    factor = atof(yytext);

    // assign damping parameters
    if (!sector_specified) {
        for (int h = 0; h < MAX_SECTOR_RANK; h++) {
            for (int p = 0; p < MAX_SECTOR_RANK; p++) {
                opts->damping[h][p].enabled = 1;
                opts->damping[h][p].stop = last_step;
                opts->damping[h][p].factor = factor;
            }
        }
    }
    else {
        opts->damping[h][p].enabled = 1;
        opts->damping[h][p].stop = last_step;
        opts->damping[h][p].factor = factor;
    }
}


/**
 * Syntax:
 * diis [ off || <integer> || triples ]
 */
void directive_diis(cc_options_t *opts)
{
    static char *msg = "wrong specification of DIIS!\n"
                       "Possible arguments: no arguments, 'off' or integer value";
    static char *msg2 = "wrong specification of DIIS!\n"
                        "dimension of the DIIS subspace must be >= 2";

    int token_type = next_token();
    if (token_type == END_OF_LINE) {
        opts->diis_enabled = 1;
        opts->crop_enabled = 0;
        put_back(token_type);
        return;
    }
    else if (token_type == TT_WORD) {
        str_tolower(yytext);
        if (strcmp(yytext, "off") == 0) {
            opts->diis_enabled = 0;
        }
        else if (strcmp(yytext, "triples") == 0) {
            opts->diis_enabled = 1;
            opts->diis_triples = 1;
            opts->crop_enabled = 0;
        }
        else {
            yyerror(msg);
        }
    }
    else if (token_type == TT_INTEGER) {
        int dim = atoi(yytext);
        if (dim < 2) {
            yyerror(msg2);
        }
        opts->diis_enabled = 1;
        opts->crop_enabled = 0;
        opts->diis_dim = dim;
        return;
    }
    else {
        yyerror(msg);
    }
}


/**
 * Syntax:
 * crop [ off || <integer> ]
 */
void directive_crop(cc_options_t *opts)
{
    static char *msg = "wrong specification of CROP!\n"
                       "Possible arguments: no arguments, 'off' or integer value";
    static char *msg2 = "wrong specification of CROP!\n"
                        "dimension of the CROP subspace must be >= 2";

    int token_type = next_token();
    if (token_type == END_OF_LINE) {
        opts->crop_enabled = 1;
        opts->diis_enabled = 0;
        put_back(token_type);
        return;
    }
    else if (token_type == TT_WORD) {
        str_tolower(yytext);
        if (strcmp(yytext, "off") == 0) {
            opts->crop_enabled = 0;
        }
        else if (strcmp(yytext, "triples") == 0) {
            opts->crop_enabled = 1;
            opts->diis_enabled = 0;
            opts->crop_triples = 1;
        }
        else {
            yyerror(msg);
        }
    }
    else if (token_type == TT_INTEGER) {
        int dim = atoi(yytext);
        if (dim < 2) {
            yyerror(msg2);
        }
        opts->crop_enabled = 1;
        opts->diis_enabled = 0;
        opts->crop_dim = dim;
        return;
    }
    else {
        yyerror(msg);
    }
}


/**
 * Syntax:
 * shift_type ( none || real || realimag || imag || taylor )
 */
void directive_shift_type(cc_options_t *opts)
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

    opts->shift_type = type;
    for (int h = 0; h < MAX_SECTOR_RANK; h++) {
        for (int p = 0; p < MAX_SECTOR_RANK; p++) {
            opts->shifts[h][p].type = type;
        }
    }
}


/**
 * Syntax:
 * shift <H>h<P>p <int n> [<real shift_S1>] <real shift_S2> [<real shift_S3>]
 */
void directive_shift(cc_options_t *opts)
{
    static char *msg1 = "wrong specification of denominators shifts!\n"
                        "Label of the sector is expected";
    static char *msg2 = "wrong specification of denominators shifts!\n"
                        "An integer attenuation parameter is expected";
    static char *msg3 = "wrong specification of denominators shifts!\n"
                        "At least one real number is expected";
    int token_type;
    int h, p;
    int power;
    int count = 0;
    double s[CC_DIAGRAM_MAX_RANK];

    if (!match(TT_SECTOR)) {
        yyerror(msg1);
    }
    h = yytext[0] - '0';
    p = yytext[2] - '0';

    if (!match(TT_INTEGER)) {
        yyerror(msg2);
    }
    power = atoi(yytext);

    token_type = next_token();
    while (count < CC_DIAGRAM_MAX_RANK &&
           (token_type == TT_INTEGER || token_type == TT_FLOAT)) {

        s[count++] = atof(yytext);

        token_type = next_token();
    }
    put_back(token_type);

    if (count == 0) {
        yyerror(msg3);
    }

    // update option
    int offs = MAXIMUM2(h, p) - 1; // to skip amplitudes which are absent in this sector
    if (h == 0 && p == 0) {
        offs += 1;
    }
    opts->shifts[h][p].enabled = 1;
    opts->shifts[h][p].power = power;
    for (int i = 0; i <= count; i++) {
        opts->shifts[h][p].shifts[offs + i] = s[i];
    }
}


/**
 * Syntax:
 * reuse <list of arguments>
 * argument: one of: integrals 1-integrals 2-integrals amplitudes
 *                   0h0p 0h1p 1h0p 0h2p 2h0p 1h1p
 */
void directive_reuse(cc_options_t *opts)
{
    static char *msg = "wrong parameter of the reuse directive!\n"
                       "Possible values: integrals 1-integrals 2-integrals amplitudes 0h0p 0h1p 1h0p 0h2p 2h0p 1h1p 0h3p";
    int token_type;

    token_type = next_token();
    while (token_type != END_OF_LINE && token_type != END_OF_FILE) {
        if (token_type != TT_WORD && token_type != KEYWORD_INTEGRALS && token_type != TT_SECTOR) {
            yyerror(msg);
        }
        str_tolower(yytext);

        if (strcmp(yytext, "integrals") == 0) {
            opts->reuse_integrals_1 = 1;
            opts->reuse_integrals_2 = 1;
        }
        else if (strcmp(yytext, "1-integrals") == 0) {
            opts->reuse_integrals_1 = 1;
        }
        else if (strcmp(yytext, "2-integrals") == 0) {
            opts->reuse_integrals_2 = 1;
        }
        else if (strcmp(yytext, "amplitudes") == 0) {
            for (int h = 0; h < MAX_SECTOR_RANK; h++) {
                for (int p = 0; p < MAX_SECTOR_RANK; p++) {
                    opts->reuse_amplitudes[h][p] = 1;
                }
            }
        }
            // pattern: [nh]h[np]p
        else if (strlen(yytext) == 4 &&
                 isdigit(yytext[0]) && yytext[1] == 'h' &&
                 isdigit(yytext[2]) && yytext[3] == 'p') {

            int nh = yytext[0] - '0';
            int np = yytext[2] - '0';
            if (nh >= MAX_SECTOR_RANK || np >= MAX_SECTOR_RANK) {
                yyerror(msg);
            }

            opts->reuse_amplitudes[nh][np] = 1;
        }
        else {
            yyerror(msg);
        }

        token_type = next_token();
    }
    put_back(token_type);
}


/**
 * Syntax:
 * skip <list of FS sector symbols>
 *
 * Example:
 * skip 0h0p 0h1p
 */
void directive_skip(cc_options_t *opts)
{
    static char *msg = "wrong parameter of the 'skip' directive!\n"
                       "Fock space sector symbol is expected";
    int token_type;
    int h, p;

    token_type = next_token();
    while (token_type != END_OF_LINE && token_type != END_OF_FILE) {
        if (token_type != TT_SECTOR) {
            yyerror(msg);
        }
        str_tolower(yytext);

        parse_sector_label(yytext, &h, &p);
        opts->skip_sector[h][p] = 1;

        token_type = next_token();
    }
    put_back(token_type);
}


/**
 * Syntax:
 * flush <integer> ( iter || min || hrs )
 */
void directive_flush(cc_options_t *opts)
{
    static char *msg1 = "wrong specification of amplitude flushing\n"
                        "A positive integer number is expected";
    static char *msg2 = "wrong specification of amplitude flushing\n"
                        "Only number of iterations (iter), time in minutes (min) or hours (hrs) are allowed";
    int token_type;
    int ival;
    double factor;

    token_type = next_token();
    if (token_type != TT_INTEGER) {
        yyerror(msg1);
    }
    ival = atoi(yytext);

    next_token();
    str_tolower(yytext);
    if (strcmp(yytext, "iter") == 0) {
        opts->do_flush_iter = ival;
    }
    else if (strcmp(yytext, "min") == 0) {
    }
    else if (strcmp(yytext, "hrs") == 0) {
    }
    else {
        yyerror(msg2);
    }
}


/**
 * Syntax:
 * integrals <string "1-el int-s"> <string "2-el int-s"> [<string "prop-s int-s">]
 */
void directive_integrals(cc_options_t *opts)
{
    static char *msg = "at least one path must be specified";
    int token_type;

    token_type = next_token();
    if (token_type == END_OF_LINE || token_type == END_OF_FILE) {
        yyerror(msg);
    }
    strcpy(opts->integral_file_1, yytext);

    token_type = next_token();
    if (token_type == END_OF_LINE || token_type == END_OF_FILE) {
        put_back(token_type);
        return;
    }
    strcpy(opts->integral_file_2, yytext);

    token_type = next_token();
    if (token_type == END_OF_LINE || token_type == END_OF_FILE) {
        put_back(token_type);
        return;
    }
    strcpy(opts->integral_file_prop, yytext);
}


/**
 * Syntax:
 * gaunt [<string "gaunt int-s">]
 */
void directive_gaunt(cc_options_t *opts)
{
    static char *msg = "at least one path must be specified";
    int token_type;

    opts->gaunt_defined = 1;

    token_type = next_token();
    if (token_type == END_OF_LINE || token_type == END_OF_FILE) {
        put_back(token_type);
        return;
    }
    strcpy(opts->integral_file_gaunt, yytext);
}


/**
 * Syntax:
 * oneprop <real L_re> <real L_im> <string file_re> <string file_im>
 */
void directive_oneprop(cc_options_t *opts)
{
    static char *msg1 = "wrong parameter of the oneprop directive!\n"
                        "A real number is expected";
    static char *msg2 = "wrong parameter of the oneprop directive!\n"
                        "A path to the file is expected";
    int token_type;
    double lambda_re, lambda_im;

    // perturbation parameter
    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg1);
    }
    lambda_re = atof(yytext);

    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg1);
    }
    lambda_im = atof(yytext);

    // path to files
    token_type = next_token();
    if (token_type == END_OF_LINE || token_type == END_OF_FILE) {
        yyerror(msg2);
    }
    if (token_type == TT_QUOTE) {  // not file, but name of the property stored in the MDPROP file
        yytext[strlen(yytext) - 1] = '\0';  // remove quote
        strcpy(opts->mdprop_file[opts->n_mdprop], yytext + 1);
        opts->mdprop_lambda[opts->n_mdprop] = lambda_re + lambda_im * I;
        opts->n_mdprop++;
        return;
    }
    else {
        strcpy(opts->oneprop_file_re[opts->n_oneprop], yytext);
    }

    token_type = next_token();
    if (token_type == END_OF_LINE || token_type == END_OF_FILE) {
        yyerror(msg2);
    }
    strcpy(opts->oneprop_file_im[opts->n_oneprop], yytext);

    opts->oneprop_on = 1;
    opts->oneprop_lambda[opts->n_oneprop] = lambda_re + lambda_im * I;
    opts->n_oneprop++;
}


/**
 * Syntax:
 * twoprop <real L_re> <real L_im> <string file_name>
 */
void directive_twoprop(cc_options_t *opts)
{
    static char *msg1 = "wrong parameter of the twoprop directive!\n"
                        "A real number is expected";
    static char *msg2 = "wrong parameter of the twoprop directive!\n"
                        "A path to the file is expected";
    int token_type;
    double lambda_re, lambda_im;

    // perturbation parameter
    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg1);
    }
    lambda_re = atof(yytext);

    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg1);
    }
    lambda_im = atof(yytext);

    // path to files
    token_type = next_token();
    if (token_type == END_OF_LINE || token_type == END_OF_FILE) {
        yyerror(msg2);
    }
    else {
        strcpy(opts->twoprop_file[opts->n_twoprop], yytext);
    }

    opts->twoprop_on = 1;
    opts->twoprop_lambda[opts->n_twoprop] = lambda_re + lambda_im * I;
    opts->n_twoprop++;
}


/**
 * Syntax:
 * memory <real size> ( mb || gb )
 */
void directive_memory(cc_options_t *opts)
{
    static char *msg1 = "wrong specification of maximum memory usage!\n"
                        "A positive real number is expected";
    static char *msg2 = "wrong specification of maximum memory usage!\n"
                        "Minimum memory size to be allocated must be >= 10 Mb";
    static char *msg3 = "wrong specification of maximum memory usage!\n"
                        "Only megabytes (mb) and gigabytes (gb) units are allowed";
    int token_type;
    double double_val;
    double factor;
    double const MEM_THRESH = 10 * 1024 * 1024;  // 10 Mb

    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg1);
    }
    double_val = atof(yytext);

    next_token();
    str_tolower(yytext);
    if (strcmp(yytext, "mb") == 0) {
        factor = 1024 * 1024;
    }
    else if (strcmp(yytext, "gb") == 0) {
        factor = 1024 * 1024 * 1024;
    }
    else {
        yyerror(msg3);
    }
    double_val = double_val * factor;

    if (double_val < MEM_THRESH) {
        yyerror(msg2);
    }

    opts->max_memory_size = (size_t) double_val;
}


/**
 * Syntax:
 * tilesize <integer size>
 */
void directive_tilesize(cc_options_t *opts)
{
    static char *msg = "wrong specification of tilesize!\n"
                       "A positive integer is expected";

    if (!match(TT_INTEGER)) {
        yyerror(msg);
    }
    opts->tile_size = atoi(yytext);
}


/**
 * Syntax:
 * disk_usage <integer mode>
 */
void directive_disk_usage(cc_options_t *opts)
{
    static char *msg = "wrong specification of disk usage!\n"
                       "A positive integer is expected";
    int level;

    if (!match(TT_INTEGER)) {
        yyerror(msg);
    }
    level = atoi(yytext);

    if (level > 4) {
        opts->disk_usage_level = 4;
    }
    else {
        opts->disk_usage_level = level;
    }

    if (opts->disk_usage_level == 4) {
        opts->compress = CC_COMPRESS_LZ4;
    }
}


/**
 * Syntax:
 * ( nthreads || openmp ) <integer n_omp_threads>
 */
void directive_nthreads(cc_options_t *opts)
{
    static char *msg = "wrong specification of number of OpenMP threads!\n"
                       "A positive integer is expected";

    if (!match(TT_INTEGER)) {
        yyerror(msg);
    }
    opts->nthreads = atoi(yytext);
}


/**
 * Syntax:
 * openmp_algorithm ( internal || external )
 */
void directive_openmp_algorithm(cc_options_t *opts)
{
    static char *msg = "wrong specification of parallelization algorithm!\n"
                       "Possible values: internal, external";

    if (!match(TT_WORD)) {
        yyerror(msg);
    }
    str_tolower(yytext);
    if (strcmp(yytext, "internal") == 0) {
        opts->openmp_algorithm = CC_OPENMP_ALGORITHM_INTERNAL;
    }
    else if (strcmp(yytext, "external") == 0) {
        opts->openmp_algorithm = CC_OPENMP_ALGORITHM_EXTERNAL;
    }
    else {
        yyerror(msg);
    }
}


/**
 * Syntax:
 * arith ( complex || real )
 */
void directive_arith(cc_options_t *opts)
{
    static char *msg = "wrong specification of arithmetic!\n"
                       "Possible values: real, complex";

    if (!match(TT_WORD)) {
        yyerror(msg);
    }
    str_tolower(yytext);
    if (strcmp(yytext, "real") == 0) {
        opts->recommended_arith = CC_ARITH_REAL;
    }
    else if (strcmp(yytext, "complex") == 0) {
        opts->recommended_arith = CC_ARITH_COMPLEX;
    }
    else {
        yyerror(msg);
    }
}


/**
 * Syntax:
 * mdprop "<property-name>" [transpose]
 * Example:
 * mdprop "ZDIPLEN"
 */
void directive_mdprop(cc_options_t *opts)
{
    static char *msg_param = "wrong parameter of the mdprop directive";
    static char *msg_trans = "the 'transpose' keyword is expected";
    int token_type;

    // quotation mark + property_name + quotation mark
    token_type = next_token();
    if (token_type != TT_QUOTE) {
        yyerror(msg_param);
    }
    yytext[strlen(yytext) - 1] = '\0';  // remove quote
    cc_ms_prop_query_t *q = &opts->prop_queries[opts->n_model_space_props];
    strcpy(q->prop_name, yytext + 1);
    q->do_transpose = 0;
    q->scheme = CC_DIRECT_PROP_SCHEME_HERMITIAN;
    q->source = CC_PROP_FROM_MDPROP;
    strcpy(q->irrep_name, "");
    q->approx_denominator = CC_PROPERTIES_APPROX_MODEL_SPACE;
    q->approx_numerator = CC_PROPERTIES_APPROX_MODEL_SPACE;

    while (1) {
        token_type = next_token();

        if (token_type == TT_WORD) {

            // optional keywords
            if (strcmp(yytext, "transpose") == 0) {
                q->do_transpose = 1;
            }
            else if (strcmp(yytext, "sym") == 0) {
                token_type = next_token();
                if (!(token_type == TT_WORD || token_type == TT_INTEGER)) {
                    yyerror("expected irrep name");
                }
                strcpy(q->irrep_name, yytext);
            }
            else if (strcmp(yytext, "scheme") == 0) {
                token_type = next_token();

                if (token_type == TT_WORD && strcmp(yytext, "hermitian") == 0) {
                    q->scheme = CC_DIRECT_PROP_SCHEME_HERMITIAN;
                }
                else if (token_type == TT_WORD && strcmp(yytext, "non-hermitian") == 0) {
                    q->scheme = CC_DIRECT_PROP_SCHEME_NON_HERMITIAN;
                }
                else if (token_type == TT_WORD && strcmp(yytext, "connected") == 0) {
                    q->scheme = CC_DIRECT_PROP_SCHEME_CONNECTED;
                }
                else {
                    yyerror("expected keyword: 'hermitian', 'non-hermitian' or 'connected'");
                }
            }
            else if (strcmp(yytext, "approx") == 0) {
                /*
                 * get order of numerator
                 */
                token_type = next_token();
                if (token_type != TT_INTEGER) {
                    yyerror("integer expected (order of numerator)");
                }
                int numer_order = atoi(yytext);
                if (!(numer_order == 0 || numer_order == 1 || numer_order == 2)) {
                    yyerror("allowed orders of numerator: 0, 1, 2");
                }
                if (numer_order == 0) {
                    numer_order = CC_PROPERTIES_APPROX_MODEL_SPACE;
                }
                if (numer_order == 1) {
                    numer_order = CC_PROPERTIES_APPROX_LINEAR;
                }
                if (numer_order == 2) {
                    numer_order = CC_PROPERTIES_APPROX_QUADRATIC;
                }
                q->approx_numerator = numer_order;

                /*
                 * get order of denominator
                 */
                token_type = next_token();
                if (token_type != TT_INTEGER) {
                    yyerror("integer expected (order of denominator)");
                }
                int denom_order = atoi(yytext);
                if (!(denom_order == 0 || denom_order == 2)) {
                    yyerror("allowed orders of denominator: 0, 2");
                }
                if (denom_order == 0) {
                    denom_order = CC_PROPERTIES_APPROX_MODEL_SPACE;
                }
                if (denom_order == 2) {
                    denom_order = CC_PROPERTIES_APPROX_QUADRATIC;
                }
                q->approx_denominator = denom_order;
            }
            else {
                yyerror(msg_trans);
            }
        }
        else {
            put_back(token_type);
            break;
        }
    }

    opts->n_model_space_props++;
}


/**
 * Syntax:
 * txtprop <file_real_part> <file_imag_part> [transpose]
 * Example:
 * txtprop PropInts_Re.txt PropInts_Im.txt
 */
void directive_txtprop(cc_options_t *opts)
{
    static char *msg_param = "wrong parameter of the txtprop directive";
    static char *msg_trans = "the 'transpose' keyword is expected";
    int token_type;

    cc_ms_prop_query_t *q = &opts->prop_queries[opts->n_model_space_props];
    q->source = CC_PROP_FROM_TXTPROP;
    q->do_transpose = 0;
    q->scheme = CC_DIRECT_PROP_SCHEME_HERMITIAN;
    strcpy(q->irrep_name, "");
    q->approx_denominator = CC_PROPERTIES_APPROX_MODEL_SPACE;
    q->approx_numerator = CC_PROPERTIES_APPROX_MODEL_SPACE;


    // path to file with real part
    token_type = next_token();
    if (token_type != TT_WORD) {
        yyerror(msg_param);
    }
    strcpy(q->file_real, yytext);

    // path to file with imag part
    token_type = next_token();
    if (token_type != TT_WORD) {
        yyerror(msg_param);
    }
    strcpy(q->file_imag, yytext);


    while (1) {
        token_type = next_token();

        if (token_type == TT_WORD) {

            // optional keywords
            if (strcmp(yytext, "transpose") == 0) {
                q->do_transpose = 1;
            }
            else if (strcmp(yytext, "sym") == 0) {
                token_type = next_token();
                if (!(token_type == TT_WORD || token_type == TT_INTEGER)) {
                    yyerror("expected irrep name");
                }
                strcpy(q->irrep_name, yytext);
            }
            else if (strcmp(yytext, "scheme") == 0) {
                token_type = next_token();

                if (token_type == TT_WORD && strcmp(yytext, "hermitian") == 0) {
                    q->scheme = CC_DIRECT_PROP_SCHEME_HERMITIAN;
                }
                else if (token_type == TT_WORD && strcmp(yytext, "non-hermitian") == 0) {
                    q->scheme = CC_DIRECT_PROP_SCHEME_NON_HERMITIAN;
                }
                else if (token_type == TT_WORD && strcmp(yytext, "connected") == 0) {
                    q->scheme = CC_DIRECT_PROP_SCHEME_CONNECTED;
                }
                else {
                    yyerror("expected keyword: 'hermitian', 'non-hermitian' or connected");
                }
            }
            else if (strcmp(yytext, "approx") == 0) {
                /*
                 * get order of numerator
                 */
                token_type = next_token();
                if (token_type != TT_INTEGER) {
                    yyerror("integer expected (order of numerator)");
                }
                int numer_order = atoi(yytext);
                if (!(numer_order == 0 || numer_order == 1 || numer_order == 2)) {
                    yyerror("allowed orders of numerator: 0, 1, 2");
                }
                if (numer_order == 0) {
                    numer_order = CC_PROPERTIES_APPROX_MODEL_SPACE;
                }
                if (numer_order == 1) {
                    numer_order = CC_PROPERTIES_APPROX_LINEAR;
                }
                if (numer_order == 2) {
                    numer_order = CC_PROPERTIES_APPROX_QUADRATIC;
                }
                q->approx_numerator = numer_order;

                /*
                 * get order of denominator
                 */
                token_type = next_token();
                if (token_type != TT_INTEGER) {
                    yyerror("integer expected (order of denominator)");
                }
                int denom_order = atoi(yytext);
                if (!(denom_order == 0 || denom_order == 2)) {
                    yyerror("allowed orders of denominator: 0, 2");
                }
                if (denom_order == 0) {
                    denom_order = CC_PROPERTIES_APPROX_MODEL_SPACE;
                }
                if (denom_order == 2) {
                    denom_order = CC_PROPERTIES_APPROX_QUADRATIC;
                }
                q->approx_denominator = denom_order;
            }
            else {
                yyerror(msg_trans);
            }
        }
        else {
            put_back(token_type);
            break;
        }
    }


    opts->n_model_space_props++;
}


/**
 * Syntax:
 * analyt_prop <property-name>
 * Example:
 * analyt_prop ZDIPLEN
 */
void directive_analyt_prop(cc_options_t *opts)
{
    static char *msg_param = "wrong parameter of the analyt_prop directive";

    int token_type = next_token();
    if (token_type != TT_WORD) {
        yyerror(msg_param);
    }

    strcpy(opts->analyt_prop_files[opts->n_analyt_prop], yytext);
    opts->n_analyt_prop++;
}


/**
 * Syntax:
 * density <sector-label> [lambda|expect] [irrep:state]

 * Example:
 * density 0h0p
 * (lambda by default)
 * or
 * density 0h0p expect
 * or
 * density 0h1p 1/2+:1
 * (the 1st state of symmetry 1/2+)
 */
void directive_density(cc_options_t *opts)
{
    static char *msg_param = "wrong parameter of the 'density' directive";
    static char *msg_wrong_sector = "density matrix calculation is not implemented for this sector";

    int token_type = next_token();
    if (token_type != TT_SECTOR) {
        yyerror(msg_param);
    }

    /*
     * parse sector label
     */
    int sect_h = 0, sect_p = 0;
    parse_sector_label(yytext, &sect_h, &sect_p);
    opts->selects[opts->n_select].sect_h = sect_h;
    opts->selects[opts->n_select].sect_p = sect_p;
    if ((sect_h == 0 && sect_p == 0) ||
        (sect_h == 0 && sect_p == 1) ||
        (sect_h == 1 && sect_p == 0) ||
        (sect_h == 0 && sect_p == 2) ||
        (sect_h == 0 && sect_p == 3)) {
        // OK
    }
    else {
        yyerror(msg_wrong_sector);
    }

    /*
     * defaults
     */
    opts->calc_density[sect_h][sect_p] = CC_DENSITY_MATRIX_LAMBDA;

    /*
     * optional method of calculation: "lambda" or "expect" (0h0p only)
     */
    if (sect_h == 0 && sect_p == 0) {
        token_type = next_token();
        if (token_type == TT_WORD) {
            str_tolower(yytext);
            if (strcmp(yytext, "lambda") == 0) {
                    //opts->calc_density_0h0p = CC_DENSITY_MATRIX_LAMBDA;
                opts->calc_density[0][0] = CC_DENSITY_MATRIX_LAMBDA;
            }
            else if (strcmp(yytext, "expect") == 0) {
                    //opts->calc_density_0h0p = CC_DENSITY_MATRIX_EXPECTATION;
                opts->calc_density[0][0] = CC_DENSITY_MATRIX_EXPECTATION;
            }
            else {
                yyerror(msg_param);
            }
        }
        else {
            put_back(token_type);
        }
        return;
    }

    /*
     * number of electronic state of interest
     * mandatory for the 0h1p, 1h0p, 0h2p and 0h3p sectors
     */
    token_type = next_token();

    if (token_type == TT_WORD && strcmp(yytext, "all") == 0) {
        opts->density_num_states[sect_h][sect_p] = -1;
        return;
    }

    if (token_type != TT_ELEC_STATE) {
        yyerror("symbol of electronic state is required");
    }

    cc_denmat_query_t *query = opts->density_target_states[sect_h][sect_p] + opts->density_num_states[sect_h][sect_p];
    query->sect1[0] = sect_h;
    query->sect1[1] = sect_p;
    parse_state_spec(yytext, query->rep1_name, &query->state1);
    query->state1 -= 1;
    opts->density_num_states[sect_h][sect_p] += 1;

    // try to read the second state, maybe the transition density matrix is specified

    token_type = next_token();
    if (token_type == TT_HYPHEN) {
        token_type = next_token();
        if (token_type != TT_ELEC_STATE) {
            yyerror("symbol of electronic state is required");
        }
        query->sect2[0] = sect_h;
        query->sect2[1] = sect_p;
        parse_state_spec(yytext, query->rep2_name, &query->state2);
        query->state2 -= 1;
    }
    else {
        put_back(token_type);
        strcpy(query->rep2_name, query->rep1_name);
        query->state2 = query->state1;
    }
}


/**
 * Syntax:
 * lambda <sector-label> irrep:state

 * Example:
 * lambda 0h1p 1/2+:1
 * (the 1st state of symmetry 1/2+)
 */
void directive_lambda(cc_options_t *opts)
{
    static char *msg_param = "wrong parameter of the 'lambda' directive";
    static char *msg_wrong_sector = "lambda equations are not implemented for this sector";

    int token_type = next_token();
    if (token_type != TT_SECTOR) {
        yyerror(msg_param);
    }

    /*
     * parse sector label
     */
    int sect_h = 0, sect_p = 0;
    parse_sector_label(yytext, &sect_h, &sect_p);
    opts->selects[opts->n_select].sect_h = sect_h;
    opts->selects[opts->n_select].sect_p = sect_p;
    if (sect_h == 0 && sect_p == 1) {
        // OK
    }
    else {
        yyerror(msg_wrong_sector);
    }

    /*
     * defaults
     */
    if (sect_h == 0 && sect_p == 1) {
        opts->calc_lambda_0h1p = 1;
    }

    /*
     * number of electronic state of interest
     * mandatory for the 0h1p, 0h2p and 0h3p sectors
     */

    // 0h1p
    if (sect_h == 0 && sect_p == 1) {
        token_type = next_token();

        if (token_type != TT_ELEC_STATE) {
            yyerror("symbol of electronic state is required");
        }

        cc_denmat_query_t *query = opts->lambda_0h1p_states + opts->lambda_0h1p_num_states;
        query->sect1[0] = 0;
        query->sect1[1] = 1;
        parse_state_spec(yytext, query->rep1_name, &query->state1);
        query->state1 -= 1;
        opts->lambda_0h1p_num_states++;
    }
}


/**
 * Syntax:
 * overlap <sector-label>
 *
 * Example:
 * overlap 0h1p 0h2p
 */
void directive_overlap(cc_options_t *opts)
{
    static char *msg_param = "wrong parameter of the 'overlap' directive";
    static char *msg_wrong_sector = "cannot calculate overlap in this sector (not implemented)";

    int n_sectors = 0;
    int token_type = next_token();

    while (token_type != END_OF_LINE && token_type != END_OF_FILE) {
        if (token_type != TT_SECTOR) {
            yyerror(msg_param);
        }

        int sect_h = 0, sect_p = 0;
        str_tolower(yytext);
        parse_sector_label(yytext, &sect_h, &sect_p);

        // check if the code for overlap integrals is implemented for the required FS sector
        if (sect_h == 0 && sect_p == 0) {
            yyerror("use the 'density' keyword to calculate overlap in the 0h0p sector");
        }
        if ((sect_h == 1 && sect_p == 0) ||
            (sect_h == 0 && sect_p == 1) ||
            (sect_h == 1 && sect_p == 1) ||
            (sect_h == 2 && sect_p == 0) ||
            (sect_h == 0 && sect_p == 2) ||
            (sect_h == 3 && sect_p == 0) ||
            (sect_h == 0 && sect_p == 3)) {
            // OK
        }
        else {
            yyerror(msg_wrong_sector);
        }

        n_sectors++;
        opts->calc_overlap[sect_h][sect_p] = 1;

        token_type = next_token();
    }
    put_back(token_type);

    if (n_sectors == 0) {
        yyerror("at least one sector is required");
    }
}


/**
 * Syntax:
 * interface (dirac | pyscf)
 */
void directive_interface(cc_options_t *opts)
{
    static char *msg = "unknown interface";
    int token_type;

    token_type = next_token();
    if (token_type == TT_WORD) {
        str_tolower(yytext);
        if (strcmp(yytext, "dirac") == 0) {
            cc_opts->int_source = CC_INTEGRALS_DIRAC;
        }
        else if (strcmp(yytext, "pyscf") == 0) {
            cc_opts->int_source = CC_INTEGRALS_PYSCF;
        }
        else {
            yyerror(msg);
        }
    }
    else {
        yyerror(msg);
    }
}


void read_space_specification(cc_space_t *space)
{
    static char *msg = "wrong specification of the space";
    int token_type;

    token_type = next_token();
    if (token_type == TT_WORD) {
        str_tolower(yytext);
        if (strcmp(yytext, "energy") != 0) {
            yyerror(msg);
        }
        token_type = next_token();
        if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
            yyerror(msg);
        }
        space->emin = atof(yytext);
        token_type = next_token();
        if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
            yyerror(msg);
        }
        space->emax = atof(yytext);
        space->deftype = CC_SPACE_BY_ENERGY;
    }
    else if (token_type == TT_INTEGER) {
        space->deftype = CC_SPACE_BY_TOTAL;
        space->total = atoi(yytext);
        return;
    }
    else if (token_type == TT_ELEC_STATE) {
        space->deftype = CC_SPACE_BY_IRREPS;
        do {
            char rep_name[32];
            int norb;
            int written = 0;

            parse_state_spec(yytext, rep_name, &norb);

            // try to find entry
            for (int i = 0; i < CC_MAX_NUM_IRREPS; i++) {
                if (strcmp(space->rep_names[i], rep_name) == 0) {
                    space->dim[i] = norb;
                    written = 1;
                    break;
                }
            }
            // if entry not found -- write to the first empty slot
            if (!written) {
                for (int i = 0; i < CC_MAX_NUM_IRREPS; i++) {
                    if (*(space->rep_names[i]) == '\0') {
                        space->dim[i] = norb;
                        strcpy(space->rep_names[i], rep_name);
                        written = 1;
                        break;
                    }
                }
            }
            token_type = next_token();
        } while (token_type == TT_ELEC_STATE);
        put_back(token_type);
    }
    else {
        yyerror(msg);
    }
}


/**
 * Syntax (example):
 * model
 *   ccsdt
 *   0h0p t3 act2act != 0
 *   0h1p t3 eps_window -1 +1 != 0
 *   ...
 * end
 *
 * (or)
 *
 * model (ccsd||ccsdt||...)
 */
void directive_model(cc_options_t *opts)
{
    int token_type;
    int sect_h, sect_p;
    int spaces_defined = 0;

    // "end of line" or "ccsd... etc" after the opening keyword
    token_type = next_token();
    if (token_type != END_OF_LINE) {
        parse_cc_model(opts);
        return;
    }

    while (1) {

        // sector of the cluster operator
        token_type = next_token();
        if (token_type == KEYWORD_END) {
            return;
        }
        if (token_type == TT_WORD) {
            parse_cc_model(opts);
            token_type = next_token();
            if (token_type != END_OF_LINE) {
                yyerror("end of line is expected");
            }
            continue;
        }

        if (token_type != TT_SECTOR) {
            yyerror("Fock space sector symbol is expected");
        }
        parse_sector_label(yytext, &sect_h, &sect_p);
        opts->selects[opts->n_select].sect_h = sect_h;
        opts->selects[opts->n_select].sect_p = sect_p;

        // excitation level of the cluster operator
        token_type = next_token();
        if (token_type == TT_WORD && strcmp(yytext, "t1") == 0) {
            opts->selects[opts->n_select].rank = 2;
        }
        else if (token_type == TT_WORD && strcmp(yytext, "t2") == 0) {
            opts->selects[opts->n_select].rank = 4;
        }
        else if (token_type == TT_WORD && strcmp(yytext, "t3") == 0) {
            opts->selects[opts->n_select].rank = 6;
        }
        else {
            yyerror("tag [t1|t2|t3] is expected");
        }

        // selection rule
        token_type = next_token();
        // short form
        if (token_type == TT_EQ || token_type == TT_NEQ) {
            opts->selects[opts->n_select].task =
                (token_type == TT_EQ) ? CC_SELECTION_SET_ZERO : CC_SELECTION_SET_ZERO_EXCEPT;
            token_type = next_token();
            if (token_type != TT_INTEGER || atoi(yytext) != 0) {
                yyerror("'0' (zero) is expected");
            }
            opts->selects[opts->n_select].rule = CC_SELECTION_ALL;
            opts->n_select++;
            // each subdirective must end with end of line!
            token_type = next_token();
            if (token_type != END_OF_LINE) {
                yyerror("end of line is expected");
            }
            continue;
        }
            // full form
        else if (token_type == TT_WORD && strcmp(yytext, "all") == 0) {
            opts->selects[opts->n_select].rule = CC_SELECTION_ALL;
        }
        else if (token_type == TT_WORD && strcmp(yytext, "spectator") == 0) {
            opts->selects[opts->n_select].rule = CC_SELECTION_SPECTATOR;
        }
        else if (token_type == TT_WORD && strcmp(yytext, "act_to_act") == 0) {
            opts->selects[opts->n_select].rule = CC_SELECTION_ACT_TO_ACT;
        }
        else if (token_type == TT_WORD && strcmp(yytext, "max2inact") == 0) {
            opts->selects[opts->n_select].rule = CC_SELECTION_MAX_2_INACT;
        }
        else if (token_type == TT_WORD &&
                 (strcmp(yytext, "exc_window") == 0 ||
                  strcmp(yytext, "eps_window") == 0)) {
            if (strcmp(yytext, "exc_window") == 0) {
                opts->selects[opts->n_select].rule = CC_SELECTION_EXC_WINDOW;
            }
            else {
                opts->selects[opts->n_select].rule = CC_SELECTION_EPS_WINDOW;
            }
            token_type = next_token();
            if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                yyerror("float number is expected");
            }
            opts->selects[opts->n_select].e1 = atof(yytext);
            token_type = next_token();
            if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                yyerror("float number is expected");
            }
            opts->selects[opts->n_select].e2 = atof(yytext);
        }
        else {
            yyerror("keyword [all|spectator|act_to_act|exc_window|eps_window] is expected");
        }

        // what to do: set to zero or not?
        token_type = next_token();
        if (token_type == TT_EQ) {
            opts->selects[opts->n_select].task = CC_SELECTION_SET_ZERO;
            token_type = next_token();
            if (token_type != TT_INTEGER || atoi(yytext) != 0) {
                yyerror("'0' (zero) is expected");
            }
        }
        else if (token_type == TT_NEQ) {
            opts->selects[opts->n_select].task = CC_SELECTION_SET_ZERO_EXCEPT;
            token_type = next_token();
            if (token_type != TT_INTEGER || atoi(yytext) != 0) {
                yyerror("'0' (zero) is expected");
            }
        }
        else {
            yyerror("'=' or '!=' sign is expected");
        }

        opts->n_select++;

        // each subdirective must end with end of line!
        token_type = next_token();
        if (token_type != END_OF_LINE) {
            yyerror("end of line is expected");
        }
    }
}


/**
 * Syntax:
 * restrict_t3 <double lower_bound> <double upper_bound>
 */
void directive_restrict_t3(cc_options_t *opts)
{
    static char *msg = "wrong specification of spinor space for triples!\n"
                       "Two reals are expected";
    double lower_bound;
    double upper_bound;
    int token_type;

    // lower bound for spinor energies
    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg);
    }
    lower_bound = atof(yytext);

    // upper bound for spinor energies
    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg);
    }
    upper_bound = atof(yytext);

    opts->do_restrict_t3 = 1;
    opts->restrict_t3_bounds[0] = lower_bound;
    opts->restrict_t3_bounds[1] = upper_bound;
}


void directive_compress_triples(cc_options_t *opts)
{
    static char *msg1 = "wrong specification of the triples compression parameters!\n"
                        "real number is expected";
    static char *msg2 = "wrong specification of the triples compression parameters!\n"
                        "allowed data types are 'float' or 'double'";
    int token_type;
    int data_type;
    double zero_thresh;

    // threshold
    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg1);
    }
    zero_thresh = atof(yytext);

    // type: float or double
    next_token();
    str_tolower(yytext);
    if (strcmp(yytext, "float") == 0) {
        data_type = CC_FLOAT;
    }
    else if (strcmp(yytext, "double") == 0) {
        data_type = CC_DOUBLE;
    }
    else {
        yyerror(msg2);
    }

    opts->do_compress_triples = 1;
    opts->compress_triples_thresh = zero_thresh;
    opts->compress_triples_data_type = data_type;
}


void directive_spinor_labels(cc_options_t *opts)
{
    static char *err_end_of_line = "end of line is expected";
    static char *err_keyword = "unknown keyword";
    static char *err_integer = "integer number is expected";
    static char *err_string = "label is expected";

    // end of line after the opening keyword is required
    int token_type = next_token();
    if (token_type != END_OF_LINE) {
        yyerror(err_end_of_line);
    }

    // loop over sub-directives
    token_type = next_token();
    while (token_type != END_OF_FILE) {

        if (token_type == TT_INTEGER) {
            int spinor_number = atoi(yytext);
            if (spinor_number <= 0 || spinor_number > CC_MAX_SPINORS) {
                errquit("spinor number must be in range [1,%d]", CC_MAX_SPINORS);
            }

            token_type = next_token();
            if (token_type != END_OF_LINE && token_type != END_OF_FILE) {
                opts->spinor_labels[spinor_number - 1] = cc_strdup(yytext);
            }
            else {
                errquit(err_string);
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

        /* go to the next directive */
        token_type = next_token();
    }
}


/*******************************************************************************
                     PARSER OF INPUT FILES - GLOBAL VARIABLES
 ******************************************************************************/

// flag: was the last token returned to the lexer?
static int pushed_back = 0;

// name of the input file with language to be processed
static char input_file_name[CC_MAX_PATH_LENGTH];


/*******************************************************************************
                     PARSER OF INPUT FILES - HELPER FUNCTIONS
 ******************************************************************************/


void set_input_file_name(char *file_name)
{
    strcpy(input_file_name, file_name);
}


void print_line_by_num(char *file_name, int required_lineno)
{
    FILE *f;
    int c;
    int lineno;

    f = fopen(file_name, "r");
    if (f == NULL) {
        errquit("print_line_by_num(): cannot open file '%s'", file_name);
    }

    lineno = 1;
    while ((c = fgetc(f)) != EOF) {
        if (lineno == required_lineno) {
            printf(" ");
            do {
                putchar(c);
                c = fgetc(f);
            } while (c != '\n' && c != EOF);
            printf("\n");
            fclose(f);
            return;
        }
        if (c == '\n') {
            lineno++;
        }
    }

    errquit("print_line_by_num(): no line number %d in '%s'", required_lineno, file_name);
    fclose(f);
}


void yyerror(char *s)
{
    int lineno = yylineno;

    // lineno = yylineno-1 if last token was END_OF_LINE
    if (strcmp(yytext, "\n") == 0) {
        lineno = yylineno - 1;
    }
    fclose(yyin);

    printf("%s:%d:%d: error: %s\n\n", input_file_name, lineno, yycol, s);
    print_line_by_num(input_file_name, lineno);

    // print pointer to the erroneous token
    printf(" ");
    for (int i = 0; i < yycol; i++) {
        printf(" ");
    }
    printf("^\n\n");
    printf("execution terminated.\n");

    exit(1);
}


int next_token()
{
    if (pushed_back) {
        int token_type = pushed_back;
        pushed_back = 0;
        return token_type;
    }
    return yylex();
}


void put_back(int token_type)
{
    pushed_back = token_type;
}


int match(int required_type)
{
    int token_type = next_token();
    if (token_type != required_type) {
        return 0;
    }
    return 1;
}


void str_tolower(char *s)
{
    for (; *s; s++) {
        *s = tolower(*s);
    }
}


int parse_state_spec(char *buf, char *rep_name, int *state_no)
{
    // parse first state descriptor
    char *p = strchr(buf, ']');
    size_t len = p - buf - 1;
    strcpy(rep_name, buf + 1);
    rep_name[len] = '\0';
    *state_no = atoi(p + 2);
}


/**
 * Parses sector specification: "<int>h<int>p" -> (h,p)
 * example: 0h1p -> (0,1)
 *
 * @note this function assumes 's' was classified by the lexer
 * as the TT_SECTOR-type token
 */
void parse_sector_label(char *s, int *h, int *p)
{
    *h = s[0] - '0';
    *p = s[2] - '0';
}


double match_float_number()
{
    int token_type = next_token();
    if (token_type != TT_FLOAT && token_type != TT_INTEGER) {
        yyerror("float number is expected");
    }

    return atof(yytext);
}


double match_positive_float_number()
{
    int token_type = next_token();
    if (token_type != TT_FLOAT && token_type != TT_INTEGER) {
        yyerror("positive float number is expected");
    }

    double num = atof(yytext);

    if (num <= 0.0) {
        yyerror("positive float number is expected");
    }

    return num;
}


int match_positive_integer()
{
    int token_type = next_token();
    if (token_type != TT_INTEGER) {
        yyerror("integer positive number is expected");
    }
    int number = atoi(yytext);
    if (number <= 0) {
        yyerror("positive number is expected");
    }

    return number;
}


int match_non_negative_integer()
{
    int token_type = next_token();
    if (token_type != TT_INTEGER) {
        yyerror("integer positive number is expected");
    }

    int number = atoi(yytext);

    return number;
}


void match_end_of_line()
{
    int token_type = next_token();
    if (token_type != END_OF_LINE) {
        yyerror("end of line is expected");
    }
}

