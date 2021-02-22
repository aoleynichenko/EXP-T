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
 * readinp.c
 *
 * Reads CC input file.
 * Lexical analysis is performed with the automatically generated code
 * (by lex, see expt.l)
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lexer.h"
#include "platform.h"
#include "error.h"
#include "options.h"

#define MAX_LINE_LEN 1024
#define MAXIMUM2(x, y) (((x) > (y)) ? (x) : (y))

// functions defined in this file:
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
void directive_maxiter(cc_options_t *opts);
void directive_conv(cc_options_t *opts);
void directive_div_thresh(cc_options_t *opts);
void directive_damping(cc_options_t *opts);
void directive_diis(cc_options_t *opts);
void directive_shifttype(cc_options_t *opts);
void directive_orbshift(cc_options_t *opts);
void directive_orbshift00(cc_options_t *opts);
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
void directive_arith(cc_options_t *opts);
void directive_mdprop(cc_options_t *opts);
void directive_txtprop(cc_options_t *opts);
void directive_interface(cc_options_t *opts);

void yyerror(char *s);
int next_token();
void put_back(int token_type);
int match(int required_type);
void str_tolower(char *s);
int parse_state_spec(char *buf, char *rep_name, int *state_no);
void read_space_specification(cc_space_t *space);
void set_input_file_name(char *file_name);
void parse_sector(char *s, int *h, int *p);


/*******************************************************************************
 * readinp
 *
 * Read options from input file, store them in object 'opts'
 ******************************************************************************/
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
            case KEYWORD_NOHERMIT:
                opts->do_hermit = 0;
                break;
            case KEYWORD_MSTDM:
                opts->do_diplen_tdm = 1;
                break;
            case KEYWORD_NATORB:
                directive_natorb(opts);
                break;
            case KEYWORD_NROOTS:
                directive_nroots(opts);
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
            case KEYWORD_SHIFTTYPE:
                directive_shifttype(opts);
                break;
            case KEYWORD_ORBSHIFT:
                directive_orbshift(opts);
                break;
            case KEYWORD_ORBSHIFT00:
                directive_orbshift00(opts);
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
            case KEYWORD_CUDA:
                opts->cuda_enabled = 1;
                break;
            case KEYWORD_OPENMP:
            case KEYWORD_NTHREADS:
                directive_nthreads(opts);
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
            case KEYWORD_INTERFACE:
                directive_interface(opts);
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

    // triples are available in the development version only
#ifndef VERSION_DEVEL
    if (opts->cc_model > CC_MODEL_CCSD) {
        errquit("CC models including triples are not part of the public release");
    }
#endif
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
                       "Possible values: low, medium, high, debug";

    if (!match(TT_WORD)) {
        yyerror(msg);
    }
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

    parse_sector(yytext, &opts->sector_h, &opts->sector_p);
}


/**
 * Syntax:
 * one of: ( ccs || ccd || ccsd || ccsd(t) || ccsd+t(3) || ccsdt-1 || ccsdt-1' || ccsdt-2 || ccsdt-3 )
 */
void parse_cc_model(cc_options_t *opts)
{
    static char *msg = "wrong coupled cluster model!\n"
                       "Possible values: ccs, ccd, ccsd, ccsd+t(3)=ccsd(t), ccsdt-1a, ccsdt-1b=ccsdt-1, ccsdt-2, ccsdt-3";

    /*if (!match(TT_WORD)) {
        yyerror(msg);
    }*/
    str_tolower(yytext);
    if (strcmp(yytext, "ccs") == 0) {
        opts->cc_model = CC_MODEL_CCS;
        strcpy(opts->cc_model_str, "CCS");
    }
    else if (strcmp(yytext, "ccd") == 0) {
        opts->cc_model = CC_MODEL_CCD;
        strcpy(opts->cc_model_str, "CCD");
    }
    else if (strcmp(yytext, "ccsd") == 0) {
        opts->cc_model = CC_MODEL_CCSD;
        strcpy(opts->cc_model_str, "CCSD");
    }
    else if (strcmp(yytext, "ccsd+t(3)") == 0) {
        opts->cc_model = CC_MODEL_CCSD_T3;
        strcpy(opts->cc_model_str, "CCSD(T)/CCSD+T(3)");
    }
    else if (strcmp(yytext, "ccsd+t*(3)") == 0) {
        opts->cc_model = CC_MODEL_CCSD_T3_STAR;
        strcpy(opts->cc_model_str, "CCSD(T)/CCSD+T*(3)");
    }
    else if (strcmp(yytext, "ccsd+t(4)") == 0) {
        opts->cc_model = CC_MODEL_CCSD_T4;
        strcpy(opts->cc_model_str, "CCSD(T)/CCSD+T(4)");
    }
    else if (strcmp(yytext, "ccsd+t*(4)") == 0) {
        opts->cc_model = CC_MODEL_CCSD_T4_STAR;
        strcpy(opts->cc_model_str, "CCSD(T)/CCSD+T*(4)");
    }
    else if (strcmp(yytext, "ccsd(t)") == 0) {
        opts->cc_model = CC_MODEL_CCSD_T3;
        strcpy(opts->cc_model_str, "CCSD(T)/CCSD+T(3)");
    }
    else if (strcmp(yytext, "ccsdt-1a") == 0) {
        opts->cc_model = CC_MODEL_CCSDT_1A;
        strcpy(opts->cc_model_str, "CCSDT-1a");
    }
    else if (strcmp(yytext, "ccsdt-1b") == 0) {
        opts->cc_model = CC_MODEL_CCSDT_1B;
        strcpy(opts->cc_model_str, "CCSDT-1b");
    }
    else if (strcmp(yytext, "ccsdt-1") == 0) {
        opts->cc_model = CC_MODEL_CCSDT_1B;
        strcpy(opts->cc_model_str, "CCSDT-1b");
    }
    else if (strcmp(yytext, "ccsdt-1'") == 0) {
        opts->cc_model = CC_MODEL_CCSDT_1B_PRIME;
        strcpy(opts->cc_model_str, "CCSDT-1b'");
    }
    else if (strcmp(yytext, "ccsdt-2") == 0) {
        opts->cc_model = CC_MODEL_CCSDT_2;
        strcpy(opts->cc_model_str, "CCSDT-2");
    }
    else if (strcmp(yytext, "ccsdt-3") == 0) {
        opts->cc_model = CC_MODEL_CCSDT_3;
        strcpy(opts->cc_model_str, "CCSDT-3");
    }
    else if (strcmp(yytext, "ccsdt") == 0) {
        opts->cc_model = CC_MODEL_CCSDT;
        strcpy(opts->cc_model_str, "CCSDT");
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
    static char *msg = "wrong specification of irrep occupation numbers!\n"
                       "A list of non-negative integers separated by spaces is expected";
    int token_type;
    int count = 0;

    opts->nelec_defined = 1;

    token_type = next_token();
    while (token_type != END_OF_LINE && token_type != END_OF_FILE) {
        if (token_type != TT_INTEGER) {
            yyerror(msg);
        }

        // note: TT_INTEGER means that the value is already non-negative
        opts->nelec[count++] = atoi(yytext);

        token_type = next_token();
    }
    put_back(token_type);
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
            opts->actsp_min = atof(yytext);

            // upper limit
            token_type = next_token();
            if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                yyerror(msg1);
            }
            opts->actsp_max = atof(yytext);

            opts->actsp_defined = 1;
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
    opts->actsp_defined = 4;
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
        opts->actsp_defined = 2;
        return;
    }
    else if (token_type == TT_ELEC_STATE) {
        do {
            char rep_name[32];
            int nact;
            int written = 0;

            parse_state_spec(yytext, rep_name, &nact);
            // try to find entry
            for (int i = 0; i < CC_MAX_NREP; i++) {
                cc_active_spec_t *sp = &opts->active_specs[i];
                if (strcmp(sp->rep_name, rep_name) == 0) {
                    sp->nactp = nact;
                    written = 1;
                    break;
                }
            }
            // if entry not found -- write to the first empty slot
            if (!written) {
                for (int i = 0; i < CC_MAX_NREP; i++) {
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
        opts->actsp_defined = 3;
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
        opts->actsp_defined = 2;
        return;
    }
    else if (token_type == TT_ELEC_STATE) {
        do {
            char rep_name[32];
            int nact;
            int written = 0;

            parse_state_spec(yytext, rep_name, &nact);
            // try to find entry
            for (int i = 0; i < CC_MAX_NREP; i++) {
                cc_active_spec_t *sp = &opts->active_specs[i];
                if (strcmp(sp->rep_name, rep_name) == 0) {
                    sp->nacth = nact;
                    written = 1;
                }
            }
            // if entry not found -- write to the first empty slot
            if (!written) {
                for (int i = 0; i < CC_MAX_NREP; i++) {
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
        opts->actsp_defined = 3;
    }
    else {
        yyerror(msg);
    }
}


/**
 * Syntax:
 * degen_thresh <real thresh>
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
    opts->degen_thresh = thresh;
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
    opts->conv = thresh;
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
 * diis [ off || <integer> ]
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
        opts->diis_dim = dim;
        return;
    }
    else {
        yyerror(msg);
    }
}


/**
 * Syntax:
 * shifttype ( none || real || realimag || imag || taylor )
 */
void directive_shifttype(cc_options_t *opts)
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
 * orbshift [<integer power default 3>] <real shift s>
 */
void directive_orbshift(cc_options_t *opts)
{
    static char *msg = "wrong specification of orbital shifts!";
    static char buf[256];
    int token_type;

    // attenuation parameter or shift ?
    token_type = next_token();
    if (token_type == TT_INTEGER || token_type == TT_FLOAT) {
        strcpy(buf, yytext);
    }
    else {
        yyerror(msg);
    }

    token_type = next_token();
    if (token_type == TT_FLOAT) {
        // orbshift <n> <sh>
        opts->do_orbshift = 1;
        opts->orbshift = atof(yytext);
        opts->orbshift_power = atoi(buf);
    }
    else {
        // orbshift <sh>
        opts->do_orbshift = 1;
        opts->orbshift = atof(buf);
        opts->orbshift_power = 3;
        put_back(token_type);
    }
}


/**
 * Syntax:
 * orbshift00 [<integer power default 3>] <real shift s>
 */
void directive_orbshift00(cc_options_t *opts)
{
    static char *msg = "wrong specification of orbital shifts in the 0h0p sector!";
    static char buf[256];
    int token_type;

    token_type = next_token();
    if (token_type == TT_INTEGER) {
        opts->orbshift_0h0p_power = atoi(yytext);
    }
    else {
        yyerror(msg);
    }

    int i = 0;
    token_type = next_token();
    while (token_type == TT_INTEGER || token_type == TT_FLOAT) {
        opts->orbshift_0h0p[i++] = atof(yytext);
        token_type = next_token();
    }
    put_back(token_type);

    if (i == 1) { // one parameter for all amplitudes independent on the number of valence lines
        opts->do_orbshift_0h0p = 1;
    }
    else { // individual shift parameter for each number of valence lines
        opts->do_orbshift_0h0p = 2;
    }

    // attenuation parameter or shift ?
    /*token_type = next_token();
    if (token_type == TT_INTEGER || token_type == TT_FLOAT) {
        strcpy(buf, yytext);
    }
    else {
        yyerror(msg);
    }

    token_type = next_token();
    if (token_type == TT_FLOAT) {
        // orbshift00 <n> <sh>
        opts->do_orbshift_0h0p = 1;
        opts->orbshift_0h0p = atof(yytext);
        opts->orbshift_0h0p_power = atoi(buf);
    }
    else {
        // orbshift00 <sh>
        opts->do_orbshift_0h0p = 1;
        opts->orbshift_0h0p = atof(buf);
        opts->orbshift_0h0p_power = 3;
        put_back(token_type);
    }*/
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

    // disable orbital shifts
    opts->do_orbshift = 0;
    opts->do_orbshift_0h0p = 0;
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
                       "Possible values: integrals 1-integrals 2-integrals amplitudes 0h0p 0h1p 1h0p 0h2p 2h0p 1h1p";
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
            opts->reuse_0h0p = 1;
            opts->reuse_0h1p = 1;
            opts->reuse_1h1p = 1;
            opts->reuse_0h2p = 1;
            opts->reuse_1h0p = 1;
            opts->reuse_2h0p = 1;
            opts->reuse_0h3p = 1;
        }
        else if (strcmp(yytext, "0h0p") == 0) {
            opts->reuse_0h0p = 1;
        }
        else if (strcmp(yytext, "0h1p") == 0) {
            opts->reuse_0h1p = 1;
        }
        else if (strcmp(yytext, "1h0p") == 0) {
            opts->reuse_1h0p = 1;
        }
        else if (strcmp(yytext, "1h1p") == 0) {
            opts->reuse_1h1p = 1;
        }
        else if (strcmp(yytext, "0h2p") == 0) {
            opts->reuse_0h2p = 1;
        }
        else if (strcmp(yytext, "2h0p") == 0) {
            opts->reuse_2h0p = 1;
        }
        else if (strcmp(yytext, "0h3p") == 0) {
            opts->reuse_0h3p = 1;
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

        parse_sector(yytext, &h, &p);
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
    cc_ms_prop_query_t *q = &opts->prop_queries[opts->n_ms_props];
    strcpy(q->prop_name, yytext + 1);
    q->do_transpose = 0;
    q->source = CC_PROP_FROM_MDPROP;

    // optional "transpose" keyword
    token_type = next_token();
    if (token_type == TT_WORD) {
        if (strcmp(yytext, "transpose") == 0) {
            q->do_transpose = 1;
        }
        else {
            yyerror(msg_trans);
        }
    }
    else {
        put_back(token_type);
    }

    opts->n_ms_props++;
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

    cc_ms_prop_query_t *q = &opts->prop_queries[opts->n_ms_props];
    q->source = CC_PROP_FROM_TXTPROP;
    q->do_transpose = 0;

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

    // optional "transpose" keyword
    token_type = next_token();
    if (token_type == TT_WORD) {
        if (strcmp(yytext, "transpose") == 0) {
            q->do_transpose = 1;
        }
        else {
            yyerror(msg_trans);
        }
    }
    else {
        put_back(token_type);
    }

    opts->n_ms_props++;
}


/**
 * Syntax:
 * interface ( dirac || tm2c )
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
        else if (strcmp(yytext, "tm2c") == 0) {
            cc_opts->int_source = CC_INTEGRALS_TM2C;
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
    static char *msg = "wrong specification of the space of spinors!";
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
            for (int i = 0; i < CC_MAX_NREP; i++) {
                if (strcmp(space->rep_names[i], rep_name) == 0) {
                    space->dim[i] = norb;
                    written = 1;
                    break;
                }
            }
            // if entry not found -- write to the first empty slot
            if (!written) {
                for (int i = 0; i < CC_MAX_NREP; i++) {
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
 *   0h0p triples off     # set T3 = 0
 *   0h1p doubles -5 20   # T2 != 0 only for excitations with energy in the given range
 *   0h2p triples spectator   # only spectator triples in 0h2p are nonzero
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
        parse_sector(yytext, &sect_h, &sect_p);
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
 * no_inner_core_corr <real thresh>
 */
void directive_no_inner_core_corr(cc_options_t *opts)
{
    static char *msg = "wrong specification of no_inner_core_corr!\n"
                       "A real number is expected";
    double thresh;
    int token_type;

    token_type = next_token();
    if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
        yyerror(msg);
    }

    thresh = atof(yytext);

    opts->do_remove_inner_core_corr = 1;
    opts->inner_core_thresh = thresh;
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
void parse_sector(char *s, int *h, int *p)
{
    *h = s[0] - '0';
    *p = s[2] - '0';
}
