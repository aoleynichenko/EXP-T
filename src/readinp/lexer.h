/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2020 The EXP-T developers.
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

#ifndef LEXER_H_INCLUDED
#define LEXER_H_INCLUDED

#include <stdio.h>

extern FILE *yyin;
extern char *yytext;
extern int yylineno;
extern int yycol;

int yylex();

enum {
    /* token types */
    TT_QUOTE = 0,
    TT_SECTOR,
    TT_STAR,
    TT_HYPHEN,
    TT_ELEC_STATE,
    TT_INTEGER,
    TT_FLOAT,
    TT_WORD,
    /* keywords */
    KEYWORD_TITLE,
    KEYWORD_PRINT,
    KEYWORD_SCRATCH_DIR,
    KEYWORD_SECTOR,
    KEYWORD_MODEL,
    KEYWORD_OCC,
    KEYWORD_OCC_IRREPS,
    KEYWORD_ACTIVE,
    KEYWORD_NACTP,
    KEYWORD_NACTH,
    KEYWORD_DEGEN_THRESH,
    KEYWORD_NOHERMIT,
    KEYWORD_DLTDM,
    KEYWORD_NATORB,
    KEYWORD_NROOTS,
    KEYWORD_MAXITER,
    KEYWORD_CONV,
    KEYWORD_DAMPING,
    KEYWORD_SHIFTTYPE,
    KEYWORD_ORBSHIFT,
    KEYWORD_ORBSHIFT00,
    KEYWORD_SHIFT,
    KEYWORD_REUSE,
    KEYWORD_FLUSH,
    KEYWORD_INTEGRALS,
    KEYWORD_ONEPROP,
    KEYWORD_TWOPROP,
    KEYWORD_MEMORY,
    KEYWORD_MEMORY_MB,
    KEYWORD_MEMORY_GB,
    KEYWORD_TILESIZE,
    KEYWORD_DISK_USAGE,
    KEYWORD_COMPRESS,
    KEYWORD_NTHREADS,
    KEYWORD_OPENMP,
    KEYWORD_CUDA,
    KEYWORD_ARITH,
    KEYWORD_PROP,
    KEYWORD_SINGLES,
    KEYWORD_CORE,
    KEYWORD_VIRTUAL,
    KEYWORD_END,
    KEYWORD_HPTT,
    KEYWORD_DIIS,
    KEYWORD_GAUNT,
    KEYWORD_BREIT,
    KEYWORD_X2CMMF,
    KEYWORD_NOINNER,
    /* special */
    END_OF_LINE,
    END_OF_FILE
};

#endif /* LEXER_H_INCLUDED */
