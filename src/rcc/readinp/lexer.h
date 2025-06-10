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
    TT_EQ,
    TT_NEQ,
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
    KEYWORD_HERMIT,
    KEYWORD_MSTDM,
    KEYWORD_NATORB,
    KEYWORD_NROOTS,
    KEYWORD_ROOTS_CUTOFF,
    KEYWORD_MAXITER,
    KEYWORD_CONV,
    KEYWORD_DIV_THRESH,
    KEYWORD_DAMPING,
    KEYWORD_SHIFT_TYPE,
    KEYWORD_SHIFT,
    KEYWORD_REUSE,
    KEYWORD_FLUSH,
    KEYWORD_INTEGRALS,
    KEYWORD_ONEPROP,
    KEYWORD_TWOPROP,
    KEYWORD_MEMORY,
    KEYWORD_TILESIZE,
    KEYWORD_DISK_USAGE,
    KEYWORD_COMPRESS,
    KEYWORD_COMPRESS_TRIPLES,
    KEYWORD_NTHREADS,
    KEYWORD_OPENMP,
    KEYWORD_OPENMP_ALGORITHM,
    KEYWORD_CUDA,
    KEYWORD_ARITH,
    KEYWORD_MDPROP,
    KEYWORD_TXTPROP,
    KEYWORD_END,
    KEYWORD_SELECT,
    KEYWORD_IH_IMMS,
    KEYWORD_DIIS,
    KEYWORD_CROP,
    KEYWORD_GAUNT,
    KEYWORD_BREIT,
    KEYWORD_X2CMMF,
    KEYWORD_NEW_SORTING,
    KEYWORD_SKIP,
    KEYWORD_INTERFACE,
    KEYWORD_RESTRICT_T3,
    KEYWORD_FLUSH_AMPLITUDES_TXT,
    KEYWORD_ANALYT_PROP,
    KEYWORD_DENSITY,
    KEYWORD_LAMBDA,
    KEYWORD_OVERLAP,
    KEYWORD_SPINOR_LABELS,
    KEYWORD_HUGHES_KALDOR_1H2P,
    KEYWORD_HUGHES_KALDOR_2H1P,
    KEYWORD_USE_ORB_ENERGIES,
    KEYWORD_RECALC_ORB_ENERGIES,
    /* tensor trains */
    KEYWORD_USE_TT_CCSD,
    KEYWORD_TT_SVD_TOL,
    KEYWORD_TT_CHOLESKY_TOL,
    KEYWORD_TT_MULT_PPPP,
    KEYWORD_TT_DIIS,
    KEYWORD_TENSOR_TRAINS,
    KEYWORD_GOLDSTONE,
    /* special */
    END_OF_LINE,
    END_OF_FILE
};

#endif /* LEXER_H_INCLUDED */
