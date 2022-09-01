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

#ifndef EXPT_LEXER_H_INCLUDED
#define EXPT_LEXER_H_INCLUDED

#include <stdio.h>

extern FILE *yyin;
extern char *yytext;
extern int yylineno;
extern int yycol;

int yylex();

enum {
    /* token types */
    TT_QUOTE = 0,
    TT_STAR,
    TT_HYPHEN,
    TT_INTEGER,
    TT_FLOAT,
    TT_WORD,
    /* keywords */
    KEYWORD_GEOMETRY,
    KEYWORD_BASIS,
    KEYWORD_ECP,
    KEYWORD_SO,
    KEYWORD_END,
    /* special */
    END_OF_LINE,
    END_OF_FILE
};

#endif /* EXPT_LEXER_H_INCLUDED */
