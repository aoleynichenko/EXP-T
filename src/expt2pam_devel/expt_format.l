%{

#include "expt_lexer.h"

#define MAX_TOKEN_LEN 1024

void inc_col();
int yycol = 0;
int prev_len = -1;

%}

%option case-insensitive

%%
"#"[^\n]*       { inc_col(); } /* single-line comments beginning with '#' */
[ \t\r]         { inc_col(); } /* ignore whitespace */
geometry        { inc_col(); return KEYWORD_GEOMETRY;     }
basis           { inc_col(); return KEYWORD_BASIS;        }
ecp             { inc_col(); return KEYWORD_ECP;          }
so              { inc_col(); return KEYWORD_SO;           }
end             { inc_col(); return KEYWORD_END;          }
\".*\"          { inc_col(); return TT_QUOTE;             }
[-]             { inc_col(); return TT_HYPHEN;            }
[*]             { inc_col(); return TT_STAR;              }
[0-9]*          { inc_col(); return TT_INTEGER;           }
E[0-9]*         { inc_col(); return TT_WORD;              }  /* for E120 symbols etc */
[+-]?[0-9]*(\.[0-9]*)?([eE][+-]?[0-9]+)?   { inc_col(); return TT_FLOAT; }
([^ *\t\n\[\]])+ { inc_col(); return TT_WORD;     }
\n              { inc_col(); yylineno++; return END_OF_LINE; }
<<EOF>>         { return END_OF_FILE; }
%%

/* just dummy */
int yywrap() {
   return 1;
}

/**
 * updates the number of the current column
 */
void inc_col()
{
    // update yycol
    if (prev_len == -1) {
        yycol = 0;
    }
    else {
        yycol += prev_len;
    }

    // update prev_len
    if (strcmp(yytext, "\n") == 0) {
        prev_len = -1;
    }
    else {
        prev_len = strlen(yytext);
    }
}
