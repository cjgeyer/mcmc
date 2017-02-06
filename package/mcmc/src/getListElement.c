
#include <R.h>
#include <Rinternals.h>
#include "myutil.h"

SEXP getListElement(SEXP list, char *str)
{
    SEXP elmt = R_NilValue;
    SEXP names = getAttrib(list, R_NamesSymbol);
    int i;

    if (names == R_NilValue)
        return R_NilValue;

    for (i = 0; i < length(list); i++)
        if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
            elmt = VECTOR_ELT(list, i);
            break;
        }
    return elmt;
}

