
#include <R.h>
#include <Rinternals.h>
#include "myutil.h"

int
isAllFinite(SEXP foo)
{
    int d, i;
    int result = TRUE;

    if (! isReal(foo))
        error("argument must be real");

    d = LENGTH(foo);
    for (i = 0; i < d; i++)
        result &= R_finite(REAL(foo)[i]);
    return result;
}

