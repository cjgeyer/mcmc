
#include <R.h>
#include <Rinternals.h>
#include "myutil.h"

int
getScalarInteger(SEXP foo, char *argname)
{
    if (! isNumeric(foo))
        error("argument \"%s\" must be numeric", argname);
    if (LENGTH(foo)  != 1)
        error("argument \"%s\" must be scalar", argname);
    if (isInteger(foo)) {
        return INTEGER(foo)[0];
    } else {
        SEXP bar = coerceVector(foo, INTSXP);
        return INTEGER(bar)[0];
    }
}

