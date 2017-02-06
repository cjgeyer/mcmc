
#include <R.h>
#include <Rinternals.h>
#include "myutil.h"

int
getScalarLogical(SEXP foo, char *argname)
{
    if (! isLogical(foo))
        error("argument \"%s\" must be logical", argname);
    if (LENGTH(foo)  != 1)
        error("argument \"%s\" must be scalar", argname);
    return LOGICAL(foo)[0];
}

