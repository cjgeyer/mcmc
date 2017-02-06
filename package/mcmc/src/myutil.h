
#ifndef MCMC_MYUTIL_H
#define MCMC_MYUTIL_H

#include <R.h>
#include <Rinternals.h>

SEXP getListElement(SEXP list, char *str);
int getScalarInteger(SEXP foo, char *argname);
int getScalarLogical(SEXP foo, char *argname);
int isAllFinite(SEXP foo);

#endif /* MCMC_MYUTIL_H */

