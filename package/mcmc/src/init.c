
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "mcmc.h"

static R_NativePrimitiveArgType olbm_types[7] =
    {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, LGLSXP};

static R_CMethodDef cMethods[] = {
    {"olbm", (DL_FUNC) &olbm, 7, olbm_types},
    {NULL, NULL, 0, NULL}
};

static R_CallMethodDef callMethods[]  = {
    {"metrop", (DL_FUNC) &metrop, 10},
    {"temper", (DL_FUNC) &temper, 12},
    {"initseq", (DL_FUNC) &initseq, 1},
    {NULL, NULL, 0}
};

void attribute_visible R_init_mcmc(DllInfo *info)
{
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}

