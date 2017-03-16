
#include <R.h>
#include <Rinternals.h>
#include "myutil.h"
#include "mcmc.h"

SEXP initseq(SEXP x)
{
    SEXP xreal;

    if (! isNumeric(x))
        error("argument must be numeric");

    PROTECT(xreal = coerceVector(x, REALSXP));
    if (! isAllFinite(x))
        error("all elements of argument must be finite");
    int len = LENGTH(xreal);

    double *buff = (double *) R_alloc(len / 2, sizeof(double));

    int i;
    double gamma_zero = 0.0;   /* for gcc -Wall -Wextra */

    for (i = 0; i < len / 2; ++i) {

        int lag1 = 2 * i;
        double gam1 = 0.0;
        for (int j = 0; j + lag1 < len; ++j)
            gam1 += REAL(xreal)[j] * REAL(xreal)[j + lag1];
        gam1 /= len;

        if (i == 0)
            gamma_zero = gam1;

        int lag2 = lag1 + 1;
        double gam2 = 0.0;
        for (int j = 0; j + lag2 < len; ++j)
            gam2 += REAL(xreal)[j] * REAL(xreal)[j + lag2];
        gam2 /= len;

        buff[i] = gam1 + gam2;
        if (buff[i] < 0.0) {
            buff[i] = 0.0;
            ++i;
            break;
        }
    }

    SEXP gamma_pos, gamma_dec, gamma_con;

    PROTECT(gamma_pos = allocVector(REALSXP, i));
    for (int j = 0; j < i; ++j)
        REAL(gamma_pos)[j] = buff[j];

    for (int j = 1; j < i; ++j)
        if (buff[j] > buff[j - 1])
            buff[j] = buff[j - 1];

    PROTECT(gamma_dec = allocVector(REALSXP, i));
    for (int j = 0; j < i; ++j)
        REAL(gamma_dec)[j] = buff[j];

    for (int j = i - 1; j > 0; --j)
        buff[j] -= buff[j - 1];

    /* Pool Adjacent Violators Algorithm (PAVA) */
    double *puff = (double *) R_alloc(i, sizeof(double));
    int *nuff = (int *) R_alloc(i, sizeof(int));
    int nstep = 0;
    for (int j = 1; j < i; ++j) {
        puff[nstep] = buff[j];
        nuff[nstep] = 1;
        ++nstep;
        while(nstep > 1 && puff[nstep - 1] / nuff[nstep - 1]
            < puff[nstep - 2] / nuff[nstep - 2]) {
            puff[nstep - 2] += puff[nstep - 1];
            nuff[nstep - 2] += nuff[nstep - 1];
            --nstep;
        }
    }

    for (int jstep = 0, j = 1; jstep < nstep; ++jstep) {
        double muff = puff[jstep] / nuff[jstep];
        for (int k = 0; k < nuff[jstep]; ++j, ++k)
            buff[j] = buff[j - 1] + muff;
    }

    PROTECT(gamma_con = allocVector(REALSXP, i));
    for (int j = 0; j < i; ++j)
        REAL(gamma_con)[j] = buff[j];

    double var_pos = 0.0;
    double var_dec = 0.0;
    double var_con = 0.0;
    for (int j = 0; j < i; ++j) {
        var_pos += REAL(gamma_pos)[j];
        var_dec += REAL(gamma_dec)[j];
        var_con += REAL(gamma_con)[j];
    }
    var_pos *= 2.0;
    var_dec *= 2.0;
    var_con *= 2.0;
    var_pos -= gamma_zero;
    var_dec -= gamma_zero;
    var_con -= gamma_zero;

    SEXP result, resultnames;
    PROTECT(result = allocVector(VECSXP, 7));
    PROTECT(resultnames = allocVector(STRSXP, 7));
    SET_VECTOR_ELT(result, 0, ScalarReal(gamma_zero));
    SET_STRING_ELT(resultnames, 0, mkChar("gamma0"));
    SET_VECTOR_ELT(result, 1, gamma_pos);
    SET_STRING_ELT(resultnames, 1, mkChar("Gamma.pos"));
    SET_VECTOR_ELT(result, 2, gamma_dec);
    SET_STRING_ELT(resultnames, 2, mkChar("Gamma.dec"));
    SET_VECTOR_ELT(result, 3, gamma_con);
    SET_STRING_ELT(resultnames, 3, mkChar("Gamma.con"));
    SET_VECTOR_ELT(result, 4, ScalarReal(var_pos));
    SET_STRING_ELT(resultnames, 4, mkChar("var.pos"));
    SET_VECTOR_ELT(result, 5, ScalarReal(var_dec));
    SET_STRING_ELT(resultnames, 5, mkChar("var.dec"));
    SET_VECTOR_ELT(result, 6, ScalarReal(var_con));
    SET_STRING_ELT(resultnames, 6, mkChar("var.con"));
    namesgets(result, resultnames);

    UNPROTECT(6);
    return result;
}

