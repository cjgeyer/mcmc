
/*
*
* mcmc and MCMC package for R
* Copyright (c) 2005 Charles J. Geyer
*
* All rights reserved.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, and/or sell copies of the
* Software, and to permit persons to whom the Software is furnished to do so,
* provided that the above copyright notice(s) and this permission notice appear
* in all copies of the Software and that both the above copyright notice(s) and
* this permission notice appear in supporting documentation.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY RIGHTS.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS NOTICE
* BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL DAMAGES,
* OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
* WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
* ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*
* Except as contained in this notice, the name of a copyright holder shall
* not be used in advertising or otherwise to promote the sale, use or other
* dealings in this Software without prior written authorization of the
* copyright holder.
*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "myutil.h"

static void proposal_setup(SEXP scale, int d);

static void propose(SEXP state, SEXP proposal);

static double logh(SEXP func, SEXP state, SEXP rho);

static int out_setup(SEXP func, SEXP rho, SEXP state);

static void outfun(SEXP state, SEXP buffer);

SEXP metrop(SEXP func1, SEXP initial, SEXP nbatch, SEXP blen, SEXP nspac,
    SEXP scale, SEXP func2, SEXP debug, SEXP rho1, SEXP rho2)
{
    int int_nbatch, int_blen, int_nspac, int_debug;
    SEXP state, proposal;
    int dim_state, dim_out;
    SEXP result, resultnames, acceptance_rate, path,
        save_initial, save_final;
    double *batch_buffer;
    SEXP out_buffer;
    int ibatch, jbatch, ispac;

    int i, k;
    double acceptances = 0.0;
    double tries = 0.0;

    double current_log_dens;

    if (! isFunction(func1))
        error("argument \"func1\" must be function");
    if (! isEnvironment(rho1))
        error("argument \"rho1\" must be environment");

    if (! isNumeric(initial))
        error("argument \"initial\" must be numeric");
    if (! isNumeric(nbatch))
        error("argument \"nbatch\" must be numeric");
    if (! isNumeric(blen))
        error("argument \"blen\" must be numeric");
    if (! isNumeric(nspac))
        error("argument \"nspac\" must be numeric");
    if (! isNumeric(scale))
        error("argument \"scale\" must be numeric");

    if (! isLogical(debug))
        error("argument \"debug\" must be logical");

    int_nbatch = getScalarInteger(nbatch);
    int_blen = getScalarInteger(blen);
    int_nspac = getScalarInteger(nspac);

    int_debug = getScalarLogical(debug);

    if (int_nbatch <= 0)
        error("argument \"nbatch\" must be positive");
    if (int_blen <= 0)
        error("argument \"blen\" must be positive");
    if (int_nspac <= 0)
        error("argument \"nspac\" must be positive");

    PROTECT(state = coerceVector(duplicate(initial), REALSXP));
    if (! isAllFinite(state))
        error("all elements of \"state\" must be finite");
    dim_state = LENGTH(state);

    PROTECT(proposal = allocVector(REALSXP, dim_state));
    proposal_setup(scale, dim_state);

    dim_out = out_setup(func2, rho2, state);
    batch_buffer = (double *) R_alloc(dim_out, sizeof(double));
    PROTECT(out_buffer = allocVector(REALSXP, dim_out));

     if (! int_debug) {
         PROTECT(result = allocVector(VECSXP, 4));
         PROTECT(resultnames = allocVector(STRSXP, 4));
     } else {
         PROTECT(result = allocVector(VECSXP, 8));
         PROTECT(resultnames = allocVector(STRSXP, 8));
     }
     PROTECT(acceptance_rate = allocVector(REALSXP, 1));
     SET_VECTOR_ELT(result, 0, acceptance_rate);
     PROTECT(path = allocMatrix(REALSXP, dim_out, int_nbatch));
     SET_VECTOR_ELT(result, 1, path);
     PROTECT(save_initial = duplicate(state));
     SET_VECTOR_ELT(result, 2, save_initial);
     UNPROTECT(3);
     SET_STRING_ELT(resultnames, 0, mkChar("accept"));
     SET_STRING_ELT(resultnames, 1, mkChar("batch"));
     SET_STRING_ELT(resultnames, 2, mkChar("initial"));
     SET_STRING_ELT(resultnames, 3, mkChar("final"));
     if (int_debug) {
         SEXP spath, ppath, gpath, upath;
         int nn = int_nbatch * int_blen * int_nspac;
         PROTECT(spath = allocMatrix(REALSXP, dim_state, nn));
         SET_VECTOR_ELT(result, 4, spath);
         PROTECT(ppath = allocMatrix(REALSXP, dim_state, nn));
         SET_VECTOR_ELT(result, 5, ppath);
         PROTECT(gpath = allocVector(REALSXP, nn));
         SET_VECTOR_ELT(result, 6, gpath);
         PROTECT(upath = allocVector(REALSXP, nn));
         SET_VECTOR_ELT(result, 7, upath);
         UNPROTECT(4);
         SET_STRING_ELT(resultnames, 4, mkChar("current"));
         SET_STRING_ELT(resultnames, 5, mkChar("proposal"));
         SET_STRING_ELT(resultnames, 6, mkChar("log.green"));
         SET_STRING_ELT(resultnames, 7, mkChar("u"));
     }
     namesgets(result, resultnames);
     UNPROTECT(1);

     GetRNGstate();

     current_log_dens = logh(func1, state, rho1);
     if (current_log_dens == R_NegInf)
         error("log unnormalized density -Inf at initial state");

     /* REVISED DOWN TO HERE */

     for (ibatch = 0, k = 0; ibatch < int_nbatch; ibatch++) {

         int j;

         for (i = 0; i < dim_out; i++)
             batch_buffer[i] = 0.0;

         for (jbatch = 0; jbatch < int_blen; jbatch++) {

             double proposal_log_dens;

             for (ispac = 0; ispac < int_nspac; ispac++) {

                 int accept;
                 double u = -1.0; /* impossible return from unif_rand() */

                 /* Note: should never happen! */
                 if (current_log_dens == R_NegInf)
                     error("log density -Inf at current state");

                 propose(state, proposal);

                 proposal_log_dens = logh(func1, proposal, rho1);

                 accept = FALSE;
                 if (proposal_log_dens != R_NegInf) {
                     if (proposal_log_dens > current_log_dens) {
                         accept = TRUE;
                     } else {
                         double green = exp(proposal_log_dens
                             - current_log_dens);
                         u = unif_rand();
                         accept = u < green;
                     }
                 }

                 if (int_debug) {
                     int l = ispac + int_nspac * (jbatch + int_blen * ibatch);
                     int lbase = l * dim_state;
                     SEXP spath = VECTOR_ELT(result, 4);
                     SEXP ppath = VECTOR_ELT(result, 5);
                     SEXP gpath = VECTOR_ELT(result, 6);
                     SEXP upath = VECTOR_ELT(result, 7);
                     int lj;
                     for (lj = 0; lj < dim_state; lj++) {
                         REAL(spath)[lbase + lj] = REAL(state)[lj];
                         REAL(ppath)[lbase + lj] = REAL(proposal)[lj];
                     }
                     REAL(gpath)[l] = proposal_log_dens - current_log_dens;
                     if (u == -1.0)
                         REAL(upath)[l] = NA_REAL;
                     else
                         REAL(upath)[l] = u;
                 }

                 if (accept) {
                     int jj;
                     for (jj = 0; jj < dim_state; jj++)
                         REAL(state)[jj] = REAL(proposal)[jj];
                     current_log_dens = proposal_log_dens;
                     acceptances++;
                 }
                 tries++;
             } /* end of inner loop (one iteration) */

             outfun(state, out_buffer);
             for (j = 0; j < dim_out; j++)
                 batch_buffer[j] += REAL(out_buffer)[j];

         } /* end of middle loop (one batch) */

         for (j = 0; j < dim_out; j++, k++)
             REAL(path)[k] = batch_buffer[j] / int_blen;

     } /* end of outer loop */

     PutRNGstate();

     REAL(acceptance_rate)[0] = acceptances / tries;

     PROTECT(save_final = coerceVector(state, REALSXP));
     SET_VECTOR_ELT(result, 3, save_final);

     UNPROTECT(5);
     return result;
}

static double logh(SEXP func, SEXP state, SEXP rho)
{
     SEXP call, result, foo;
     double bar;

     PROTECT(call = lang2(func, state));
     PROTECT(result = eval(call, rho));
     if (! isNumeric(result))
         error("logh: result of function call must be numeric");
     if (LENGTH(result) != 1)
         error("logh: result of function call must be scalar");
     PROTECT(foo = coerceVector(result, REALSXP));
     bar = REAL(foo)[0];
     UNPROTECT(3);
     if (bar == R_PosInf)
         error("logh: func returned +Inf");
     if (R_IsNaN(bar) || R_IsNA(bar))
         error("logh: func returned NA or NaN");
     /* Note: -Inf is allowed */
     return bar;
}

static double *scale_factor;
static double scale_factor_buffer;
static int scale_option;
static int state_dimension;
#define CONSTANT   1
#define DIAGONAL   2
#define FULL       3

static void proposal_setup(SEXP scale, int d)
{
    SEXP foo;

    state_dimension = d;

    PROTECT(foo = coerceVector(scale, REALSXP));
    if (isMatrix(scale)) {
        SEXP bar;
        PROTECT(bar = getAttrib(scale, R_DimSymbol));
        if (INTEGER(bar)[0] == d && INTEGER(bar)[1] == d) {
            int i;
            scale_factor = (double *) R_alloc(d * d, sizeof(double));
            for (i = 0; i < d * d; i++)
                scale_factor[i] = REAL(foo)[i];
            scale_option = FULL;
        } else {
            error("dimensions of \"scale\" matrix not d by d");
        }
        UNPROTECT(1);
    } else if (LENGTH(foo) == d) {
        int i;
        scale_factor = (double *) R_alloc(d, sizeof(double));
        for (i = 0; i < d; i++)
            scale_factor[i] = REAL(foo)[i];
        scale_option = DIAGONAL;
    } else if (LENGTH(foo) == 1) {
        scale_factor = &scale_factor_buffer;
        scale_factor[0] = REAL(foo)[0];
        scale_option = CONSTANT;
    } else {
        error("length of \"scale\" vector not d or 1");
    }
    UNPROTECT(1);
}

static void propose(SEXP state, SEXP proposal)
{
    int d = state_dimension;
    int i, j, k;

    if (scale_option == 0)
        error("attempt to call propose without setup");

    if (LENGTH(state) != d || LENGTH(proposal) != d)
        error("State or proposal length different from initialization\n");

    switch (scale_option) {
        case CONSTANT:
            for (j = 0; j < d; j++)
                REAL(proposal)[j] = REAL(state)[j]
                    + scale_factor[0] * norm_rand();
            break;
        case DIAGONAL:
            for (j = 0; j < d; j++)
                REAL(proposal)[j] = REAL(state)[j]
                    + scale_factor[j] * norm_rand();
            break;
        case FULL:
            for (j = 0; j < d; j++)
                REAL(proposal)[j] = REAL(state)[j];

            for (i = 0, k = 0; i < d; i++) {
                double u = norm_rand();
                for (j = 0; j < d; j++)
                    REAL(proposal)[j] += scale_factor[k++] * u;
            }
            break;
        default:
            error("bogus scaling option\n");
    }
}

static SEXP out_func;
static SEXP out_env;
static int *out_index;
static int out_option;
static int out_dimension;
static int out_state_dimension;
#define OUT_FUNCTION   1
#define OUT_INDEX      2
#define OUT_IDENTITY   3

static int out_setup(SEXP func, SEXP rho, SEXP state)
{
    out_state_dimension = LENGTH(state);

    if (func == R_NilValue) {
        out_option = OUT_IDENTITY;
        out_dimension = out_state_dimension;
        out_func = R_NilValue;
        out_env = R_NilValue;
    } else if (isFunction(func)) {
        if (! isEnvironment(rho))
            error("out_setup: argument \"rho\" must be environment");
        out_option = OUT_FUNCTION;
        out_func = func;
        out_env = rho;
        out_dimension = LENGTH(eval(lang2(func, state), rho));
    } else if (isLogical(func)) {
        int i;
        if (LENGTH(func) != out_state_dimension)
            error("is.logical(outfun) & (length(outfun) != length(initial))");
        out_option = OUT_INDEX;
        out_index = (int *) R_alloc(out_state_dimension, sizeof(int));
        for (i = 0, out_dimension = 0; i < out_state_dimension; i++) {
            out_index[i] = LOGICAL(func)[i];
            out_dimension += out_index[i];
        }
    } else if (isNumeric(func)) {
        SEXP foo;
        int foolen, i;
        int foopos = 0;
        int fooneg = 0;
        PROTECT(foo = coerceVector(func, REALSXP));
        foolen = LENGTH(foo);
        for (i = 0; i < foolen; i++) {
            double foodble = REAL(foo)[i];
            int fooint = foodble;
            int fooabs = fooint > 0 ? fooint : (- fooint);

            if (foodble == 0)
                error("is.numeric(outfun) & any(outfun == 0)");
            if (foodble != fooint)
                error("is.numeric(outfun) & any(outfun != as.integer(outfun))");
            if (fooabs > out_state_dimension)
                error("is.numeric(outfun) & any(abs(outfun) > length(initial)");

            if (foodble > 0)
                foopos++;
            else if (foodble < 0)
                fooneg++;
        }

        if ((foopos > 0) && (fooneg > 0))
            error("is.numeric(outfun) & any(outfun > 0) & any(outfun < 0)");

        out_option = OUT_INDEX;
        out_index = (int *) R_alloc(out_state_dimension, sizeof(int));
        if (foopos > 0) {
            for (i = 0; i < out_state_dimension; i++)
                out_index[i] = FALSE;
            for (i = 0; i < foolen; i++) {
                 int fooint = REAL(foo)[i];
                 out_index[fooint - 1] = TRUE;
            }
        } else /* (fooneg > 0) */ {
            for (i = 0; i < out_state_dimension; i++)
                out_index[i] = TRUE;
            for (i = 0; i < foolen; i++) {
                 int fooint = REAL(foo)[i];
                 int fooabs = (- fooint);
                 out_index[fooabs - 1] = FALSE;
            }
        }
        for (i = 0, out_dimension = 0; i < out_state_dimension; i++)
            out_dimension += out_index[i];
        UNPROTECT(1);
    }
    return out_dimension;
}

static void outfun(SEXP state, SEXP buffer)
{
    int j, k;

    if (out_option == 0)
        error("attempt to call outfun without setup");

    if (LENGTH(state) != out_state_dimension)
        error("outfun: state length different from initialization");
    if (! isReal(buffer))
        error("outfun: buffer must be real");
    if (LENGTH(buffer) != out_dimension)
        error("outfun: buffer length different from initialization");

    switch (out_option) {
        case OUT_IDENTITY:
            for (j = 0; j < out_state_dimension; j++)
                REAL(buffer)[j] = REAL(state)[j];
            break;
        case OUT_INDEX:
            for (j = 0, k = 0; j < out_state_dimension; j++)
                if (out_index[j])
                    REAL(buffer)[k++] = REAL(state)[j];
            break;
        case OUT_FUNCTION:
            {
                SEXP call, result, foo;

                PROTECT(call = lang2(out_func, state));
                PROTECT(result = eval(call, out_env));
                if (! isNumeric(result))
                    error("outfun: result of function call must be numeric");
                PROTECT(foo = coerceVector(result, REALSXP));
                if (! isAllFinite(foo))
                    error("outfun returned vector with non-finite element");
                if (LENGTH(foo) != out_dimension)
                    error("outfun return vector length changed from initial");
                for (k = 0; k < out_dimension; k++)
                    REAL(buffer)[k] = REAL(foo)[k];
                UNPROTECT(3);
            }
            break;
        default:
            error("bogus out option\n");
    }
}

