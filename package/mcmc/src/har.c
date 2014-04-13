
/*
*
* mcmc and MCMC package for R
* Copyright (c) 2014 Charles J. Geyer
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

void propose(SEXP state, SEXP proposal, SEXP amat, SEXP bvec,
    double *z, double *smax, double *smin, double *u);

static double logh(SEXP func, SEXP state, SEXP rho);

static int out_setup(SEXP func, SEXP rho, SEXP state);

static void outfun(SEXP state, SEXP buffer);

SEXP har(SEXP func1, SEXP initial, SEXP nbatch, SEXP blen, SEXP nspac,
    SEXP amat, SEXP bvec, SEXP func2, SEXP debug, SEXP rho1, SEXP rho2)
{
    int int_nbatch, int_blen, int_nspac, int_debug;
    SEXP state, proposal;
    int dim_state, dim_out;
    SEXP result, resultnames, path, save_initial, save_final;
    double *batch_buffer;
    SEXP out_buffer;

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
    if (! isNumeric(amat))
        error("argument \"amat\" must be numeric");
    if (! isNumeric(bvec))
        error("argument \"bvec\" must be numeric");

    if (! isLogical(debug))
        error("argument \"debug\" must be logical");

    int_nbatch = getScalarInteger(nbatch, "nbatch");
    int_blen = getScalarInteger(blen, "blen");
    int_nspac = getScalarInteger(nspac, "nspac");

    int_debug = getScalarLogical(debug, "debug");

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

    dim_out = out_setup(func2, rho2, state);
    batch_buffer = (double *) R_alloc(dim_out, sizeof(double));
    PROTECT(out_buffer = allocVector(REALSXP, dim_out));

     if (! int_debug) {
         PROTECT(result = allocVector(VECSXP, 3));
         PROTECT(resultnames = allocVector(STRSXP, 3));
     } else {
         PROTECT(result = allocVector(VECSXP, 11));
         PROTECT(resultnames = allocVector(STRSXP, 11));
     }
     PROTECT(path = allocMatrix(REALSXP, dim_out, int_nbatch));
     SET_VECTOR_ELT(result, 0, path);
     PROTECT(save_initial = duplicate(state));
     SET_VECTOR_ELT(result, 1, save_initial);
     UNPROTECT(2);
     SET_STRING_ELT(resultnames, 0, mkChar("batch"));
     SET_STRING_ELT(resultnames, 1, mkChar("initial"));
     SET_STRING_ELT(resultnames, 2, mkChar("final"));
     if (int_debug) {
         SEXP spath, ppath, zpath, u1path, u2path, s1path, s2path, gpath;
         int nn = int_nbatch * int_blen * int_nspac;
         PROTECT(spath = allocMatrix(REALSXP, dim_state, nn));
         SET_VECTOR_ELT(result, 3, spath);
         PROTECT(ppath = allocMatrix(REALSXP, dim_state, nn));
         SET_VECTOR_ELT(result, 4, ppath);
         PROTECT(zpath = allocMatrix(REALSXP, dim_state, nn));
         SET_VECTOR_ELT(result, 5, zpath);
         PROTECT(u1path = allocVector(REALSXP, nn));
         SET_VECTOR_ELT(result, 6, u1path);
         PROTECT(u2path = allocVector(REALSXP, nn));
         SET_VECTOR_ELT(result, 7, u2path);
         PROTECT(s1path = allocVector(REALSXP, nn));
         SET_VECTOR_ELT(result, 8, s1path);
         PROTECT(s2path = allocVector(REALSXP, nn));
         SET_VECTOR_ELT(result, 9, s2path);
         PROTECT(gpath = allocVector(REALSXP, nn));
         SET_VECTOR_ELT(result, 10, gpath);
         UNPROTECT(8);
         SET_STRING_ELT(resultnames, 3, mkChar("current"));
         SET_STRING_ELT(resultnames, 4, mkChar("proposal"));
         SET_STRING_ELT(resultnames, 5, mkChar("z"));
         SET_STRING_ELT(resultnames, 6, mkChar("u1"));
         SET_STRING_ELT(resultnames, 7, mkChar("u2"));
         SET_STRING_ELT(resultnames, 8, mkChar("s1"));
         SET_STRING_ELT(resultnames, 9, mkChar("s2"));
         SET_STRING_ELT(resultnames, 10, mkChar("log.green"));
     }
     namesgets(result, resultnames);
     UNPROTECT(1);

     GetRNGstate();

     current_log_dens = logh(func1, state, rho1);
     if (current_log_dens == R_NegInf)
         error("log unnormalized density -Inf at initial state");

     for (int ibatch = 0, k = 0; ibatch < int_nbatch; ibatch++) {

         int j;

         for (int i = 0; i < dim_out; i++)
             batch_buffer[i] = 0.0;

         for (int jbatch = 0; jbatch < int_blen; jbatch++) {

             double proposal_log_dens;

             for (int ispac = 0; ispac < int_nspac; ispac++) {

                 /* Note: should never happen! */
                 if (current_log_dens == R_NegInf)
                     error("log density -Inf at current state");

                 double u1 = R_NaReal;
                 double u2 = R_NaReal;
                 double smax = R_NaReal;
                 double smin = R_NaReal;
                 double z[dim_state];

                 propose(state, proposal, amat, bvec, z, &smax, &smin, &u1);

                 proposal_log_dens = logh(func1, proposal, rho1);

                 int accept = FALSE;
                 if (proposal_log_dens != R_NegInf) {
                     if (proposal_log_dens >= current_log_dens) {
                         accept = TRUE;
                     } else {
                         double green = exp(proposal_log_dens
                             - current_log_dens);
                         u2 = unif_rand();
                         accept = u2 < green;
                     }
                 }

                 if (int_debug) {
                     int l = ispac + int_nspac * (jbatch + int_blen * ibatch);
                     int lbase = l * dim_state;
                     SEXP spath = VECTOR_ELT(result, 3);
                     SEXP ppath = VECTOR_ELT(result, 4);
                     SEXP zpath = VECTOR_ELT(result, 5);
                     SEXP u1path = VECTOR_ELT(result, 6);
                     SEXP u2path = VECTOR_ELT(result, 7);
                     SEXP s1path = VECTOR_ELT(result, 8);
                     SEXP s2path = VECTOR_ELT(result, 9);
                     SEXP gpath = VECTOR_ELT(result, 10);
                     for (int lj = 0; lj < dim_state; lj++) {
                         REAL(spath)[lbase + lj] = REAL(state)[lj];
                         REAL(ppath)[lbase + lj] = REAL(proposal)[lj];
                         REAL(zpath)[lbase + lj] = z[lj];
                     }
                     REAL(u1path)[l] = u1;
                     REAL(u2path)[l] = u2;
                     REAL(s1path)[l] = smin;
                     REAL(s2path)[l] = smax;
                     REAL(gpath)[l] = proposal_log_dens - current_log_dens;
                 }

                 if (accept) {
                     for (int jj = 0; jj < dim_state; jj++)
                         REAL(state)[jj] = REAL(proposal)[jj];
                     current_log_dens = proposal_log_dens;
                 }
             } /* end of inner loop (one iteration) */

             outfun(state, out_buffer);
             for (j = 0; j < dim_out; j++)
                 batch_buffer[j] += REAL(out_buffer)[j];

         } /* end of middle loop (one batch) */

         for (j = 0; j < dim_out; j++, k++)
             REAL(path)[k] = batch_buffer[j] / int_blen;

     } /* end of outer loop */

     PutRNGstate();

     PROTECT(save_final = coerceVector(state, REALSXP));
     SET_VECTOR_ELT(result, 2, save_final);

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

void propose(SEXP state, SEXP proposal, SEXP amat, SEXP bvec,
    double *z, double *smax_out, double *smin_out, double *u_out)
{
    int d = LENGTH(state);
    int n = nrows(amat);
    if (LENGTH(proposal) != d)
        error("propose: length(state) != length(proposal)");
    if (ncols(amat) != d)
        error("propose: ncol(amat) != length(state)");
    if (n != LENGTH(bvec))
        error("propose: nrow(amat) != length(bvec)");
    double *a = REAL(amat);
    double *b = REAL(bvec);
    double *x = REAL(state);

    for (int i = 0; i < d; i++) {
        z[i] = norm_rand();
    }

    double smax = R_PosInf;
    double smin = R_NegInf;

    for (int i = 0; i < n; i++) {

        double ax = 0.0;
        double az = 0.0;
        for (int j = 0; j < d; j++) {
            ax += a[i + j * n] * x[j];
            az += a[i + j * n] * z[j];
        }
        double bound = (b[i] - ax) / az;
        if (az > 0 && bound < smax)
                smax = bound;
        if (az < 0 && bound > smin)
                smin = bound;
    }

    double u = unif_rand();

    for (int i = 0; i < d; i++)
        REAL(proposal)[i] = x[i] + (u * smin + (1.0 - u) * smax) * z[i];

    *smax_out = smax;
    *smin_out = smin;
    *u_out = u;
}

static SEXP out_func;
static SEXP out_env;
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
    } else {
        error("out_setup: argument \"func\" must be function or NULL");
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

