
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "myutil.h"
#include "mcmc.h"

static void proposal_setup(SEXP scale, int d);

static void propose(SEXP state, SEXP proposal, double *z);

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
        save_initial, save_final, acceptance_rate_batches;
    double *batch_buffer;
    SEXP out_buffer;

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
    proposal_setup(scale, dim_state);

    dim_out = out_setup(func2, rho2, state);
    batch_buffer = (double *) R_alloc(dim_out, sizeof(double));
    PROTECT(out_buffer = allocVector(REALSXP, dim_out));

     if (! int_debug) {
         PROTECT(result = allocVector(VECSXP, 5));
         PROTECT(resultnames = allocVector(STRSXP, 5));
     } else {
         PROTECT(result = allocVector(VECSXP, 11));
         PROTECT(resultnames = allocVector(STRSXP, 11));
     }
     PROTECT(acceptance_rate = allocVector(REALSXP, 1));
     SET_VECTOR_ELT(result, 0, acceptance_rate);
     PROTECT(path = allocMatrix(REALSXP, dim_out, int_nbatch));
     SET_VECTOR_ELT(result, 1, path);
     PROTECT(save_initial = duplicate(state));
     SET_VECTOR_ELT(result, 2, save_initial);
     /* cannot set final yet because we haven't got it yet
        (final value at end of run).
        See third to last statement of this function. */
     PROTECT(acceptance_rate_batches = allocVector(REALSXP, int_nbatch));
     SET_VECTOR_ELT(result, 4, acceptance_rate_batches);
     UNPROTECT(4);
     SET_STRING_ELT(resultnames, 0, mkChar("accept"));
     SET_STRING_ELT(resultnames, 1, mkChar("batch"));
     SET_STRING_ELT(resultnames, 2, mkChar("initial"));
     SET_STRING_ELT(resultnames, 3, mkChar("final"));
     SET_STRING_ELT(resultnames, 4, mkChar("accept.batch"));
     if (int_debug) {
         SEXP spath, ppath, gpath, upath, zpath, apath;
         int nn = int_nbatch * int_blen * int_nspac;
         PROTECT(spath = allocMatrix(REALSXP, dim_state, nn));
         SET_VECTOR_ELT(result, 5, spath);
         PROTECT(ppath = allocMatrix(REALSXP, dim_state, nn));
         SET_VECTOR_ELT(result, 6, ppath);
         PROTECT(gpath = allocVector(REALSXP, nn));
         SET_VECTOR_ELT(result, 7, gpath);
         PROTECT(upath = allocVector(REALSXP, nn));
         SET_VECTOR_ELT(result, 8, upath);
         PROTECT(zpath = allocMatrix(REALSXP, dim_state, nn));
         SET_VECTOR_ELT(result, 9, zpath);
         PROTECT(apath = allocVector(LGLSXP, nn));
         SET_VECTOR_ELT(result, 10, apath);
         UNPROTECT(6);
         SET_STRING_ELT(resultnames, 5, mkChar("current"));
         SET_STRING_ELT(resultnames, 6, mkChar("proposal"));
         SET_STRING_ELT(resultnames, 7, mkChar("log.green"));
         SET_STRING_ELT(resultnames, 8, mkChar("u"));
         SET_STRING_ELT(resultnames, 9, mkChar("z"));
         SET_STRING_ELT(resultnames, 10, mkChar("debug.accept"));
     }
     namesgets(result, resultnames);
     UNPROTECT(1);

     GetRNGstate();

     current_log_dens = logh(func1, state, rho1);
     if (current_log_dens == R_NegInf)
         error("log unnormalized density -Inf at initial state");

     for (int ibatch = 0, k = 0; ibatch < int_nbatch; ibatch++) {

         double acceptances_this_batch = 0.0;
         double tries_this_batch = 0.0;

         for (int i = 0; i < dim_out; i++)
             batch_buffer[i] = 0.0;

         for (int jbatch = 0; jbatch < int_blen; jbatch++) {

             double proposal_log_dens;

             for (int ispac = 0; ispac < int_nspac; ispac++) {

                 int accept;
                 double u = -1.0; /* impossible return from unif_rand() */
                 double z[dim_state]; /* buffer for output of norm_rand() */

                 /* Note: should never happen! */
                 if (current_log_dens == R_NegInf)
                     error("log density -Inf at current state");

                 propose(state, proposal, z);

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
                     SEXP spath = VECTOR_ELT(result, 5);
                     SEXP ppath = VECTOR_ELT(result, 6);
                     SEXP gpath = VECTOR_ELT(result, 7);
                     SEXP upath = VECTOR_ELT(result, 8);
                     SEXP zpath = VECTOR_ELT(result, 9);
                     SEXP apath = VECTOR_ELT(result, 10);
                     for (int lj = 0; lj < dim_state; lj++) {
                         REAL(spath)[lbase + lj] = REAL(state)[lj];
                         REAL(ppath)[lbase + lj] = REAL(proposal)[lj];
                         REAL(zpath)[lbase + lj] = z[lj];
                     }
                     REAL(gpath)[l] = proposal_log_dens - current_log_dens;
                     if (u == -1.0)
                         REAL(upath)[l] = NA_REAL;
                     else
                         REAL(upath)[l] = u;
                     LOGICAL(apath)[l] = accept;
                 }

                 if (accept) {
                     for (int jj = 0; jj < dim_state; jj++)
                         REAL(state)[jj] = REAL(proposal)[jj];
                     current_log_dens = proposal_log_dens;
                     acceptances++;
                     acceptances_this_batch++;
                 }
                 tries++;
                 tries_this_batch++;
             } /* end of inner loop (one iteration) */

             outfun(state, out_buffer);
             for (int j = 0; j < dim_out; j++)
                 batch_buffer[j] += REAL(out_buffer)[j];

         } /* end of middle loop (one batch) */

         for (int j = 0; j < dim_out; j++, k++)
             REAL(path)[k] = batch_buffer[j] / int_blen;

         REAL(acceptance_rate_batches)[ibatch] =
             acceptances_this_batch / tries_this_batch;

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
            scale_factor = (double *) R_alloc(d * d, sizeof(double));
            for (int i = 0; i < d * d; i++)
                scale_factor[i] = REAL(foo)[i];
            scale_option = FULL;
        } else {
            error("dimensions of \"scale\" matrix not d by d");
        }
        UNPROTECT(1);
    } else if (LENGTH(foo) == d) {
        scale_factor = (double *) R_alloc(d, sizeof(double));
        for (int i = 0; i < d; i++)
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

static void propose(SEXP state, SEXP proposal, double *z)
{
    int d = state_dimension;

    if (scale_option == 0)
        error("attempt to call propose without setup");

    if (LENGTH(state) != d || LENGTH(proposal) != d)
        error("State or proposal length different from initialization\n");

    for (int j = 0; j < d; j++)
        z[j] = norm_rand();

    switch (scale_option) {
        case CONSTANT:
            for (int j = 0; j < d; j++)
                REAL(proposal)[j] = REAL(state)[j]
                    + scale_factor[0] * z[j];
            break;
        case DIAGONAL:
            for (int j = 0; j < d; j++)
                REAL(proposal)[j] = REAL(state)[j]
                    + scale_factor[j] * z[j];
            break;
        case FULL:
            for (int j = 0; j < d; j++)
                REAL(proposal)[j] = REAL(state)[j];

            for (int i = 0, k = 0; i < d; i++) {
                double u = z[i];
                for (int j = 0; j < d; j++)
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
        if (LENGTH(func) != out_state_dimension)
            error("is.logical(outfun) & (length(outfun) != length(initial))");
        out_option = OUT_INDEX;
        out_index = (int *) R_alloc(out_state_dimension, sizeof(int));
        out_dimension = 0;
        for (int i = 0; i < out_state_dimension; i++) {
            out_index[i] = LOGICAL(func)[i];
            out_dimension += out_index[i];
        }
    } else if (isNumeric(func)) {
        SEXP foo;
        int foopos = 0;
        int fooneg = 0;
        PROTECT(foo = coerceVector(func, REALSXP));
        int foolen = LENGTH(foo);
        for (int i = 0; i < foolen; i++) {
            double foodble = REAL(foo)[i];
            if (ISNAN(foodble))
                error("NA or NaN index for outfun");
            if (! R_FINITE(foodble))
                error("-Inf or Inf index for outfun");
            int fooint = foodble;
            int fooabs = fooint >= 0 ? fooint : (- fooint);

            if (fooint == 0)
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
            for (int i = 0; i < out_state_dimension; i++)
                out_index[i] = FALSE;
            for (int i = 0; i < foolen; i++) {
                 int fooint = REAL(foo)[i];
                 out_index[fooint - 1] = TRUE;
            }
        } else /* (fooneg > 0) */ {
            for (int i = 0; i < out_state_dimension; i++)
                out_index[i] = TRUE;
            for (int i = 0; i < foolen; i++) {
                 int fooint = REAL(foo)[i];
                 int fooabs = (- fooint);
                 out_index[fooabs - 1] = FALSE;
            }
        }
        out_dimension = 0;
        for (int i = 0; i < out_state_dimension; i++)
            out_dimension += out_index[i];
        UNPROTECT(1);
    } else {
        error("outfun must be NULL, a function, a numeric vector,"
            " or a logical vector");
    }
    return out_dimension;
}

static void outfun(SEXP state, SEXP buffer)
{
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
            for (int j = 0; j < out_state_dimension; j++)
                REAL(buffer)[j] = REAL(state)[j];
            break;
        case OUT_INDEX:
            for (int j = 0, k = 0; j < out_state_dimension; j++)
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
                for (int k = 0; k < out_dimension; k++)
                    REAL(buffer)[k] = REAL(foo)[k];
                UNPROTECT(3);
            }
            break;
        default:
            error("bogus out option\n");
    }
}

