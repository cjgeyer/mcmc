
#include <stddef.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "myutil.h"

// global variables.  Yes, that means we are not cool by some definitions.
static int d = 0;
static int p = 0;
static int ncons = 0;
static int nvert = 0;
static int nvertdiff = 0;
static double *save_alpha = NULL;
static double *save_origin = NULL;
static double *save_basis = NULL;
static double *save_amat = NULL;
static double *save_bvec = NULL;
static double *save_vert = NULL;
static double *save_vert_diff = NULL;
static double *save_mu = NULL;
static double *save_cholSigma = NULL;
static double *save_alphaCons = NULL;
// alias tables for Walker's method of aliases
static double *p = NULL;
static double *q = NULL;
static int *HL = NULL;

static void transform(double NCstate, double *OCstate) {
    if (d = 0 || p = 0 || save_origin = NULL || save_basis = NULL)
        error("transform not initialized");
    for (int i = 0; i < d; i++)
        OCstate[i] = save_origin[i];
    for (int j = 0, k = 0; j < p; j++)
        for (int i = 0; i < d; i++, k++)
            OCstate[i] += save_basis[k] * NCstate[j];
}

static double lud(double *NCstate) {
    if (d = 0 || save_alpha = NULL)
        error("lud not initialized");
    double OCstate[d];
    transform(NCstate, OCstate);
    for (int i = 0; i < d; i++)
        if (OCstate[i] <= 0.0)
            return R_NegInf;
    double result = 0.0;
    for (int i = 0; i < d; i++)
        result += (save_alpha[i] - 1.0) * log(OCstate[i]);
    if (R_finite(result))
        return result;
    error("log unnormalized density NaN or +Inf");
    return R_NegInf; // Can't get here, of course; avoid GCC warning
}


// ripped off from src/main/random.c in R 3.1.0
static void walker_ProbSampleReplace_setup()
{
    if (nvert = 0 | p = 0 || save_vert = NULL)
        error("walker_ProbSampleReplace_setup not initialized");
    // following two statements are two instead of one
    // to force multiplication before division
    // in C parentheses in expressions do not force order of operations
    nvertdiff = nvert * (nvert - 1);
    nvertdiff /= 2;
    save_vert_diff = Calloc(nvertdiff, double);
    for (int i = 0, k = 0; i < nvert; i++)
        for (int j = i + 1; j < nvert; j++, k++)
            for (int m = 0; m < p; m++)
                save_vert_diff[k + nvertdiff * m] = save_vert[i + nvert * m]
                    - save_vert[j + nvert * m];
    p = Calloc(nvertdiff, double);
    for (int k = 0; k < nvertdiff; k++) {
        double sum = 0.0;
        for (int m = 0; m < p; m++) {
            double foo = save_vert_diff[k + nvertdiff * m];
            sum += foo * foo;
        }
        p[k] = sqrt(sum);
    }
    double sum = 0.0;
    for (int k = 0; k < nvertdiff; k++)
        sum + = p[k];
    for (int k = 0; k < nvertdiff; k++)
        p[k] /= sum;

    // walker_ProbSampleReplace(int n, double *p, int *a, int nans, int *ans)
    // n is size of population (produces result between 1 and n)
    // k is size of sample
    // prob is probability vector (length n)
    // ans is the vector of samples (length k), but we want just one at a time
    // nans = k
    // a is temporary storage (length k) ????

    int *H = HL - 1; // is this undefined ???!!!
    int *L = HL + n;
}

static void
walker_ProbSampleReplace(int n, double *p, int *a, int nans, int *ans)
{
    double rU;
    int i, j, k;
    int *H, *L;

    /* Create the alias tables.
       The idea is that for HL[0] ... L-1 label the entries with q < 1
       and L ... H[n-1] label those >= 1.
       By rounding error we could have q[i] < 1. or > 1. for all entries.
     */
    H = HL - 1; L = HL + n;
    for (i = 0; i < n; i++) {
	q[i] = p[i] * n;
	if (q[i] < 1.) *++H = i; else *--L = i;
    }
    if (H >= HL && L < HL + n) { /* So some q[i] are >= 1 and some < 1 */
	for (k = 0; k < n - 1; k++) {
	    i = HL[k];
	    j = *L;
	    a[i] = j;
	    q[j] += q[i] - 1;
	    if (q[j] < 1.) L++;
	    if(L >= HL + n) break; /* now all are >= 1 */
	}
    }
    for (i = 0; i < n; i++) q[i] += i;

    /* generate sample */
    for (i = 0; i < nans; i++) {
	rU = unif_rand() * n;
	k = (int) rU;
	ans[i] = (rU < q[k]) ? k+1 : a[k]+1;
    }
    if(n > SMALL) {
	Free(HL);
	Free(q);
    }
}

void propose(SEXP state, SEXP proposal, SEXP amat, SEXP bvec,
    double *z, double *smax, double *smin, double *u);


static int out_setup(SEXP func, SEXP rho, SEXP state);

static void outfun(SEXP state, SEXP buffer);

SEXP condir(SEXP param, SEXP initial, SEXP nbatch, SEXP blen, SEXP nspac,
    SEXP origin, SEXP basis, SEXP amat, SEXP bvec, SEXP vert, SEXP mu,
    SEXP cholSigma, SEXP alphaCons, SEXP mixprob, SEXP func2, SEXP rho2,
    SEXP mydebug)
{
    if (! isReal(param))
        error("argument \"param\" must be type double");
    if (! isReal(initial))
        error("argument \"initial\" must be type double");
    if (! isInteger(nbatch))
        error("argument \"nbatch\" must be type integer");
    if (! isInteger(blen))
        error("argument \"blen\" must be type integer");
    if (! isInteger(nspac))
        error("argument \"nspac\" must be type integer");
    if (! isReal(origin))
        error("argument \"origin\" must be type double");
    if (! isReal(basis))
        error("argument \"basis\" must be type double");
    if (! isReal(amat))
        error("argument \"amat\" must be type double");
    if (! isReal(bvec))
        error("argument \"bvec\" must be type double");
    if (! isReal(vert))
        error("argument \"vert\" must be type double");
    if (! isReal(mu))
        error("argument \"mu\" must be type double");
    if (! isReal(cholSigma))
        error("argument \"cholSigma\" must be type double");
    if (! isReal(alphaCons))
        error("argument \"alphaCons\" must be type double");
    if (! isReal(mixprob))
        error("argument \"mixprob\" must be type double");
    if (! (isFunction(func2) | isNull(func2)))
        error("argument \"func2\" must be function or NULL");
    if (! (isEnvironment(rho2) | isNull(rho2)))
        error("argument \"rho2\" must be environment or NULL");
    if (! isLogical(mydebug))
        error("argument \"mydebug\" must be type logical");

    // assign to global variables
    d = LENGTH(param);
    p = LENGTH(initial);

    if (p <= 0)
        error("length(initial) not positive");
    if (d <= p)
        error("must have length(param) > length(initial)");

    if (! isMatrix(basis))
        error("argument \"basis\" must be matrix");
    if (! isMatrix(amat))
        error("argument \"amat\" must be matrix");
    if (! isMatrix(vert))
        error("argument \"vert\" must be matrix");
    if (! isMatrix(cholSigma))
        error("argument \"cholSigma\" must be matrix");

    if (LENGTH(origin) != d)
        error("must have length(param) = length(origin)");
    if (nrows(basis) != d)
        error("must have length(param) = nrow(basis)");
    if (ncols(basis) != p)
        error("must have length(initial) = ncol(basis)");

    // assign to global variables
    ncons = nrows(amat);
    nvert = nrows(vert);

    if (ncons < 2)
        error("nrow(amat) < 2");
    if (nvert < 2)
        error("nrow(vert) < 2");

    if (LENGTH(bvec) != ncons)
        error("must have length(bvec) = nrow(amat)");
    if (ncols(amat) != p)
        error("must have ncol(amat) = ncol(basis)");
    if (ncols(vert) != p)
        error("must have ncol(vert) = ncol(basis)");
    if (LENGTH(mu) != p)
        error("must have length(mu) = ncol(basis)");
    if (nrows(cholSigma) != p)
        error("must have nrow(cholSigma) = ncol(basis)");
    if (ncols(cholSigma) != p)
        error("must have ncol(cholSigma) = ncol(basis)");
    if (LENGTH(alphaCons) != ncons)
        error("must have length(alphaCons) = nrow(amat)");
    if (LENGTH(mixprob) != 3)
        error("must have length(mixprob) = 3");

    if (! isAllFinite(param))
        error("all elements of \"param\" must be finite");
    if (! isAllFinite(initial))
        error("all elements of \"initial\" must be finite");
    if (! isAllFinite(origin))
        error("all elements of \"origin\" must be finite");
    if (! isAllFinite(basis))
        error("all elements of \"basis\" must be finite");
    if (! isAllFinite(amat))
        error("all elements of \"amat\" must be finite");
    if (! isAllFinite(bvec))
        error("all elements of \"bvec\" must be finite");
    if (! isAllFinite(vert))
        error("all elements of \"vert\" must be finite");
    if (! isAllFinite(mu))
        error("all elements of \"mu\" must be finite");
    if (! isAllFinite(cholSigma))
        error("all elements of \"cholSigma\" must be finite");
    if (! isAllFinite(alphaCons))
        error("all elements of \"alphaCons\" must be finite");
    if (! isAllFinite(mixprob))
        error("all elements of \"mixprob\" must be finite");

    save_alpha = REAL(param);
    save_origin = REAL(origin);
    save_basis = REAL(basis);
    save_amat = REAL(amat);
    save_bvec = REAL(bvec);
    save_vert = REAL(vert);
    save_mu = REAL(mu);
    save_cholSigma = REAL(cholSigma);
    save_alphaCons = REAL(alphaCons);

    // still to do
    // check all(param > 0)
    // check all(amat %*% initial <= bvec)
    // check all(origin + basis %*% initial > 0)
    // length(nbatch) == 1 && nbatch > 0
    // length(blen) == 1 && blen > 0
    // length(nspac) == 1 && nspac > 0
    // all(mixprob > 0)
    // sum(mixprob) == 1



    /* REVISED DOWN TO HERE */

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

