
#include <stddef.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "myutil.h"

// global variables.  Yes, that means we are not cool by some definitions.
static int dimOC = 0;
static int dimNC = 0;
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
static int walker_n = 0;
static double *walker_q = NULL;
static int *walker_j = NULL;

static void transform(double *NCstate, double *OCstate) {
    if (dimOC == 0 || dimNC == 0 || save_origin == NULL || save_basis == NULL)
        error("transform not initialized");
    for (int i = 0; i < dimOC; i++)
        OCstate[i] = save_origin[i];
    for (int j = 0, k = 0; j < dimNC; j++)
        for (int i = 0; i < dimOC; i++, k++)
            OCstate[i] += save_basis[k] * NCstate[j];
}

static double ludfun(double *NCstate) {
    if (dimOC == 0 || save_alpha == NULL)
        error("ludfun not initialized");
    double OCstate[dimOC];
    transform(NCstate, OCstate);
    for (int i = 0; i < dimOC; i++)
        if (OCstate[i] <= 0.0)
            return R_NegInf;
    double result = 0.0;
    for (int i = 0; i < dimOC; i++)
        result += (save_alpha[i] - 1.0) * log(OCstate[i]);
    if (result == R_PosInf)
        error("log unnormalized density +Inf");
    if (R_IsNaN(result))
        error("log unnormalized density NaN");
    return result;
}

static int is_feasible(double *NCstate) {
    if (ncons == 0 || dimNC == 0 || save_amat == NULL | save_bvec == NULL)
        error("is_feasible not initialized");
    double product[ncons];
    for (int j = 0, k = 0; j < dimNC; j++)
        for (int i = 0; i < ncons; i++, k++)
            product[i] += save_amat[k] * NCstate[j];
    for (int i = 0; i < ncons; i++)
        if (product[i] > save_bvec[i])
            return 0;
    return 1;
}

// Walker's method of aliases
//
// sources
//     Devroye (1986).  Non-Uniform Random Variate Generation.  Springer.
//     Kronmal and Peterson (1979).  American Statistician, 33, 214-216.

static void walker_teardown()
{
    if (walker_n == 0 && walker_q == NULL && walker_j == NULL) {
        return;
    } else if (walker_n != 0 && walker_q != NULL && walker_j != NULL) {
        Free(walker_q);
        Free(walker_j);
        walker_j = NULL; 
        walker_q = NULL; 
        walker_n = 0;
        return;
    } else {
        error("walker_teardown: bad state");
    }
}

static void walker_setup(double *p, int n)
{
    if (n <= 0)
        error("walker_setup: n not positive");

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        // Rprintf("p[%d] = %f\n", i + 1, p[i]);
        if (p[i] < 0.0)
            error("walker_setup: p[%d] < 0.0", i + 1);
        sum += p[i];
    }

    if (walker_n == 0 && walker_q == NULL && walker_j == NULL) {
        walker_n = n;
        walker_j = Calloc(n, int);
        walker_q = Calloc(n, double);
    } else if (walker_n != 0 && walker_q != NULL && walker_j != NULL) {
        error("walker_setup: alias tables already set up ???");
    } else {
        error("walker_setup: bad state");
    }

    int *smaller = Calloc(n, int);
    int *greater = Calloc(n, int);
    int nsmaller = 0;
    int ngreater = 0;

    for (int i = 0; i < n; i++) {
        double foo = p[i] / sum * n;
        walker_q[i] = foo;
        if (foo < 1.0) {
            smaller[nsmaller++] = i;
        } else {
            greater[ngreater++] = i;
        }
    }
    while (nsmaller > 0 && ngreater > 0) {
        // Kronmal and Peterson assert that we do not need the test for
        // ngreater > 0 because greater is always nonempty at the bottom
        // of the loop, but this assumes infinite precision arithmetic;
        // can have an index moved from greater to smaller when its q is
        // >= 1 in infinite-precision arithmetic but < 1 n computer arithmetic
        // if ngreater == 0 it must be the case that all ilow in smaller have
        // q[ilow] nearly equal to 1.0, in which case we can consider them
        // all in greater and quit
        int ilow = smaller[nsmaller-1];
        int ihig = greater[ngreater-1];
        walker_j[ilow] = ihig;
        walker_q[ihig] -= (1.0 - walker_q[ilow]);
        if (walker_q[ihig] < 1.0) {
            ngreater--;
            smaller[nsmaller-1] = ihig;
        } else {
            nsmaller--;
        }
    }
    Free(smaller);
    Free(greater);
}

static int walker_rand()
{
    if (walker_n == 0 || walker_q == NULL || walker_j == NULL)
        error("walker_rand: alias tables not properly set up");

    // use runif rather than unif_rand because it guarantees 0 < u < 1
    double u = runif(0.0, 1.0);
    int i = u * walker_n;
    double v = runif(0.0, 1.0);
    if (v <= walker_q[i])
        return i;
    else
        return walker_j[i];
}

// more global variables
static SEXP outfun_func;
static SEXP outfun_env;
static int outfun_result_dimension;

static void outfun_setup(SEXP func, SEXP rho, SEXP state)
{
    if (dimNC == 0 || dimOC == 0)
        error("outfun_setup: not initialized");
    if (LENGTH(state) != dimNC)
        error("outfun_setup: state wrong dimension");
    if (func == R_NilValue) {
        outfun_func = R_NilValue;
        outfun_env = R_NilValue;
        outfun_result_dimension = dimOC;
    } else if (isFunction(func) && isEnvironment(rho)) {
        outfun_func = func;
        outfun_env = rho;
        outfun_result_dimension = LENGTH(eval(lang2(func, state), rho));
    } else if (isFunction(func)) {
        error("out_setup: argument \"rho\" must be environment");
    } else {
        error("out_setup: argument \"func\" must be function or NULL");
    }
}

static void outfun(double *state, double *result)
{
    if (isFunction(outfun_func)) {
        SEXP sexp_state, sexp_result, call, foo;
        PROTECT(sexp_state = allocVector(REALSXP, dimNC));
        for (int i = 0; i < dimNC; i++)
            REAL(sexp_state)[i] = state[i];
        PROTECT(call = lang2(outfun_func, sexp_state));
        PROTECT(sexp_result = eval(call, outfun_env));
        if (! isNumeric(sexp_result))
            error("outfun: result of function call must be numeric");
        PROTECT(foo = coerceVector(sexp_result, REALSXP));
        if (! isAllFinite(foo))
            error("outfun returned vector with non-finite element");
        if (LENGTH(foo) != outfun_result_dimension)
            error("outfun return vector length changed from initial");
        for (i = 0; i < outfun_result_dimension; i++)
            result[i] = REAL(foo)[i];
        UNPROTECT(4);
    } else if (outfun_func == R_NilValue) {
        transform(state, result);
    } else {
        error("outfun not initialized properly");
    }
}


    /* REVISED DOWN TO HERE */

void propose(SEXP state, SEXP proposal, SEXP amat, SEXP bvec,
    double *z, double *smax, double *smin, double *u);



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
    if (! IS_SCALAR(nbatch, INTSXP))
        error("argument \"nbatch\" must be scalar integer");
    if (! IS_SCALAR(blen, INTSXP))
        error("argument \"blen\" must be scalar integer");
    if (! IS_SCALAR(nspac, INTSXP))
        error("argument \"nspac\" must be scalar integer");
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
    if (! IS_SCALAR(mydebug, LGLSXP))
        error("argument \"mydebug\" must be type logical");

    // assign to global variables
    dimOC = LENGTH(param);
    dimNC = LENGTH(initial);

    if (dimNC <= 0)
        error("length(initial) not positive");
    if (dimOC <= dimNC)
        error("must have length(param) > length(initial)");

    if (! isMatrix(basis))
        error("argument \"basis\" must be matrix");
    if (! isMatrix(amat))
        error("argument \"amat\" must be matrix");
    if (! isMatrix(vert))
        error("argument \"vert\" must be matrix");
    if (! isMatrix(cholSigma))
        error("argument \"cholSigma\" must be matrix");

    if (LENGTH(origin) != dimOC)
        error("must have length(param) = length(origin)");
    if (nrows(basis) != dimOC)
        error("must have length(param) = nrow(basis)");
    if (ncols(basis) != dimNC)
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
    if (ncols(amat) != dimNC)
        error("must have ncol(amat) = ncol(basis)");
    if (ncols(vert) != dimNC)
        error("must have ncol(vert) = ncol(basis)");
    if (LENGTH(mu) != dimNC)
        error("must have length(mu) = ncol(basis)");
    if (nrows(cholSigma) != dimNC)
        error("must have nrow(cholSigma) = ncol(basis)");
    if (ncols(cholSigma) != dimNC)
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

    int_nbatch = INTEGER(nbatch)[0];
    int_blen = INTEGER(blen)[0];
    int_nspac = INTEGER(nspac)[0];
    int_debug = LOGICAL(mydebug)[0];

    if (int_nbatch <= 0)
        error("argument \"nbatch\" must be positive");
    if (int_blen <= 0)
        error("argument \"blen\" must be positive");
    if (int_nspac <= 0)
        error("argument \"nspac\" must be positive");

    // assign to global variables
    save_alpha = REAL(param);
    save_origin = REAL(origin);
    save_basis = REAL(basis);
    save_amat = REAL(amat);
    save_bvec = REAL(bvec);
    save_vert = REAL(vert);
    save_mu = REAL(mu);
    save_cholSigma = REAL(cholSigma);
    save_alphaCons = REAL(alphaCons);

    for (int i = 0; i < dimOC; i++)
        if (save_alpha[i] <= 0)
            error("param not all strictly positive");

    double *double_mixprob = REAL(mixprob);
    for (int i = 0; i < 3; i++)
        if (double_mixprob[i] < 0.0)
            error("mixprob not all nonnegative");
    {
        double sum == 0;
        for (int i = 0; i < 3; i++)
            sum += double_mixprob[i];
        if (fabs(sum - 1.0) > 1e-14)
            error("mixprob does not sum to one");
    }
    // Warning: if we ever allow other than length 3 for mixprob this is FUBAR
    double cumsum_mixprob[3];
    cumsum_mixprob[0] = double_mixprob[0];
    cumsum_mixprob[1] = double_mixprob[1] + cumsum_mixprob[9];
    cumsum_mixprob[2] = 1.0;

    if (! is_feasible(REAL(initial)))
        stop("initial does not satisfy constraints");
    double current_log_dens = ludfun(REAL(initial));
    if (! R_finite(current_log_dens))
        stop("ludfun(initial) infinite or NaN");

    double state[dimNC];
    double proposal[dimNC];
    for (int i = 0; i < dimNC; i++)
        state[i] = REAL(initial)[i];

    /* REVISED DOWN TO HERE */

    int dim_out;
    SEXP result, resultnames, path, save_initial, save_final;
    double *batch_buffer;
    SEXP out_buffer;





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

