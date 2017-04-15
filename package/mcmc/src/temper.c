
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "myutil.h"
#include "mcmc.h"

#ifdef BLEAT
#include <stdio.h>
#endif /* BLEAT */

static void propose(SEXP coproposal, SEXP proposal, SEXP scale, double *z);

static double logh(SEXP func, SEXP state, SEXP rho);

static SEXP outfun(SEXP func, SEXP state, SEXP rho);

static void check_valid_scale(SEXP scale, int i, int ncomp, int nx);

SEXP temper(SEXP func1, SEXP initial, SEXP neighbors, SEXP nbatch,
    SEXP blen, SEXP nspac, SEXP scale, SEXP func2, SEXP debug,
    SEXP parallel, SEXP rho1, SEXP rho2)
{
    if (! isFunction(func1))
        error("argument \"func1\" must be function");
    if (! isEnvironment(rho1))
        error("argument \"rho1\" must be environment");

    int is_parallel = getScalarLogical(parallel, "parallel");
    int is_debug = getScalarLogical(debug, "debug");

    if (! isLogical(neighbors))
        error("argument \"neighbors\" must be logical");
    if (! isMatrix(neighbors))
        error("argument \"neighbors\" must be matrix");
    if (nrows(neighbors) != ncols(neighbors))
        error("argument \"neighbors\" must have same row and column dimension");
    int ncomp = nrows(neighbors);
    if (ncomp <= 1)
        error("must have at least two components");
    for (int i = 0; i < ncomp; i++)
        for (int j = 0; j < ncomp; j++)
            if (LOGICAL(neighbors)[i + ncomp * j] != LOGICAL(neighbors)[j + ncomp * i])
                error("argument \"neighbors\" must be symmetric matrix");

    if (! isReal(initial))
        error("argument \"initial\" must be real");
    if (! isAllFinite(initial))
        error("argument \"initial\" must have all elements finite");
    int nx;
    if (is_parallel) {
        if (! isMatrix(initial))
            error("argument \"initial\" must be matrix in parallel case");
        if (nrows(initial) != ncomp)
            error("row dim of args \"initial\" and \"neighbors\" must be same in parallel case");
        nx = ncols(initial);
    } else /* serial */ {
        if (! (LENGTH(initial) > 1))
            error("argument \"initial\" must have length > 1 in serial case");
        double reali = REAL(initial)[0];
        int i = reali;
        if (i != reali)
            error("1st elem of argument \"initial\" must be integer in serial case");
        if (i <= 0 || i > ncomp)
            error("1st elem of argument \"initial\" must be in 1, ..., k in serial case");
        nx = LENGTH(initial) - 1;
    }

    int int_nbatch = getScalarInteger(nbatch, "nbatch");
    if (int_nbatch <= 0)
        error("argument \"nbatch\" must be positive");

    int int_blen = getScalarInteger(blen, "blen");
    if (int_blen <= 0)
        error("argument \"blen\" must be positive");

    int int_nspac = getScalarInteger(nspac, "nspac");
    if (int_nspac <= 0)
        error("argument \"nspac\" must be positive");

    if (isNewList(scale)) {
        if (LENGTH(scale) != ncomp)
            error("argument \"scale\" must have length k if list");
        for (int i = 0; i < ncomp; i++) {
            SEXP fred = VECTOR_ELT(scale, i);
            check_valid_scale(fred, i, ncomp, nx);
        }
    } else /* scale not list */ {
        check_valid_scale(scale, -1, ncomp, nx);
    }

    int no_outfun = isNull(func2);
    if (! no_outfun) {
        if (! isFunction(func2))
            error("argument \"outfun\" must be function");
        if (! isEnvironment(rho2))
            error("argument \"rho2\" must be environment");
    }

    double current_log_dens[ncomp];
    if (is_parallel) {
        SEXP fred;
        PROTECT(fred = allocVector(REALSXP, nx + 1));
        for (int i = 0; i < ncomp; i++) {
            REAL(fred)[0] = i + 1;
            for (int j = 0; j < nx; j++)
                REAL(fred)[j + 1] = REAL(initial)[i + ncomp * j];
            current_log_dens[i] = logh(func1, fred, rho1);
#ifdef BLATHER
            fprintf(stderr, "current_log_dens[%d] = %e\n", i, current_log_dens[i]);
            for (int j = 0; j < nx; j++)
                fprintf(stderr, "    state[%d, %d] = %e\n",
                    i, j, REAL(initial)[i + ncomp * j]);
            for (int j = 0; j <= nx; j++)
                fprintf(stderr, "    fred[%d] = %e\n", j, REAL(fred)[j]);
            fprintf(stderr, "    logh(func1, fred, rho1)) = %e\n",
                        logh(func1, fred, rho1));
#endif /* BLATHER */
            if (current_log_dens[i] == R_NegInf)
                error("log unnormalized density -Inf at initial state");
        }
        UNPROTECT(1);
    } else /* serial */ {
        for (int i = 0; i < ncomp; i++)
            current_log_dens[i] = R_NaN;
        int i = REAL(initial)[0] - 1;
        current_log_dens[i] = logh(func1, initial, rho1);
        if (current_log_dens[i] == R_NegInf)
            error("log unnormalized density -Inf at initial state");
    }

    /* at this point, all arguments have been checked for validity */

    // current_log_dens saves info that cuts in half
    // the number of invocations of log unnormalized density

    SEXP state, proposal, coproposal;
    PROTECT(state = duplicate(initial));
    PROTECT(proposal = allocVector(REALSXP, nx + 1));
    PROTECT(coproposal = allocVector(REALSXP, nx + 1));

    int nout;
    if (no_outfun) {
        nout = is_parallel ? nx * ncomp : nx;
    } else /* has outfun */ {
        nout = LENGTH(outfun(func2, state, rho2));
    }

    int niter = int_nbatch * int_blen * int_nspac;

    // TO DO LIST
    //
    // regular output (if serial)
    //
    //     batch      nbatch x nout matrix
    //     ibatch     nbatch x ncomp matrix
    //     acceptx    vector of length ncomp
    //     accepti    ncomp x ncomp matrix
    //     initial    copy of initial state
    //     final      final state
    //
    // regular output (if parallel)
    //
    //     batch      (if outfun) nbatch x nout matrix
    //                (if no outfun) nbatch x ncomp x nx array
    //     acceptx    vector of length ncomp
    //     accepti    ncomp x ncomp matrix
    //     initial    copy of initial state
    //     final      final state
    //
    // debug output (if not parallel)
    //
    //     which         vector of length niter (TRUE if within-component)
    //     unif_which    vector of length niter (uniform for deciding which)
    //     state         niter x (nx + 1) matrix (state before)
    //     proposal      niter x (nx + 1) matrix
    //     log_hastings  niter vector
    //     unif_hastings niter vector (uniform for deciding acceptance)
    //     acceptd       niter vector (TRUE if accept)
    //     norm          niter x nx matrix (std normals)
    //     unif_choose   niter vector (uniform for choosing neighbor)
    //
    // debug output (if parallel)
    //
    //     which         vector of length niter (TRUE if within-component)
    //     unif_which    vector of length niter (uniform for deciding which)
    //     state         niter x ncomp x nx array (state before)
    //     coproposal    niter x (nx + 1) matrix
    //     proposal      niter x (nx + 1) matrix
    //     log_hastings  niter vector
    //     unif_hastings niter vector (uniform for deciding acceptance)
    //     acceptd       niter vector (TRUE if accept)
    //     norm          niter x nx matrix (std normals)
    //     unif_choose   niter x 2 matrix (uniforms for choosing components
    //                       to update)
    //
    //     for within-component move coproposal and proposal have natural
    //         meaning
    //     for swap move we have 2 coproposals (i, x_i) and (j, x_j)
    //         and 2 proposals (i, x_j) and (j, x_i) but since the information
    //         here is quite redundant we just store (i, x_i) in "coproposal"
    //         and (j, x_j) in "proposal" -- the checker can figure it out

    int len_result_regular = is_parallel ? 5 : 6;
    int len_result_debug = is_parallel ? 10 : 9;
    int len_result = len_result_regular;
    len_result += is_debug ? len_result_debug : 0;

#ifdef BLEAT
    fprintf(stderr, "len_result = %d\n", len_result);
#endif /* BLEAT */

    SEXP result, resultnames, acceptx, accepti, batch, ibatch,
        save_initial, save_final, debug_which, debug_unif_which,
        debug_state, debug_coproposal, debug_proposal,
        debug_log_hastings, debug_unif_hastings, debug_acceptd,
        debug_norm, debug_unif_choose;

    PROTECT(result = allocVector(VECSXP, len_result));
    PROTECT(resultnames = allocVector(STRSXP, len_result));
    namesgets(result, resultnames);
    UNPROTECT(1);

    if (no_outfun && is_parallel)
        PROTECT(batch = alloc3DArray(REALSXP, int_nbatch, ncomp, nx));
    else
        PROTECT(batch = allocMatrix(REALSXP, int_nbatch, nout));
    SET_VECTOR_ELT(result, 0, batch);
    SET_STRING_ELT(resultnames, 0, mkChar("batch"));
    UNPROTECT(1);

    PROTECT(acceptx = allocVector(REALSXP, ncomp));
    SET_VECTOR_ELT(result, 1, acceptx);
    SET_STRING_ELT(resultnames, 1, mkChar("acceptx"));
    UNPROTECT(1);

    PROTECT(accepti = allocMatrix(REALSXP, ncomp, ncomp));
    SET_VECTOR_ELT(result, 2, accepti);
    SET_STRING_ELT(resultnames, 2, mkChar("accepti"));
    UNPROTECT(1);

    PROTECT(save_initial = duplicate(initial));
    SET_VECTOR_ELT(result, 3, save_initial);
    SET_STRING_ELT(resultnames, 3, mkChar("initial"));
    UNPROTECT(1);

    SET_STRING_ELT(resultnames, 4, mkChar("final"));
    // at end need to duplicate state as save_final and copy to result[4]

    if (! is_parallel) {
        PROTECT(ibatch = allocMatrix(REALSXP, int_nbatch, ncomp));
        SET_VECTOR_ELT(result, 5, ibatch);
        SET_STRING_ELT(resultnames, 5, mkChar("ibatch"));
        UNPROTECT(1);
    }

    if (is_debug) {

        PROTECT(debug_which = allocVector(LGLSXP, niter));
        SET_VECTOR_ELT(result, len_result_regular + 0, debug_which);
        SET_STRING_ELT(resultnames, len_result_regular + 0, mkChar("which"));
        UNPROTECT(1);

        PROTECT(debug_unif_which = allocVector(REALSXP, niter));
        SET_VECTOR_ELT(result, len_result_regular + 1, debug_unif_which);
        SET_STRING_ELT(resultnames, len_result_regular + 1,
            mkChar("unif.which"));
        UNPROTECT(1);

        if (is_parallel)
            PROTECT(debug_state = alloc3DArray(REALSXP, niter, ncomp, nx));
        else
            PROTECT(debug_state = allocMatrix(REALSXP, niter, nx + 1));
        SET_VECTOR_ELT(result, len_result_regular + 2, debug_state);
        SET_STRING_ELT(resultnames, len_result_regular + 2, mkChar("state"));
        UNPROTECT(1);

        PROTECT(debug_log_hastings = allocVector(REALSXP, niter));
        SET_VECTOR_ELT(result, len_result_regular + 3, debug_log_hastings);
        SET_STRING_ELT(resultnames, len_result_regular + 3,
            mkChar("log.hastings"));
        UNPROTECT(1);

        PROTECT(debug_unif_hastings = allocVector(REALSXP, niter));
        SET_VECTOR_ELT(result, len_result_regular + 4, debug_unif_hastings);
        SET_STRING_ELT(resultnames, len_result_regular + 4,
            mkChar("unif.hastings"));
        UNPROTECT(1);

        PROTECT(debug_proposal = allocMatrix(REALSXP, niter, nx + 1));
        SET_VECTOR_ELT(result, len_result_regular + 5, debug_proposal);
        SET_STRING_ELT(resultnames, len_result_regular + 5,
            mkChar("proposal"));
        UNPROTECT(1);

        PROTECT(debug_acceptd = allocVector(LGLSXP, niter));
        SET_VECTOR_ELT(result, len_result_regular + 6, debug_acceptd);
        SET_STRING_ELT(resultnames, len_result_regular + 6,
            mkChar("acceptd"));
        UNPROTECT(1);

        PROTECT(debug_norm = allocMatrix(REALSXP, niter, nx));
        SET_VECTOR_ELT(result, len_result_regular + 7, debug_norm);
        SET_STRING_ELT(resultnames, len_result_regular + 7,
            mkChar("norm"));
        UNPROTECT(1);

        if (is_parallel)
            PROTECT(debug_unif_choose = allocMatrix(REALSXP, niter, 2));
        else
            PROTECT(debug_unif_choose = allocVector(REALSXP, niter));
        SET_VECTOR_ELT(result, len_result_regular + 8, debug_unif_choose);
        SET_STRING_ELT(resultnames, len_result_regular + 8,
            mkChar("unif.choose"));
        UNPROTECT(1);

        if (is_parallel) {
            PROTECT(debug_coproposal = allocMatrix(REALSXP, niter, nx + 1));
            SET_VECTOR_ELT(result, len_result_regular + 9, debug_coproposal);
            SET_STRING_ELT(resultnames, len_result_regular + 9,
                mkChar("coproposal"));
            UNPROTECT(1);
        }
    }

    // at this point, entire output structure (SEXP result) is set up, except
    // for aforementioned need to duplicate final state and put in result[4]

    GetRNGstate();

    // need buffers for acceptance rate(s)

    double acceptx_numer[ncomp];
    double acceptx_denom[ncomp];
    double accepti_numer[ncomp][ncomp];
    double accepti_denom[ncomp][ncomp];
    for (int i = 0; i < ncomp; i++) {
        acceptx_numer[i] = 0;
        acceptx_denom[i] = 0;
        for (int j = 0; j < ncomp; j++) {
            accepti_numer[i][j] = 0;
            accepti_denom[i][j] = 0;
        }
    }

    // need neighbor counts and neighbors
    // note: the_neighbors uses zero-origin indexing both for indexing
    //     and values

    double n_neighbors[ncomp];
    for (int i = 0; i < ncomp; i++) {
        n_neighbors[i] = 0;
        for (int j = 0; j < ncomp; j++)
            n_neighbors[i] += LOGICAL(neighbors)[i + ncomp * j];
    }

    double the_neighbors[ncomp][ncomp];
    for (int i = 0; i < ncomp; i++) {
        for (int j = 0, k = 0; j < ncomp; j++) {
            if (LOGICAL(neighbors)[i + ncomp * j])
                the_neighbors[i][k++] = j;
        }
    }

    // need buffers for batch means

    double batch_buff[nout];
    double ibatch_buff[ncomp];

    for (int kbatch = 0, iiter = 0; kbatch < int_nbatch; kbatch++) {

        for (int i = 0; i < nout; i++)
            batch_buff[i] = 0.0;
        for (int i = 0; i < ncomp; i++)
            ibatch_buff[i] = 0.0;

        for (int jbatch = 0; jbatch < int_blen; jbatch++) {

            for (int ispac = 0; ispac < int_nspac; ispac++, iiter++) {

#ifdef EXTRA_CHECK
#ifdef WOOF
                fprintf(stderr, "Check for validity of current_log_dens at top of inner loop\n");
#endif /* WOOF */
                if (is_parallel) {
                    for (int i = 0; i < ncomp; i++) {
                        REAL(proposal)[0] = i + 1;
                        for (int j = 0; j < nx; j++)
                            REAL(proposal)[j + 1] = REAL(state)[i + ncomp * j];
#ifdef BLATHER
            fprintf(stderr, "current_log_dens[%d] = %e\n", i, current_log_dens[i]);
            for (int j = 0; j < nx; j++)
                fprintf(stderr, "    state[%d, %d] = %e\n",
                    i, j, REAL(state)[i + ncomp * j]);
            for (int j = 0; j <= nx; j++)
                fprintf(stderr, "    proposal[%d] = %e\n", j, REAL(proposal)[j]);
            fprintf(stderr, "    logh(func1, proposal, rho1)) = %e\n",
                        logh(func1, proposal, rho1));
#endif /* BLATHER */
                        if (current_log_dens[i] != logh(func1, proposal, rho1))
                            error("current_log_dens[%d] bogus\n", i);
                    }
                } else /* serial */ {
                    for (int j = 0; j <= nx; j++)
                        REAL(proposal)[j] = REAL(state)[j];
                    int i = REAL(proposal)[0] - 1;
                    double my_actual_logh = logh(func1, proposal, rho1);
                    double my_stored_logh = current_log_dens[i];
#ifdef WOOF
                    fprintf(stderr, "icomp = %d, stored logh = %e, actual logh = %e\n",
                        i + 1, my_stored_logh, my_actual_logh);
#endif /* WOOF */
                    if (my_stored_logh != my_actual_logh)
                        error("current_log_dens[%d] bogus\n", i);
                }
#endif /* EXTRA_CHECK */

                if (is_debug) {
                    int len_state = is_parallel ? ncomp * nx : nx + 1;
                    for (int j = 0; j < len_state; j++)
                        REAL(debug_state)[iiter + niter * j] = REAL(state)[j];
                }

                double my_unif_which = unif_rand();
                int my_which = my_unif_which < 0.5;

                if (is_debug) {
                    LOGICAL(debug_which)[iiter] = my_which;
                    REAL(debug_unif_which)[iiter] = my_unif_which;
                }

                if (my_which) /* within-component update */ {

                    if (is_parallel) {

                        // note: my_i and my_j are 1-origin indexing (for R)
                        //         go from 1, ..., ncomp
                        // everything else 0-origin indexing (for C)

                        double unif_choose = unif_rand();

                        int my_i = trunc(ncomp * unif_choose) + 1;
                        if (my_i > ncomp) my_i--;
                        if (my_i <= 0 || my_i > ncomp)
                            error("Can't happen: my_i out of range");

                        REAL(coproposal)[0] = my_i;
                        for (int j = 0; j < nx; j++)
                            REAL(coproposal)[j + 1] =
                                REAL(state)[(my_i - 1) + ncomp * j];

                        double z[nx];
                        propose(coproposal, proposal, scale, z);

                        double my_coproposal_log_dens =
                            current_log_dens[my_i - 1];

#ifdef EXTRA_CHECK
                        if (my_coproposal_log_dens != logh(func1, coproposal, rho1)) {
                            fprintf(stderr, "with-in component update (parallel)\n");
                            error("saving logh didn't work right (coproposal)");
                        }
#endif /* EXTRA_CHECK */
#ifdef BLEAT
                        if (my_coproposal_log_dens == R_NegInf) {
                            fprintf(stderr, "Oopsie #1!\n");
                            fprintf(stderr, "    my_i = %d\n", my_i);
                            fprintf(stderr, "    current_log_dens[my_i - 1] = %e\n", current_log_dens[my_i - 1]);
                            fprintf(stderr, "    my_coproposal_log_dens = %e\n", my_coproposal_log_dens);
                            for (int j = 0; j <= nx; j++)
                                fprintf(stderr, "    coproposal[%d] = %e\n", j + 1, REAL(coproposal)[j]);
                            fprintf(stderr, "    logh(coproposal) = %e\n", logh(func1, coproposal, rho1));
                        }
#endif /* BLEAT */
                        if (my_coproposal_log_dens == R_NegInf)
                            error("Can't happen: log density -Inf at current state");

                        double my_new_log_dens = logh(func1, proposal, rho1);
                        double my_log_hastings = my_new_log_dens -
                            my_coproposal_log_dens;

                        if (isnan(my_log_hastings) ||
                            (isinf(my_log_hastings) && my_log_hastings > 0))
                            error("Can't happen: log hastings ratio +Inf or NaN\n");

                        int my_accept = 1;
                        double my_unif_hastings = R_NaReal;
                        if (my_log_hastings < 0.0) {
                            my_unif_hastings = unif_rand();
                            my_accept = my_unif_hastings < exp(my_log_hastings);
                        }

                        if (is_debug) {
                            for (int j = 0; j <= nx; j++)
                                REAL(debug_proposal)[iiter + niter * j] =
                                    REAL(proposal)[j];
                            for (int j = 0; j <= nx; j++)
                                REAL(debug_coproposal)[iiter + niter * j] =
                                    REAL(coproposal)[j];
                            REAL(debug_log_hastings)[iiter] = my_log_hastings;
                            REAL(debug_unif_hastings)[iiter] = my_unif_hastings;
                            LOGICAL(debug_acceptd)[iiter] = my_accept;
                            for (int j = 0; j < nx; j++)
                                REAL(debug_norm)[iiter + niter * j] = z[j];
                            REAL(debug_unif_choose)[iiter] = unif_choose;
                            REAL(debug_unif_choose)[iiter + niter] = R_NaReal;
                        }

                        if (my_accept) {
                            for (int j = 0; j < nx; j++)
                                REAL(state)[(my_i - 1) + ncomp * j] =
                                    REAL(proposal)[j + 1];
                            current_log_dens[my_i - 1] = my_new_log_dens;
                            acceptx_numer[my_i - 1]++;
                        }
                        acceptx_denom[my_i - 1]++;

                    } else /* serial */ {

                        int my_i = REAL(state)[0];
                        if (my_i <= 0 || my_i > ncomp)
                            error("Can't happen: my_i out of range");

                        REAL(coproposal)[0] = my_i;
                        for (int j = 0; j < nx; j++)
                            REAL(coproposal)[j + 1] = REAL(state)[j + 1];

                        double z[nx];
                        propose(coproposal, proposal, scale, z);
                        double my_new_log_dens = logh(func1, proposal, rho1);
                        double my_old_log_dens = current_log_dens[my_i - 1];

#ifdef BLEAT
                        if (my_old_log_dens == R_NegInf) {
                            fprintf(stderr, "Oopsie #2!\n");
                        }
#endif /* BLEAT */
                        if (my_old_log_dens == R_NegInf)
                            error("Can't happen: log density -Inf at current state");
#ifdef EXTRA_CHECK
                        if (my_old_log_dens != logh(func1, coproposal, rho1)) {
                            fprintf(stderr, "with-in component update (serial)\n");
                            error("saving logh didn't work right (coproposal)");
                        }
#endif /* EXTRA_CHECK */

                        double my_log_hastings =
                            my_new_log_dens - my_old_log_dens;

                        if (isnan(my_log_hastings) ||
                            (isinf(my_log_hastings) && my_log_hastings > 0)) {
#ifdef WOOF
                                fprintf(stderr, "my_old_log_dens = %e\n", my_old_log_dens);
                                fprintf(stderr, "my_new_log_dens = %e\n", my_old_log_dens);
                                fprintf(stderr, "my_i = %d\n", my_i);
#endif /* WOOF */
                                error("Can't happen: log hastings ratio +Inf or NaN\n");
                        }

                        int my_accept = 1;
                        double my_unif_hastings = R_NaReal;
                        if (my_log_hastings < 0.0) {
                            my_unif_hastings = unif_rand();
                            my_accept = my_unif_hastings < exp(my_log_hastings);
                        }

                        if (is_debug) {
                            for (int j = 0; j <= nx; j++)
                                REAL(debug_proposal)[iiter + niter * j] =
                                    REAL(proposal)[j];
                            REAL(debug_log_hastings)[iiter] = my_log_hastings;
                            REAL(debug_unif_hastings)[iiter] = my_unif_hastings;
                            LOGICAL(debug_acceptd)[iiter] = my_accept;
                            for (int j = 0; j < nx; j++)
                                REAL(debug_norm)[iiter + niter * j] = z[j];
                            REAL(debug_unif_choose)[iiter] = R_NaReal;
                        }

                        if (my_accept) {
                            for (int j = 0; j <= nx; j++)
                                REAL(state)[j] = REAL(proposal)[j];
                            current_log_dens[my_i - 1] = my_new_log_dens;
                            acceptx_numer[my_i - 1]++;
                        }
                        acceptx_denom[my_i - 1]++;

                    }

                } else /* jump/swap update */ {

                    if (is_parallel) {

                        double unif_choose_one = unif_rand();
                        double unif_choose_two = unif_rand();

                        int my_i = trunc(ncomp * unif_choose_one) + 1;
                        if (my_i > ncomp) my_i--;
                        if (my_i <= 0 || my_i > ncomp)
                            error("Can't happen: my_i out of range");

                        REAL(coproposal)[0] = my_i;
                        for (int j = 0; j < nx; j++)
                            REAL(coproposal)[j + 1] =
                                REAL(state)[(my_i - 1) + ncomp * j];

                        int my_i_neighbors = n_neighbors[my_i - 1];

                        int foo = trunc(my_i_neighbors * unif_choose_two) + 1;
                        if (foo > my_i_neighbors) foo--;

                        int my_j = the_neighbors[my_i - 1][foo - 1] + 1;
#ifdef BLEAT
                        fprintf(stderr, "my_i = %d, my_i_neighbors = %d, foo = %d\n", my_i, my_i_neighbors, foo);
                        fprintf(stderr, "(parallel) ncomp = %d, my_j = %d\n", ncomp, my_j);
#endif /* BLEAT */
                        if (my_j <= 0 || my_j > ncomp)
                            error("Can't happen: my_j out of range");

                        REAL(proposal)[0] = my_j;
                        for (int j = 0; j < nx; j++)
                            REAL(proposal)[j + 1] =
                                REAL(state)[(my_j - 1) + ncomp * j];

                        double my_coproposal_log_dens =
                            current_log_dens[my_i - 1];

#ifdef EXTRA_CHECK
                        if (my_coproposal_log_dens != logh(func1, coproposal, rho1)) {
                            fprintf(stderr, "swap component update (parallel)\n");
                            error("saving logh didn't work right (coproposal)");
                        }
#endif /* EXTRA_CHECK */
#ifdef BLEAT
                        if (my_coproposal_log_dens == R_NegInf) {
                            fprintf(stderr, "Oopsie #3!\n");
                            fprintf(stderr, "    my_i = %d\n", my_i);
                            fprintf(stderr, "    current_log_dens[my_i - 1] = %e\n", current_log_dens[my_i - 1]);
                            fprintf(stderr, "    my_coproposal_log_dens = %e\n", my_coproposal_log_dens);
                            for (int j = 0; j <= nx; j++)
                                fprintf(stderr, "    coproposal[%d] = %e\n", j + 1, REAL(coproposal)[j]);
                            fprintf(stderr, "    logh(coproposal) = %e\n", logh(func1, coproposal, rho1));
                        }
#endif /* BLEAT */
                        if (my_coproposal_log_dens == R_NegInf)
                            error("Can't happen: log density -Inf at current state");

                        double my_proposal_log_dens =
                            current_log_dens[my_j - 1];

#ifdef EXTRA_CHECK
                        if (my_proposal_log_dens != logh(func1, proposal, rho1)) {
                            fprintf(stderr, "swap component update (parallel)\n");
                            error("saving logh didn't work right (proposal)");
                        }
#endif /* EXTRA_CHECK */
#ifdef BLEAT
                        if (my_proposal_log_dens == R_NegInf) {
                            fprintf(stderr, "Oopsie #4!\n");
                        }
#endif /* BLEAT */
                        if (my_proposal_log_dens == R_NegInf)
                            error("Can't happen: log density -Inf at current state");

                        if (is_debug) {
                            for (int j = 0; j <= nx; j++)
                                REAL(debug_proposal)[iiter + niter * j] =
                                    REAL(proposal)[j];
                            for (int j = 0; j <= nx; j++)
                                REAL(debug_coproposal)[iiter + niter * j] =
                                    REAL(coproposal)[j];
                        }

                        // proposal and coproposal now saved and logh evaluated
                        // for them, can clobber to evaluate for swap

                        REAL(proposal)[0] = my_i;
                        REAL(coproposal)[0] = my_j;
                        double my_swapped_coproposal_log_dens =
                            logh(func1, coproposal, rho1);
                        double my_swapped_proposal_log_dens =
                            logh(func1, proposal, rho1);
                        double my_log_hastings = my_swapped_proposal_log_dens +
                            my_swapped_coproposal_log_dens -
                            my_proposal_log_dens - my_coproposal_log_dens;

                        if (isnan(my_log_hastings) ||
                            (isinf(my_log_hastings) && my_log_hastings > 0))
                            error("Can't happen: log hastings ratio +Inf or NaN\n");

                        int my_accept = 1;
                        double my_unif_hastings = R_NaReal;
                        if (my_log_hastings < 0.0) {
                            my_unif_hastings = unif_rand();
                            my_accept = my_unif_hastings < exp(my_log_hastings);
                        }

                        if (is_debug) {
                            REAL(debug_log_hastings)[iiter] = my_log_hastings;
                            REAL(debug_unif_hastings)[iiter] = my_unif_hastings;
                            LOGICAL(debug_acceptd)[iiter] = my_accept;
                            for (int j = 0; j < nx; j++)
                                REAL(debug_norm)[iiter + niter * j] = R_NaReal;
                            REAL(debug_unif_choose)[iiter] = unif_choose_one;
                            REAL(debug_unif_choose)[iiter + niter] =
                                unif_choose_two;
                        }

                        if (my_accept) {
                            for (int j = 0; j < nx; j++)
                                REAL(state)[(my_j - 1) + ncomp * j] =
                                    REAL(coproposal)[j + 1];
                            for (int j = 0; j < nx; j++)
                                REAL(state)[(my_i - 1) + ncomp * j] =
                                    REAL(proposal)[j + 1];
                            current_log_dens[my_i - 1] =
                                my_swapped_proposal_log_dens;
                            current_log_dens[my_j - 1] =
                                my_swapped_coproposal_log_dens;
                            accepti_numer[my_i - 1][my_j - 1]++;
                        }
                        accepti_denom[my_i - 1][my_j - 1]++;

                    } else /* serial */ {

                        int my_i = REAL(state)[0];
                        if (my_i <= 0 || my_i > ncomp)
                            error("Can't happen: my_i out of range");

                        int my_i_neighbors = n_neighbors[my_i - 1];

                        double unif_choose = unif_rand();
                        int foo = trunc(my_i_neighbors * unif_choose) + 1;
                        if (foo > my_i_neighbors) foo--;

                        int my_j = the_neighbors[my_i - 1][foo - 1] + 1;
                        int my_j_neighbors = n_neighbors[my_j - 1];

#ifdef BLEAT
                        fprintf(stderr, "(serial) ncomp = %d, my_j = %d, my_i_neighbors = %d, foo = %d\n", ncomp, my_j, my_i_neighbors, foo);
                        fprintf(stderr, "    unif_choose = %f\n", unif_choose);
                        for (int j = 0; j < ncomp; j++)
                            fprintf(stderr, "    LOGICAL(neighbors)[(my_i - 1) + ncomp * %d] = %d\n", j, LOGICAL(neighbors)[(my_i - 1) + ncomp * j]);
#endif /* BLEAT */
                        if (my_j <= 0 || my_j > ncomp)
                            error("Can't happen: my_j out of range");

                        REAL(proposal)[0] = my_j;
                        for (int j = 0; j < nx; j++)
                            REAL(proposal)[j + 1] = REAL(state)[j + 1];

#ifdef WOOF
                        fprintf(stderr, "got to here, about to call logh\n");
                        fprintf(stderr, "    REAL(proposal)[0] = %e\n", REAL(proposal)[0]);
#endif /* WOOF */

                        double my_new_log_dens = logh(func1, proposal, rho1);
                        double my_old_log_dens = current_log_dens[my_i - 1];
#ifdef BLEAT
                        if (my_old_log_dens == R_NegInf) {
                            fprintf(stderr, "Oopsie #5!\n");
                        }
#endif /* BLEAT */
#ifdef EXTRA_CHECK
                        for (int j = 0; j <= nx; j++)
                            REAL(coproposal)[j] = REAL(state)[j];
                        if (my_old_log_dens != logh(func1, coproposal, rho1)) {
                            fprintf(stderr, "swap component update (serial)\n");
                            error("saving logh didn't work right (coproposal)");
                        }
#endif /* EXTRA_CHECK */
                        if (my_old_log_dens == R_NegInf)
                            error("Can't happen: log density -Inf at current state");

                        double my_log_hastings =
                            my_new_log_dens - my_old_log_dens +
                            log(my_i_neighbors) - log(my_j_neighbors);

                        if (isnan(my_log_hastings) ||
                            (isinf(my_log_hastings) && my_log_hastings > 0))
                            error("Can't happen: log hastings ratio +Inf or NaN\n");

                        int my_accept = 1;
                        double my_unif_hastings = R_NaReal;
                        if (my_log_hastings < 0.0) {
                            my_unif_hastings = unif_rand();
                            my_accept = my_unif_hastings < exp(my_log_hastings);
                        }

                        if (is_debug) {
                            for (int j = 0; j <= nx; j++)
                                REAL(debug_proposal)[iiter + niter * j] =
                                    REAL(proposal)[j];
                            REAL(debug_log_hastings)[iiter] = my_log_hastings;
                            REAL(debug_unif_hastings)[iiter] = my_unif_hastings;
                            LOGICAL(debug_acceptd)[iiter] = my_accept;
                            for (int j = 0; j < nx; j++)
                                REAL(debug_norm)[iiter + niter * j] = R_NaReal;
                            REAL(debug_unif_choose)[iiter] = unif_choose;
                        }

#ifdef WOOF_WOOF
                        fprintf(stderr, "unif_choose = %f, ", unif_choose);
                        fprintf(stderr, "REAL(debug_unif_choose)[iiter] = %f\n", REAL(debug_unif_choose)[iiter]);
#endif /* WOOF_WOOF */

                        if (my_accept) {
                            for (int j = 0; j <= nx; j++)
                                REAL(state)[j] = REAL(proposal)[j];
                            current_log_dens[my_j - 1] = my_new_log_dens;
                            accepti_numer[my_i - 1][my_j - 1]++;
                        }
                        accepti_denom[my_i - 1][my_j - 1]++;

                    }
                }

             } /* end of inner loop (one iteration) */

            if (no_outfun) {
                if (is_parallel)
                   for (int i = 0; i < nout; i++)
                       batch_buff[i] += REAL(state)[i];
                else
                   for (int i = 0; i < nout; i++)
                       batch_buff[i] += REAL(state)[i + 1];
            } else /* has outfun */ {
                SEXP fred = outfun(func2, state, rho2);
                if (LENGTH(fred) != nout)
                    error("function outfun returns results of different lengths");
                for (int i = 0; i < nout; i++)
                    batch_buff[i] += REAL(fred)[i];
            }

            if (! is_parallel)
                ibatch_buff[((int) REAL(state)[0]) - 1]++;

        } /* end of middle loop (one batch) */

        if (no_outfun && is_parallel)
            for (int i = 0; i < ncomp; i++)
                for (int j = 0; j < nx; j++)
                    REAL(batch)[kbatch + int_nbatch * (i + ncomp * j)] =
                        batch_buff[i + ncomp * j] / int_blen;
        else
            for (int i = 0; i < nout; i++)
                REAL(batch)[kbatch + int_nbatch * i] =
                    batch_buff[i] / int_blen;

        if (! is_parallel)
            for (int i = 0; i < ncomp; i++)
                REAL(ibatch)[kbatch + int_nbatch * i] =
                    ibatch_buff[i] / int_blen;

    } /* end of outer loop */

    for (int i = 0; i < ncomp; i++)
        REAL(acceptx)[i] = acceptx_numer[i] / acceptx_denom[i];

    for (int i = 0; i < ncomp; i++)
        for (int j = 0; j < ncomp; j++)
            if (LOGICAL(neighbors)[i + ncomp * j])
                REAL(accepti)[i + ncomp * j] = accepti_numer[i][j] / accepti_denom[i][j];
            else
                REAL(accepti)[i + ncomp * j] = R_NaReal;

    PutRNGstate();

    PROTECT(save_final = duplicate(state));
    SET_VECTOR_ELT(result, 4, save_final);
    UNPROTECT(5);

    return result;
}

static void check_valid_scale(SEXP scale, int i, int ncomp, int nx)
{
    if (i > ncomp)
        error("check_valid_scale: i = %d, ncomp = %d, invalid\n", i, ncomp);

    if (! isReal(scale)) {
        if (i >= 0)
            error("component %d of scale not type double", i + 1);
        else
            error("scale not type double");
    }
    if (! isAllFinite(scale)) {
        if (i >= 0)
            error("component %d of scale has non-finite element", i + 1);
        else
            error("scale has non-finite element");
    }
    if (isMatrix(scale)) {
        if (nrows(scale) != nx) {
            if (i >= 0)
                error("component %d of scale matrix with wrong row dim", i + 1);
            else
                error("scale matrix with wrong row dim");
        }
        if (ncols(scale) != nx) {
            if (i >= 0)
                error("component %d of scale matrix with wrong col dim", i + 1);
            else
                error("scale matrix with wrong col dim");
        }
    } else /* scale not matrix */ {
        if (! (LENGTH(scale) == 1 || LENGTH(scale) == nx)) {
            if (i >= 0)
                error("component %d of scale not matrix, scalar, or vector of length k", i + 1);
            else
                error("scale not matrix, scalar, or vector of length k");
        }
    }
}

static double logh(SEXP func, SEXP state, SEXP rho)
{
     SEXP call, result, foo;
     double bar;

     PROTECT(call = lang2(func, state));
     PROTECT(result = eval(call, rho));
     if (! isNumeric(result))
         error("log unnormalized density function returned non-numeric");
     if (LENGTH(result) != 1)
         error("log unnormalized density function returned non-scalar");
     PROTECT(foo = coerceVector(result, REALSXP));
     bar = REAL(foo)[0];
     UNPROTECT(3);
     if (bar == R_PosInf)
         error("log unnormalized density function returned +Inf");
     if (R_IsNaN(bar) || R_IsNA(bar))
         error("log unnormalized density function returned NA or NaN");
     /* Note: -Inf is allowed */
     return bar;
}

static SEXP outfun(SEXP func, SEXP state, SEXP rho)
{
     SEXP call, result, foo;

     PROTECT(call = lang2(func, state));
     PROTECT(result = eval(call, rho));
     if (! isNumeric(result))
         error("outfun returned non-numeric");
     PROTECT(foo = coerceVector(result, REALSXP));
     UNPROTECT(3);
     return foo;
}

static void propose(SEXP coproposal, SEXP proposal, SEXP scale, double *z)
{
    int my_i = REAL(coproposal)[0];
    int nx = LENGTH(coproposal) - 1;

    for (int j = 0; j < nx; j++)
        z[j] = norm_rand();

    if (isNewList(scale))
        scale = VECTOR_ELT(scale, my_i - 1);

    REAL(proposal)[0] = my_i;

    if (LENGTH(scale) == 1) {

        for (int j = 0; j < nx; j++)
            REAL(proposal)[j + 1] = REAL(coproposal)[j + 1] +
                REAL(scale)[0] * z[j];

    } else if (LENGTH(scale) == nx) {

        for (int j = 0; j < nx; j++)
            REAL(proposal)[j + 1] = REAL(coproposal)[j + 1] +
                REAL(scale)[j] * z[j];

    } else /* scale is nx by nx matrix */ {

        for (int j = 0; j < nx; j++)
            REAL(proposal)[j + 1] = REAL(coproposal)[j + 1];

        for (int j = 0, m = 0; j < nx; j++) {
                double u = z[j];
                for (int k = 0; k < nx; k++)
                    REAL(proposal)[k + 1] += REAL(scale)[m++] * u;
            }
    }
}

