
/*
*
* mcmc and MCMC package for R
* Copyright (c) 2009 Charles J. Geyer
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

static void propose(SEXP coproposal, SEXP proposal, SEXP scale);

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

    if (! isNumeric(scale))
        error("argument \"scale\" must be numeric");

    int no_outfun = isNull(func2);
    if (! no_outfun) {
        if (! isFunction(func2))
            error("argument \"outfun\" must be function");
        if (! isEnvironment(rho2))
            error("argument \"rho2\" must be environment");
    }

    double current_log_dens;
    int current_i;
    if (is_parallel) {
        SEXP fred;
        PROTECT(fred = allocVector(REALSXP, nx + 1));
        for (int i = 0; i < ncomp; i++) {
            current_i = i;
            REAL(fred)[0] = i;
            for (int j = 0; j < nx; j++)
                REAL(fred)[j + 1] = REAL(initial)[i + ncomp * j];
            current_log_dens = logh(func1, fred, rho1);
            if (current_log_dens == R_NegInf)
                error("log unnormalized density -Inf at initial state");
        }
        UNPROTECT(1);
    } else /* serial */ {
        current_i = REAL(initial)[0];
        current_log_dens = logh(func1, initial, rho1);
        if (current_log_dens == R_NegInf)
            error("log unnormalized density -Inf at initial state");
    }

    /* at this point, all arguments have been checked for validity */

    // current_i and current_log_dens save info that cuts in half
    // the number of invocations of log unnormalized density for serial case

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
    //
    // debug output (if parallel)
    //
    //     which         vector of length niter (TRUE if within-component)
    //     unif_which    vector of length niter (uniform for deciding which)
    //     state         niter x ncomp x (nx + 1) array (state before)
    //     coproposal    niter x (nx + 1) matrix
    //     proposal      niter x (nx + 1) matrix
    //     log_hastings  niter vector
    //     unif_hastings niter vector (uniform for deciding acceptance)
    //     acceptd       niter vector (TRUE if accept)
    //
    //     for within-component move coproposal and proposal have natural
    //         meaning
    //     for swap move we have 2 coproposals (i, x_i) and (j, x_j)
    //         and 2 proposals (i, x_j) and (j, x_i) but since the information
    //         here is quite redundant we just store (i, x_i) in "coproposal"
    //         and (j, x_j) in "proposal" -- the checker can figure it out

    int len_result_regular = is_parallel ? 5 : 6;
    int len_result_debug = is_debug ? 8 : 7;
    int len_result = len_result_regular + len_result_debug;

    SEXP result, resultnames, acceptx, accepti, batch, ibatch,
        save_initial, save_final, debug_which, debug_unif_which,
        debug_state, debug_coproposal, debug_proposal,
        debug_log_hastings, debug_unif_hastings, debug_acceptd;

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
            mkChar("unifwhich"));
        UNPROTECT(1);

        if (is_parallel)
            PROTECT(debug_state = alloc3DArray(REALSXP, niter, ncomp, nx + 1));
        else
            PROTECT(debug_state = allocMatrix(REALSXP, niter, nx + 1));
        SET_VECTOR_ELT(result, len_result_regular + 2, debug_state);
        SET_STRING_ELT(resultnames, len_result_regular + 2, mkChar("state"));
        UNPROTECT(1);

        PROTECT(debug_log_hastings = allocVector(REALSXP, niter));
        SET_VECTOR_ELT(result, len_result_regular + 3, debug_log_hastings);
        SET_STRING_ELT(resultnames, len_result_regular + 3,
            mkChar("loghastings"));
        UNPROTECT(1);

        PROTECT(debug_unif_hastings = allocVector(REALSXP, niter));
        SET_VECTOR_ELT(result, len_result_regular + 4, debug_unif_hastings);
        SET_STRING_ELT(resultnames, len_result_regular + 4,
            mkChar("unifhastings"));
        UNPROTECT(1);

        PROTECT(debug_proposal = allocMatrix(REALSXP, niter, nx + 1));
        SET_VECTOR_ELT(result, len_result_regular + 5, debug_proposal);
        SET_STRING_ELT(resultnames, len_result_regular + 5,
            mkChar("proposal"));
        UNPROTECT(1);

        if (is_parallel) {
            PROTECT(debug_coproposal = allocMatrix(REALSXP, niter, nx + 1));
            SET_VECTOR_ELT(result, len_result_regular + 5, debug_coproposal);
            SET_STRING_ELT(resultnames, len_result_regular + 5,
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

    // need buffers for batch means

    double *batch_buff = (double *) R_alloc(nout, sizeof(double));
    double ibatch_buff[ncomp];

    for (int kbatch = 0, iiter = 0; kbatch < int_nbatch; kbatch++) {

        for (int i = 0; i < nout; i++)
            batch_buff[i] = 0.0;
        for (int i = 0; i < ncomp; i++)
            ibatch_buff[i] = 0.0;

        for (int jbatch = 0; jbatch < int_blen; jbatch++) {

            for (int ispac = 0; ispac < int_nspac; ispac++, iiter++) {

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

                        int my_i = trunc(ncomp * unif_rand());
                        if (my_i == ncomp) my_i--;
                        if (my_i < 0 || my_i >= ncomp)
                            error("Can't happen: my_i out of range");

                        REAL(coproposal)[0] = my_i;
                        for (int j = 0; j < nx; j++)
                            REAL(coproposal)[j + 1] =
                                REAL(state)[my_i + ncomp * j];

                        propose(coproposal, proposal, scale);

                        double my_coproposal_log_dens = current_i == my_i ?
                            current_log_dens : logh(func1, coproposal, rho1);

                        if (my_coproposal_log_dens == R_NegInf)
                            error("Can't happen: log density -Inf at current state");

                        double my_new_log_dens = logh(func1, proposal, rho1);
                        double my_log_hastings = my_new_log_dens -
                            my_coproposal_log_dens;

                        int my_accept = 1;
                        double my_unif_hastings = R_NaReal;
                        if (my_log_hastings < 0.0) {
                            my_unif_hastings = unif_rand();
                            my_accept = exp(my_log_hastings) < my_unif_hastings;
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
                        }

                        if (my_accept) {
                            for (int j = 0; j <= nx; j++)
                                REAL(state)[my_i + ncomp * j] =
                                    REAL(proposal)[j + 1];
                            current_i = REAL(proposal)[0];
                            current_log_dens = my_new_log_dens;
                            acceptx_numer[my_i]++;
                        }
                        acceptx_denom[my_i]++;

                    } else /* serial */ {

                        int my_i = REAL(state)[0];
                        if (my_i < 0 || my_i >= ncomp)
                            error("Can't happen: my_i out of range");

                        REAL(coproposal)[0] = my_i;
                        for (int j = 0; j < nx; j++)
                            REAL(coproposal)[j + 1] = REAL(state)[j + 1];

                        propose(coproposal, proposal, scale);

                        if (current_i != my_i)
                            error("Can't happen: current_i != my_i");
                        if (current_log_dens == R_NegInf)
                            error("Can't happen: log density -Inf at current state");

                        double my_new_log_dens = logh(func1, proposal, rho1);
                        double my_log_hastings = my_new_log_dens -
                            current_log_dens;

                        int my_accept = 1;
                        double my_unif_hastings = R_NaReal;
                        if (my_log_hastings < 0.0) {
                            my_unif_hastings = unif_rand();
                            my_accept = exp(my_log_hastings) < my_unif_hastings;
                        }

                        if (is_debug) {
                            for (int j = 0; j <= nx; j++)
                                REAL(debug_proposal)[iiter + niter * j] =
                                    REAL(proposal)[j];
                            REAL(debug_log_hastings)[iiter] = my_log_hastings;
                            REAL(debug_unif_hastings)[iiter] = my_unif_hastings;
                            LOGICAL(debug_acceptd)[iiter] = my_accept;
                        }

                        if (my_accept) {
                            for (int j = 0; j <= nx; j++)
                                REAL(state)[j] = REAL(proposal)[j];
                            current_i = REAL(state)[0];
                            current_log_dens = my_new_log_dens;
                            acceptx_numer[my_i]++;
                        }
                        acceptx_denom[my_i]++;

                    }

                } else /* jump/swap update */ {

                    if (is_parallel) {

                        int my_i = trunc(ncomp * unif_rand());
                        if (my_i == ncomp) my_i--;
                        if (my_i < 0 || my_i >= ncomp)
                            error("Can't happen: my_i out of range");

                        REAL(coproposal)[0] = my_i;
                        for (int j = 0; j < nx; j++)
                            REAL(coproposal)[j + 1] =
                                REAL(state)[my_i + ncomp * j];

                        int my_i_neighbors = 0;
                        for (int j = 0; j < ncomp; j++)
                            my_i_neighbors +=
                                LOGICAL(neighbors)[my_i + ncomp * j];

                        int foo = trunc(my_i_neighbors * unif_rand());
                        if (foo == my_i_neighbors) foo--;

                        int my_j;
                        int bar = 0;
                        for (my_j = 0; my_j < ncomp; my_j++) {
                            bar += LOGICAL(neighbors)[my_i + ncomp * my_j];
                            if (bar == foo) break;
                        }
                        if (my_j < 0 || my_j >= ncomp)
                            error("Can't happen: my_j out of range");

                        REAL(proposal)[0] = my_j;
                        for (int j = 0; j < nx; j++)
                            REAL(proposal)[j + 1] =
                                REAL(state)[my_j + ncomp * j];

                        double my_coproposal_log_dens = current_i == my_i ?
                            current_log_dens : logh(func1, coproposal, rho1);

                        if (my_coproposal_log_dens == R_NegInf)
                            error("Can't happen: log density -Inf at current state");

                        double my_proposal_log_dens = current_i == my_j ?
                            current_log_dens : logh(func1, proposal, rho1);

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
                            my_swapped_proposal_log_dens - my_proposal_log_dens
                            - my_coproposal_log_dens;

                        int my_accept = 1;
                        double my_unif_hastings = R_NaReal;
                        if (my_log_hastings < 0.0) {
                            my_unif_hastings = unif_rand();
                            my_accept = exp(my_log_hastings) < my_unif_hastings;
                        }

                        if (is_debug) {
                            REAL(debug_log_hastings)[iiter] = my_log_hastings;
                            REAL(debug_unif_hastings)[iiter] = my_unif_hastings;
                            LOGICAL(debug_acceptd)[iiter] = my_accept;
                        }

                        if (my_accept) {
                            for (int j = 0; j <= nx; j++)
                                REAL(state)[my_i + ncomp * j] =
                                    REAL(coproposal)[j + 1];
                            for (int j = 0; j < nx; j++)
                                REAL(state)[my_j + ncomp * j] =
                                    REAL(proposal)[j + 1];
                            current_i = my_i;
                            current_log_dens = my_swapped_proposal_log_dens;
                            accepti_numer[my_i][my_j]++;
                        }
                        accepti_denom[my_i][my_j]++;

                    } else /* serial */ {

                        int my_i = REAL(state)[0];
                        if (my_i < 0 || my_i >= ncomp)
                            error("Can't happen: my_i out of range");

                        int my_i_neighbors = 0;
                        for (int j = 0; j < ncomp; j++)
                            my_i_neighbors +=
                                LOGICAL(neighbors)[my_i + ncomp * j];

                        int foo = trunc(my_i_neighbors * unif_rand());
                        if (foo == my_i_neighbors) foo--;

                        int my_j;
                        int bar = 0;
                        for (my_j = 0; my_j < ncomp; my_j++) {
                            bar += LOGICAL(neighbors)[my_i + ncomp * my_j];
                            if (bar == foo) break;
                        }
                        if (my_j < 0 || my_j >= ncomp)
                            error("Can't happen: my_j out of range");

                        REAL(proposal)[0] = my_j;
                        for (int j = 0; j < nx; j++)
                            REAL(proposal)[j + 1] = REAL(state)[j + 1];

                        if (current_i != my_i)
                            error("Can't happen: current_i != my_i");
                        if (current_log_dens == R_NegInf)
                            error("Can't happen: log density -Inf at current state");

                        double my_new_log_dens = logh(func1, proposal, rho1);
                        double my_log_hastings = my_new_log_dens -
                            current_log_dens;

                        int my_accept = 1;
                        double my_unif_hastings = R_NaReal;
                        if (my_log_hastings < 0.0) {
                            my_unif_hastings = unif_rand();
                            my_accept = exp(my_log_hastings) < my_unif_hastings;
                        }

                        if (is_debug) {
                            for (int j = 0; j <= nx; j++)
                                REAL(debug_proposal)[iiter + niter * j] =
                                    REAL(proposal)[j];
                            REAL(debug_log_hastings)[iiter] = my_log_hastings;
                            REAL(debug_unif_hastings)[iiter] = my_unif_hastings;
                            LOGICAL(debug_acceptd)[iiter] = my_accept;
                        }

                        if (my_accept) {
                            for (int j = 0; j <= nx; j++)
                                REAL(state)[j] = REAL(proposal)[j];
                            current_i = REAL(state)[0];
                            current_log_dens = my_new_log_dens;
                            accepti_numer[my_i][my_j]++;
                        }
                        accepti_denom[my_i][my_j]++;

                    }
                }

             } /* end of inner loop (one iteration) */

            if (no_outfun) {
                if (is_parallel)
                   for (int i = 0; i < nout; i++)
                       batch_buff[i] = REAL(state)[i];
                else
                   for (int i = 0; i < nout; i++)
                       batch_buff[i] = REAL(state)[i + 1];
            } else /* has outfun */ {
                SEXP fred = outfun(func2, state, rho2);
                if (LENGTH(fred) != nout)
                    error("function outfun returns results of different lengths");
                for (int i = 0; i < nout; i++)
                    batch_buff[i] = REAL(fred)[i];
            }

            if (! is_parallel)
                ibatch_buff[(int) REAL(state)[0]]++;

        } /* end of middle loop (one batch) */

        if (no_outfun && is_parallel)
            for (int i = 0; i < ncomp; i++)
                for (int j = 0; j < nx; j++)
                    REAL(batch)[kbatch + int_nbatch * (i + ncomp * j)] =
                        batch_buff[i + ncomp * j] / int_nbatch;
        else
            for (int i = 0; i < nout; i++)
                REAL(batch)[kbatch + int_nbatch * i] =
                    batch_buff[i] / int_nbatch;

        if (! is_parallel)
            for (int i = 0; i < ncomp; i++)
                REAL(ibatch)[kbatch + int_nbatch * i] =
                    ibatch_buff[i] / int_nbatch;

    } /* end of outer loop */

    for (int i = 0; i < ncomp; i++)
        REAL(acceptx)[i] = acceptx_numer[i] / acceptx_denom[i];

    for (int i = 0; i < ncomp; i++)
        for (int j = 0; j < ncomp; j++)
            if (LOGICAL(neighbors)[i + ncomp * j])
                REAL(accepti)[i] = accepti_numer[i][j] / accepti_denom[i][j];
            else
                REAL(accepti)[i] = R_NaReal;

    PutRNGstate();

    PROTECT(save_final = duplicate(state));
    SET_VECTOR_ELT(result, 4, save_final);
    UNPROTECT(5);

    return result;
}

static void check_valid_scale(SEXP scale, int i, int ncomp, int nx)
{
    if (! isReal(scale))
        if (i >= 0)
            error("component %d of scale not type double", i + 1);
        else
            error("scale not type double");
    if (! isAllFinite(scale))
        if (i >= 0)
            error("component %d of scale has non-finite element", i + 1);
        else
            error("scale has non-finite element");
    if (isMatrix(scale)) {
        if (nrows(scale) != nx)
            if (i >= 0)
                error("component %d of scale matrix with wrong row dim", i + 1);
            else
                error("scale matrix with wrong row dim");
        if (ncols(scale) != nx)
            if (i >= 0)
                error("component %d of scale matrix with wrong col dim", i + 1);
            else
                error("scale matrix with wrong col dim");
    } else /* scale not matrix */ {
        if (! (LENGTH(scale) == 1 && LENGTH(scale) == nx))
            if (i >= 0)
                error("component %d of scale not matrix, scalar, or vector of length k", i + 1);
            else
                error("scale not matrix, scalar, or vector of length k", i + 1);
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

static void propose(SEXP coproposal, SEXP proposal, SEXP scale)
{
    int i = REAL(coproposal)[0];
    int p = LENGTH(coproposal) - 1;

    if (isNewList(scale))
        scale = VECTOR_ELT(scale, i - 1);

    REAL(proposal)[0] = i;
    if (LENGTH(scale) == 1) {

        for (int j = 0; j < p; j++)
            REAL(proposal)[j + 1] = REAL(coproposal)[j + 1] +
                REAL(scale)[0] * norm_rand();

    } else if (LENGTH(scale) == p) {

        for (int j = 0; j < p; j++)
            REAL(proposal)[j + 1] = REAL(coproposal)[j + 1] +
                REAL(scale)[j] * norm_rand();

    } else /* scale is p by p matrix */ {

        for (int j = 0; j < p; j++)
            REAL(proposal)[j + 1] = REAL(coproposal)[j + 1];

        for (int j = 0, m = 0; j < p; j++) {
                double u = norm_rand();
                for (int k = 0; k < p; k++)
                    REAL(proposal)[k + 1] += REAL(scale)[m++] * u;
            }
    }
}

