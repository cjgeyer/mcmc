
# a1 and b1 define inequality constraints (a1 %*% x <= b1)
# inequality constraints x >= 0 are not part of a1 and b1
# a2 and b2 define inequality constraints (a2 %*% x == b2)
# equality constraint sum(x) == 1 is not part of a2 and b2
# otherwise arguments are like metrop in this package
# a1, b1, a2, and b2 can be expressed either in ordinary computer
#    numbers or in rational-numbers (character strings)
# first argument is either numeric (Dirichlet parameter vector alpha)
#    or an object of class "condir" that is the result of a previous run

condir <- function(obj, nbatch, blen = 1,
    a1, b1, a2, b2, outfun, mixprob = rep(1 / 3, 3), debug = FALSE, ...)
UseMethod("condir")

# construct to pass to C
#    vector "origin" and matrix "basis" that define the map
#            function(w) origin + basis %*% w
#        that goes from new coordinates (NC) to original coordinates (OC)
#    Amat and bvec, which define H-representation in NC
#        constraints are Amat %*% w <= bvec
#    V3, vertex set, which defines V-representation in NC
#    mu and Gamma, mean vector and Cholesky factor of variance matrix
#        of asymptotic normal distribution
#    alpha, Dirichlet parameter vector
#    alpha_cons, numeric vector the length of bvec (number of non-redundant
#        inequality constraints).  For each constraint, at most one coordinate
#        in OC is equal to zero for all x satisfying the constraint.
#        If x_i = 0 for all x satisfying the j-th constraint, then we
#        set alpha_cons[j] = alpha[i].  Otherwise, we set alpha_cons[j] = 1.

condir.default <- function(obj, nbatch, blen = 1, nspac = 1,
    a1, b1, a2, b2, outfun, mixprob = rep(1 / 3, 3), debug = FALSE, ...)
{
    if (! exists(".Random.seed")) runif(1)
    saveseed <- .Random.seed

    stopifnot(is.numeric(obj))
    stopifnot(is.finite(obj))
    stopifnot(is.vector(obj))
    stopifnot(obj > 0)
    d <- length(obj)

    stopifnot(is.numeric(nbatch))
    stopifnot(length(nbatch) == 1)
    stopifnot(round(nbatch) == nbatch)
    stopifnot(nbatch > 0)
    stopifnot(is.numeric(blen))
    stopifnot(length(blen) == 1)
    stopifnot(round(blen) == blen)
    stopifnot(blen > 0)
    stopifnot(is.numeric(nspac))
    stopifnot(length(nspac) == 1)
    stopifnot(round(nspac) == nspac)
    stopifnot(nspac > 0)
    nbatch <- as.integer(round(nbatch))
    blen <- as.integer(round(blen))
    nspac <- as.integer(round(nspac))

    stopifnot(missing(a1) == missing(b1))
    if (! missing(a1)) {
        stopifnot(is.numeric(a1) | is.character(a1))
        stopifnot(is.matrix(a1))
        stopifnot(ncol(a1) == length(obj))
        if (is.numeric(a1)) {
            stopifnot(is.finite(a1))
            a1 <- d2q(a1)
        } else {
            a1 <- q2q(a1)
            if (inherits(a1, "try-error"))
                stop("a1 of type character but not valid rational numbers")
        }
        stopifnot(is.numeric(b1) | is.character(b1))
        stopifnot(is.vector(b1))
        stopifnot(nrow(a1) == length(b1))
        if (is.numeric(b1)) {
            stopifnot(is.finite(b1))
            b1 <- d2q(b1)
        } else {
            b1 <- q2q(b1)
            if (inherits(b1, "try-error"))
                stop("b1 of type character but not valid rational numbers")
        }
        a1 <- rbind(- diag(d), a1)
        b1 <- c(rep(0, d), b1)
    } else {
        a1 <- (- diag(d))
        b1 <- rep(0, d)
    }
    stopifnot(missing(a2) == missing(b2))
    if (! missing(a2)) {
        stopifnot(is.numeric(a2) | is.character(a2))
        stopifnot(is.matrix(a2))
        stopifnot(ncol(a2) == length(obj))
        if (is.numeric(a2)) {
            stopifnot(is.finite(a2))
            a2 <- d2q(a2)
        } else {
            a2 <- q2q(a2)
            if (inherits(a2, "try-error"))
                stop("a2 of type character but not valid rational numbers")
        }
        stopifnot(is.numeric(b2) | is.character(b2))
        stopifnot(is.vector(b2))
        stopifnot(nrow(a2) == length(b2))
        if (is.numeric(b2)) {
            stopifnot(is.finite(b2))
            b2 <- d2q(b2)
        } else {
            b2 <- q2q(b2)
            if (inherits(b2, "try-error"))
                stop("b2 of type character but not valid rational numbers")
        }
        a2 <- rbind(1, a2)
        b2 <- c(1, b2)
    } else {
        a2 <- matrix(1, nrow = 1, ncol = d)
        b2 <- 1
    }

    # outfun gets checked later on

    stopifnot(is.numeric(mixprob))
    stopifnot(is.finite(mixprob))
    stopifnot(is.vector(mixprob))
    stopifnot(length(mixprob) == 3)
    stopifnot(mixprob > 0)
    stopifnot(all.equal(sum(mixprob), 1))
    mixprob <- as.double(mixprob)

    stopifnot(is.logical(debug))
    stopifnot(length(debug) == 1)

    hrep1 <- rbind(cbind(0, b1, - a1), cbind(1, b2, - a2))
    dimnames(hrep1) <- NULL
    vrep1 <- scdd(hrep1, representation = "H")$output
    if (nrow(vrep1) == 0)
        stop("empty constraint set")
    if (nrow(vrep1) == 1)
        stop("constraint set is single point")

    # find point in the relative interior of the constraint set
    v1 <- vrep1[ , - c(1, 2)]
    origin <- apply(v1, 2, function(x) qdq(qsum(x), d2q(length(x))))

    # non-redundant H-representation
    rout <- redundant(hrep1)
    hrep2 <- rout$output
    # keep track of where nonnegativity constraints are now
    nonneg.position <- rout$new.position[1:d]

    # affine hull of constraint set
    is.equality <- hrep2[ , 1] == "1"
    hrep3 <- hrep2[is.equality, , drop = FALSE]
    vrep3 <- scdd(hrep3, representation = "H")$output
    # make invertible linear transformation R^p to aff C
    is.point <- vrep3[ , 1] == "0" & vrep3[ , 2] == "1"
    is.line <- vrep3[ , 1] == "1" & vrep3[ , 2] == "0"
    if(! all(is.point | is.line))
        stop("unexpected V-representation of affine hull of constraint set")
    if (sum(is.point) != 1)
        stop("unexpected V-representation of affine hull of constraint set")
    if (sum(is.line) != d - nrow(hrep3))
        stop("unexpected V-representation of affine hull of constraint set")
    basis <- vrep3[is.line, , drop = FALSE]
    basis <- basis[ , - c(1, 2), drop = FALSE]
    basis <- t(basis)
    # invertible linear transformation is
    #     function(w) origin + basis %*% w
    # maps NC to OC
    p <- ncol(basis)

    # constraints in NC
    hrep4 <- hrep2[! is.equality, , drop = FALSE]
    b4 <- hrep4[ , 2]
    a4 <- qneg(hrep4[ , - c(1, 2), drop = FALSE])
    amat <- qmatmult(a4, basis)
    bvec <- qmq(b4, qmatmult(a4, cbind(origin)))

    # vertices in NC
    hrep5 <- makeH(amat, bvec)
    vrep5 <- scdd(hrep5)$output
    vert <- vrep5[ , - c(1, 2)]

    # convert everything from rational to double
    origin <- q2d(origin)
    basis <- q2d(basis)
    bvec <- q2d(bvec)
    amat <- q2d(amat)
    vert <- q2d(vert)

    # asymptotic normal distribution
    alpha <- as.double(obj)
    lambda <- alpha / sum(alpha)
    LambdaInv <- diag(1 / lambda)
    SigmaInv <- t(basis) %*% LambdaInv %*% basis
    cholSigmaInv <- chol(SigmaInv)
    cholSigma <- backsolve(cholSigmaInv, diag(p))
    Sigma <- cholSigma %*% t(cholSigma)
    mu <- Sigma %*% t(basis) %*% cbind(lambda - origin)

    # now where are the nonnegativity constraints?
    # see condir-dev.R in the devel directory for checks that this works
    nonneg.position.varb <- seq(along = nonneg.position)[nonneg.position != 0]
    nonneg.position.cons <- nonneg.position[nonneg.position != 0]
    idx.ineq <- (1:nrow(hrep2))[! is.equality]
    i5 <- match(nonneg.position.cons, idx.ineq)
    alpha.cons <- rep(1, nrow(amat))
    alpha.cons[i5] <- alpha[nonneg.position.varb]

    # outfun must map from NC to OC
    if (missing(outfun)) {
        func2 <- NULL
        env2 <- NULL
    } else if (is.function(outfun)) {
        func2 <- function(state) {
            userstate <- origin + as.numeric(basis %*% state)
            outfun(userstate, ...)
        }
        env2 <- environment(fun = func2)
    } else {
        stop("outfun must be function or missing")
    }

    out.time <- system.time(
    out <- .Call("condir", alpha, rep(0, p), nbatch, blen, nspac,
        origin, basis, amat, bvec, vert, mu, CholSigma, alpha.cons,
        mixprob, func2, env2, debug, PACKAGE = "mcmc")
    )

    ########## REVISED DOWN TO HERE ##########

    if (is.null(func2)) {
        foo <- out$batch
        foo <- foo %*% t(basis)
        foo <- sweep(foo, 2, origin, "+")
        out$batch <- foo
    }

    out$a1 <- a1
    out$b1 <- b1
    out$a2 <- a2
    out$b2 <- b2
    out$origin <- origin
    out$basis <- basis
    class(out) <- c("mcmc", "condir")
    return(out)
}

condir.condir <- function(obj, nbatch, blen = 1,
    a1, b1, a2, b2, outfun, mixprob = rep(1 / 3, 3), debug = FALSE, ...)
{
    if (missing(nbatch)) nbatch <- obj$nbatch
    if (missing(blen)) blen <- obj$blen
    if (missing(nspac)) nspac <- obj$nspac
    if (missing(debug)) debug <- obj$debug
    if (! missing(a1)) warning("a1 not changable, ignored")
    if (! missing(b1)) warning("b1 not changable, ignored")
    if (! missing(a2)) warning("a2 not changable, ignored")
    if (! missing(b2)) warning("b2 not changable, ignored")

    assign(".Random.seed", obj$final.seed, .GlobalEnv)
    basis <- obj$basis
    origin <- obj$origin
    amat <- obj$amat
    bvec <- obj$bvec
    if (missing(outfun)) {
        outfunNC <- obj$outfun
    } else {
        outfunNC <- function(state) {
            userstate <- origin + as.numeric(basis %*% state)
            outfun(userstate, ...)
        }
    }

    out <- condirHelperNC(obj$ludfun, obj$final, nbatch, blen, nspac,
        amat, bvec, outfunNC, debug)

    if (is.null(outfunNC)) {
        foo <- out$batch
        foo <- foo %*% t(basis)
        foo <- sweep(foo, 2, origin, "+")
        out$batch <- foo
    }

    out$a1 <- obj$a1
    out$b1 <- obj$b1
    out$a2 <- obj$a2
    out$b2 <- obj$b2
    out$origin <- origin
    out$basis <- basis
    class(out) <- c("mcmc", "condir")
    return(out)
}

# only knows about new coordinates (NC)
# ludfun is log unnormalized density with argument in NC
# initial is vector in NC
# outfun is output function with argument in NC
# constraint set (in NC) is x such that amat %*% x <= bvec
condirHelperNC <- function(ludfun, initial, nbatch, blen, nspac,
    amat, bvec, outfun, debug)
{
    if (! exists(".Random.seed")) runif(1)
    saveseed <- .Random.seed

    stopifnot(is.numeric(initial))
    stopifnot(is.finite(initial))
    stopifnot(is.numeric(nbatch))
    stopifnot(length(nbatch) == 1)
    stopifnot(nbatch == round(nbatch))
    stopifnot(nbatch >= 1)
    stopifnot(is.numeric(blen))
    stopifnot(length(blen) == 1)
    stopifnot(blen == round(blen))
    stopifnot(blen >= 1)
    stopifnot(is.numeric(nspac))
    stopifnot(length(nspac) == 1)
    stopifnot(nspac == round(nspac))
    stopifnot(nspac >= 1)
    stopifnot(is.numeric(amat))
    stopifnot(is.finite(amat))
    stopifnot(is.matrix(amat))
    stopifnot(is.numeric(bvec))
    stopifnot(is.finite(bvec))
    stopifnot(is.vector(bvec))
    stopifnot(nrow(amat) == length(bvec))
    stopifnot(ncol(amat) == length(initial))
    stopifnot(is.logical(debug))
    stopifnot(length(debug) == 1)

    stopifnot(is.function(ludfun))
    env1 <- environment(fun = ludfun)

    if (is.null(outfun)) {
        env2 <- NULL
    } else if (is.function(outfun)) {
        env2 <- environment(fun = outfun)
    } else {
        stop("outfun must be function or NULL")
    }

    out.time <- system.time(
    out <- .Call("condir", ludfun, initial, nbatch, blen, nspac,
        amat, bvec, outfun, debug, env1, env2, PACKAGE = "mcmc")
    )

    out$initial.seed <- saveseed
    out$final.seed <- .Random.seed
    out$time <- out.time
    out$ludfun <- ludfun
    out$nbatch <- nbatch
    out$blen <- blen
    out$nspac <- nspac
    out$amat <- amat
    out$bvec <- bvec
    out$outfun <- outfun
    out$batch <- t(out$batch)
    out$debug <- debug
    if (debug) {
        out$current <- t(out$current)
        out$proposal <- t(out$proposal)
        out$z <- t(out$z)
    }
    class(out) <- c("mcmc", "condir")
    return(out)
}

