
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
#    V3, vertex set, which defines V-representation in NC
#    mu and Gamma, mean vector and Cholesky factor of variance matrix
#        of asymptotic normal distribution
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

    stopifnot(is.logical(debug))
    stopifnot(length(debug) == 1)

    hrep1 <- makeH(a1, b1, a2, b2)
    hrep2 <- redundant(hrep1)$output
    # non-redundant equality constraints
    hrep3 <- hrep2[hrep2[ , 1] == "1", , drop = FALSE]
    # non-redundant inequality constraints
    hrep4 <- hrep2[hrep2[ , 1] == "0", , drop = FALSE]

    # V-representation for affine hull of constraint set
    vrep3 <- scdd(hrep3, representation = "H")$output
    is.line <- vrep3[ , 1] == "1" & vrep3[ , 2] == "0"
    is.point <- vrep3[ , 1] == "0" & vrep3[ , 2] == "1"
    if (! all(is.point | is.line))
        stop("unexpected V-representation of affine hull of constraint set")
    if (sum(is.point) != 1)
        stop("unexpected V-representation of affine hull of constraint set")
    if (sum(is.line) != d - nrow(hrep3))
        stop("unexpected V-representation of affine hull of constraint set")

    foo <- vrep3[ , - c(1, 2), drop = FALSE]
    origin <- foo[is.point, ]
    basis <- foo[is.line, , drop = FALSE]
    basis <- t(basis)
    # function(w) origin + basis %*% w
    # maps NC to OC
    # note: we do not agree with design doc (dirichlet.tex)
    # which has origin a relative interior point of constraint set in OC
    # see rip below for relative interior point in NC

    amat <- qneg(hrep4[ , - c(1, 2), drop = FALSE])
    bvec <- hrep4[ , 2]
    bvec <- qmq(bvec, qmatmult(amat, cbind(origin)))
    amat <- qmatmult(amat, basis)
    # amat %*% w <= bvec
    # is constraints expressed in NC

    hrep5 <- cbind("0", bvec, qneg(amat))
    vrep5 <- scdd(hrep5, representation = "H")$output
    if (nrow(vrep5) == 0)
        stop("empty constraint set")
    is.point <- vrep5[ , 1] == "0" & vrep5[ , 2] == "1"
    if (! all(is.point))
        stop("Can't happen: unbounded constraint set")
    if (nrow(vrep5) == 1)
        stop("one-point constraint set")
    v <- vrep5[ , - c(1, 2), drop = FALSE]
    rip <- apply(v, 2, qsum)
    rip <- qdq(rip, rep(as.character(nrow(v)), length(rip)))
    # rip is relative interior point of constraint set in NC

    origin <- q2d(origin)
    basis <- q2d(basis)
    bvec <- q2d(bvec)
    amat <- q2d(amat)
    rip <- q2d(rip)
    v <- q2d(v)

    if (missing(outfun)) {
        func2 <- NULL
    } else if (is.function(outfun)) {
        func2 <- function(state) {
            userstate <- origin + as.numeric(basis %*% state)
            outfun(userstate, ...)
        }
    } else {
        stop("outfun must be function or missing")
    }

    alpha <- as.double(obj)
    lambda <- alpha / sum(alpha)
    LambdaInv <- diag(1 / lambda)
    SigmaInv <- t(basis) %*% LambdaInv %*% basis
    p <- nrow(SigmaInv)
    cholSigmaInv <- chol(SigmaInv)
    cholSigma <- backsolve(cholSigmaInv, diag(p))
    Sigma <- cholSigma %*% t(cholSigma)
    mu <- Sigma %*% t(basis) %*% cbind(lambda - origin)

    ########## REVISED DOWN TO HERE ##########

    out <- condirHelperNC(obj, nbatch, blen, nspac,
        amat, bvec, v, rip, func2, debug)

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

