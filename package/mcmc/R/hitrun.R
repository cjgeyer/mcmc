
# a1 and b1 define inequality constraints (a1 %*% x <= b1)
# a2 and b2 define inequality constraints (a2 %*% x == b2)
# as in the makeH function in the rcdd package
# otherwise arguments are like metrop in this package
# a1 and b1 cannot be missing (necessarily results in unbounded polyhedron)
# a2 and b2 can be missing

hitrun <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, a1, b1, a2, b2, outfun, debug = FALSE, ...)
UseMethod("hitrun")

hitrun.hitrun <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, a1, b1, a2, b2, outfun, debug = FALSE, ...)
{
    if (missing(nbatch)) nbatch <- obj$nbatch
    if (missing(blen)) blen <- obj$blen
    if (missing(nspac)) nspac <- obj$nspac
    if (missing(a1)) a1 <- obj$a1
    if (missing(b1)) b1 <- obj$b1
    if (missing(a2)) a1 <- obj$a2
    if (missing(b2)) b1 <- obj$b2
    if (missing(debug)) debug <- obj$debug
    assign(".Random.seed", obj$final.seed, .GlobalEnv)
    if (missing(outfun)) {
        if (is.null(obj$outfun)) {
            hitrun.function(obj$lud, obj$final, nbatch, blen,
                nspac, scale, debug = debug, ...)
        } else {
            hitrun.function(obj$lud, obj$final, nbatch, blen,
                nspac, scale, obj$outfun, debug, ...)
        }
    } else {
        hitrun.function(obj$lud, obj$final, nbatch, blen,
            nspac, a1, b1, a2, b2, outfun, debug, ...)
    }
}

hitrun.function <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, a1, b1, a2, b2, outfun, debug = FALSE, ...)
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
    stopifnot(is.numeric(a1))
    stopifnot(is.finite(a1))
    stopifnot(is.matrix(a1))
    stopifnot(is.numeric(b1))
    stopifnot(is.finite(b1))
    stopifnot(is.vector(b1))
    stopifnot(nrow(a1) == length(b1))
    stopifnot(ncol(a1) == length(initial))
    stopifnot(missing(a2) == missing(b2))
    if (! missing(a2)) {
        stopifnot(is.numeric(a2))
        stopifnot(is.finite(a2))
        stopifnot(is.matrix(a2))
        stopifnot(is.numeric(b2))
        stopifnot(is.finite(b2))
        stopifnot(is.vector(b2))
        stopifnot(nrow(a2) == length(b2))
        stopifnot(ncol(a2) == length(initial))
    }
    stopifnot(is.logical(debug))
    stopifnot(length(debug) == 1)

    hrep1 <- makeH(a1, b1, a2, b2)
    hrep1 <- d2q(hrep1)
    hrep2 <- redundant(hrep1)$output
    hrep3 <- hrep2[hrep2[ , 1] == "1", , drop = FALSE]
    hrep4 <- hrep2[hrep2[ , 1] == "0", , drop = FALSE]

    vrep3 <- scdd(hrep3, representation = "H")$output
    is.line <- vrep3[ , 1] == "1" & vrep3[ , 2] == "0"
    is.point <- vrep3[ , 1] == "0" & vrep3[ , 2] == "1"
    if (! all(is.point | is.line))
        stop("unexpected V-representation of affine hull of constraint set")
    if (sum(is.point) != 1)
        stop("unexpected V-representation of affine hull of constraint set")
    if (sum(is.line) != length(initial) - nrow(hrep3))
        stop("unexpected V-representation of affine hull of constraint set")

    foo <- vrep3[ , - c(1, 2)]
    origin <- foo[is.point, ]
    basis <- foo[is.line, ]
    origin <- q2d(origin)
    basis <- q2d(basis)
    basis <- t(basis)

    amat <- qneg(hrep4[ , - c(1, 2)])
    bvec <- hrep4[ , 2]
    bvec <- qmq(bvec, qmatmult(amat, cbind(origin)))
    amat <- qmatmult(amat, basis)

    hrep5 <- cbind("0", bvec, qneg(amat))
    vrep5 <- scdd(hrep5, representation = "H")$output
    if (nrow(vrep5) == 0)
        stop("empty constraint set")
    is.point <- vrep5[ , 1] == "0" & vrep5[ , 2] == "1"
    if (! all(is.point))
        stop("unbounded constraint set")
    v <- vrep5[ , - c(1, 2)]
    rip <- apply(v, 2, qsum)
    rip <- qdq(rip, rep(as.character(nrow(v)), length(rip)))
    rip <- q2d(rip)
    # rip is relative interior point of constraint set (see design doc)

    bvec <- q2d(bvec)
    amat <- q2d(amat)

    func1 <- function(state) {
        userstate <- origin + as.numeric(basis %*% state)
        obj(userstate, ...)
    }
    env1 <- environment(fun = func1)

    if (missing(outfun)) {
        func2 <- NULL
        env2 <- NULL
        outfun <- NULL
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
    out <- .Call("har", func1, rip, nbatch, blen, nspac,
        amat, bvec, func2, debug, env1, env2, PACKAGE = "mcmc")
    )

    if (is.null(outfun)) {
        foo <- out$batch
        foo <- basis %*% foo
        foo <- sweep(foo, 1, origin, "+")
        out$batch <- foo
    }

    out$initial.seed <- saveseed
    out$final.seed <- .Random.seed
    out$time <- out.time
    out$lud <- obj
    out$nbatch <- nbatch
    out$blen <- blen
    out$nspac <- nspac
    out$a1 <- a1
    out$b1 <- b1
    out$a2 <- a2
    out$b2 <- b2
    out$outfun <- outfun
    out$batch <- t(out$batch)
    out$debug <- debug
    if (debug) {
        out$current <- t(out$current)
        out$proposal <- t(out$proposal)
        out$z <- t(out$z)
        out$origin <- origin
        out$basis <- basis
        out$amat <- amat
        out$bvec <- bvec
        out$rip <- rip
    }
    class(out) <- c("mcmc", "hitrun")
    return(out)
}

