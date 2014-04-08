
hitrun <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, amat, bvec, outfun, debug = FALSE, ...)
UseMethod("hitrun")

hitrun.hitrun <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, amat, bvec, outfun, debug = FALSE, ...)
{
    if (missing(nbatch)) nbatch <- obj$nbatch
    if (missing(blen)) blen <- obj$blen
    if (missing(nspac)) nspac <- obj$nspac
    if (missing(amat)) amat <- obj$amat
    if (missing(bvec)) bvec <- obj$bvec
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
            nspac, scale, outfun, debug, ...)
    }
}

hitrun.function <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, amat, bvec, outfun, debug = FALSE, ...)
{
    if (! exists(".Random.seed")) runif(1)
    saveseed <- .Random.seed
    func1 <- function(state) obj(state, ...)
    env1 <- environment(fun = func1)
    if (missing(outfun)) {
        func2 <- NULL
        env2 <- NULL
        outfun <- NULL
    } else if (is.function(outfun)) {
        func2 <- function(state) outfun(state, ...)
        env2 <- environment(fun = func2)
    } else {
        stop("outfun must be function or missing")
    }

    hrep <- makeH(amat, bvec)
    vrep <- q2d(scdd(d2q(hrep)))
    if (any(vrep[ , 2] == 0))
        stop("unbounded constraint set")
    hrep2 <- q2d(redundant(d2q(hrep)))
    hrep3 <- hrep2[hrep2[ , 1] == 1]
    vrep3 <- q2d(scdd(d2q(hrep3, rep = "H")))
    hrep4 <- hrep2[hrep2[ , 1] == 0]

    out.time <- system.time(
    out <- .Call("har", func1, initial, nbatch, blen, nspac,
        amat, bvec, func2, debug, env1, env2, PACKAGE = "mcmc")
    )
    out$initial.seed <- saveseed
    out$final.seed <- .Random.seed
    out$time <- out.time
    out$lud <- obj
    out$nbatch <- nbatch
    out$blen <- blen
    out$nspac <- nspac
    out$scale <- scale
    out$outfun <- outfun
    out$batch <- t(out$batch)
    out$debug <- debug
    if (! is.null(out$current)) out$current <- t(out$current)
    if (! is.null(out$proposal)) out$proposal <- t(out$proposal)
    class(out) <- c("mcmc", "hitrun")
    return(out)
}

