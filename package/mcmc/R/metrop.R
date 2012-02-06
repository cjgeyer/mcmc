
metrop <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, morph, debug = FALSE, ...)
UseMethod("metrop")

metrop.metropolis <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, ...)
{
    if (missing(nbatch)) nbatch <- obj$nbatch
    if (missing(blen)) blen <- obj$blen
    if (missing(nspac)) nspac <- obj$nspac
    if (missing(scale)) scale <- obj$scale
    if (missing(debug)) debug <- obj$debug
    if (missing(outfun)) outfun <- obj$outfun

    assign(".Random.seed", obj$final.seed, .GlobalEnv)
    out <- metrop.function(obj$lud, obj$final, nbatch, blen,
                           nspac, scale, outfun, debug, ...)
    return(out)
}

metrop.function <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, morph,
    debug = FALSE, ...)
{
    if (! exists(".Random.seed")) runif(1)
    saveseed <- .Random.seed

    if (missing(morph)) morph <- morph.identity()
    if (missing(outfun)) outfun <- NULL
    
    func1 <- morph$lud(obj, ...)
    env1 <- environment(fun = func1)
    func2 <- morph$outfun(outfun, ...)
    env2 <- environment(fun = func2)
    scale <- morph$scale.fun(scale)
    initial <- morph$transform(initial)
  
    out.time <- system.time(
    out <- .Call("metrop", func1, initial, nbatch, blen, nspac,
      scale, func2, debug, env1, env2)
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
    out$morph <- morph
    out$batch <- t(out$batch)
    out$debug <- debug
    
    if (! is.null(out$current)) out$current <- t(out$current)
    if (! is.null(out$proposal)) out$proposal <- t(out$proposal)
    if (! is.null(out$z)) out$z <- t(out$z)
    if (! is.null(out$morph)) out <- morph.set.object(out)
    
    class(out) <- c("mcmc", "metropolis")
    return(out)
}

