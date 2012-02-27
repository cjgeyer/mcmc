morph.metrop <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, b, r, p, center,
    ...)
UseMethod("morph.metrop")

morph.metrop.morph.metropolis <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, b, r, p, center, ...) {
  if (missing(b)) b <- obj$b
  if (missing(r)) r <- obj$r
  if (missing(p)) p <- obj$p
  if (missing(center)) center <- obj$center

  morph <- morph(b=b, r=r, p=p, center=center)
  obj$final <- obj$morph.final

  if (missing(outfun))
    outfun.save <- obj$outfun
  else
    outfun <- outfun
  if (missing(scale))
    scale.save <- obj$scale
  else
    scale.save <- scale

  if (missing(blen)) blen <- obj$blen
  if (missing(nspac)) nspac <- obj$nspac
  if (missing(debug)) debug <- obj$debug
  
  morphed.obj <- metrop.metropolis(obj,
                                   nbatch=nbatch,
                                   blen=blen,
                                   nspac=nspac,
                                   scale=scale.save,
                                   outfun=morph$outfun(outfun.save),
                                   debug=debug,
                                   ...)
  
  unmorphed.obj <- .morph.unmorph(morphed.obj, morph, outfun.save,
                                  scale.save, b, r, p, center)
  return(unmorphed.obj)
}

morph.metrop.function <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, b, r, p, center, ...) {

  # TODO move these missing statements to morph
  morph <- morph(b=b, r=r, p=p, center=center)
  if (missing(outfun)) outfun <- NULL
  outfun.save <- outfun
  scale.save <- scale
  
  morphed.obj <- metrop.function(morph$lud(obj),
                                 initial=morph$transform(initial),
                                 nbatch=nbatch,
                                 blen=blen,
                                 scale=scale,
                                 outfun=morph$outfun(outfun),
                                 debug=debug,
                                 ...)
  
  unmorphed.obj <- .morph.unmorph(morphed.obj, morph, outfun.save,
                                  scale.save, b, r, p, center)
  return(unmorphed.obj)
}

.morph.unmorph <- function(obj, morph, outfun, scale, b, r, p, center) {
  if (missing(b)) b <- NULL
  if (missing(r)) r <- NULL
  if (missing(p)) p <- NULL
  if (missing(center)) center <- NULL
  obj$morph       <- morph
  obj$morph.final <- obj$final
  obj$final       <- morph$inverse(obj$final)
  obj$outfun      <- outfun
  obj$scale       <- scale
  obj$b           <- b
  obj$r           <- r
  obj$p           <- p
  obj$center      <- center
  class(obj) <- c("mcmc", "morph.metropolis")
  return(obj)
}
