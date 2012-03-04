morph.metrop <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, morph,
    ...)
UseMethod("morph.metrop")

morph.metrop.morph.metropolis <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, morph, ...) {
  if (missing(morph)) morph <- obj$morph
  obj$final <- obj$morph.final

  if (missing(outfun))
    outfun.save <- obj$outfun
  else
    outfun.save <- outfun
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
  
  unmorphed.obj <- .morph.unmorph(morphed.obj, morph, outfun.save, scale.save)
  return(unmorphed.obj)
}

morph.metrop.function <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, morph, ...) {

  if (missing(morph)) morph <- morph.identity()
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
  
  unmorphed.obj <- .morph.unmorph(morphed.obj, morph, outfun.save, scale.save)
  return(unmorphed.obj)
}

.morph.unmorph <- function(obj, morph, outfun, scale) {
  obj$morph       <- morph
  obj$morph.final <- obj$final
  obj$final       <- morph$inverse(obj$final)
  obj$outfun      <- outfun
  obj$scale       <- scale
  class(obj) <- c("mcmc", "morph.metropolis")
  return(obj)
}
