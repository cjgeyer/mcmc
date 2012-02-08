morph.metrop <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, b, r, p, center=0,
    ...)
UseMethod("morph.metrop")

morph.metrop.morph.metropolis <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, b, r, p, center=0, ...) {
  if (any(!missing(b), !missing(r), !missing(p), !missing(center))) {
    obj$morph <- morph(b=b, r=r, p=p, center=center)
  }
  morph <- obj$morph
  
  outfun.save <- ifelse(missing(outfun), obj$outfun, outfun)
  scale.save <- ifelse(missing(scale), obj$scale, scale)
  
  morphed.obj <- metrop.metropolis(morphed.obj,
                                   nbatch=nbatch,
                                   blen=blen,
                                   nspac=nspac,
                                   scale=morph$scale.fun(scale.save),
                                   outfun=morph$outfun(outfun.save),
                                   debug=debug,
                                   ...)
  
  unmorphed.obj <- .morph.unmorph(morphed.obj, morph, outfun.save,
                                  scale.save)
  return(unmorphed.obj)
}

morph.metrop.function <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, b, r, p, center=0, ...) {

  if (missing(b)) b <- NULL
  if (missing(r)) r <- NULL
  if (missing(p)) p <- NULL
  morph <- morph(b=b, r=r, p=p)
  #morph <- morph(f=c.sub.f, f.inv=c.sub.finv, logjacobian=c.sub.logj)
  if (missing(outfun)) outfun <- NULL
  outfun.save <- outfun
  scale.save <- scale
  
  metrop.out <- metrop.function(morph$lud(obj),
                                initial=morph$transform(initial),
                                nbatch=nbatch,
                                blen=blen,
                                scale=morph$scale.fun(scale),
                                outfun=morph$outfun(outfun),
                                debug=debug,
                                ...)
  
  unmorphed.obj <- .morph.unmorph(metrop.out, morph, outfun.save,
                                  scale.save)
  return(unmorphed.obj)
}

.morph.unmorph <- function(obj, morph, outfun, scale) {
  obj$final <- morph$inverse(obj$final)
  obj$outfun <- outfun
  obj$scale <- scale
  class(obj) <- c("mcmc", "morph.metropolis")
  return(obj)
}

morph.metropC <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, b, r, p, center=0, ...) {

  if (missing(b)) b <- NULL
  if (missing(r)) r <- NULL
  if (missing(p)) p <- NULL
  #morph <- morph(b=b, r=r, p=p)
  morph <- morph(f=c.sub.f, f.inv=c.sub.finv, logjacobian=c.sub.logj)
  if (missing(outfun)) outfun <- NULL
  outfun.save <- outfun
  scale.save <- scale
  
  metrop.out <- metrop.function(morph$lud(obj),
                                initial=morph$transform(initial),
                                nbatch=nbatch,
                                blen=blen,
                                scale=morph$scale.fun(scale),
                                outfun=morph$outfun(outfun),
                                debug=debug,
                                ...)
  
  unmorphed.obj <- .morph.unmorph(metrop.out, morph, outfun.save,
                                  scale.save)
  return(unmorphed.obj)
}
