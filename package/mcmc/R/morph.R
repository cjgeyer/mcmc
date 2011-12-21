isotropic <- function(f) {
  f(1); # force evaluation of f
  function(x) {
    x.norm <- sqrt(sum(x * x))
    if (x.norm == 0)
      rep(0, length(x))
    else
      f(x.norm) * x / x.norm
  }
}

exponential.superexponetial <- function(p, r) {
  stopifnot(p > 2)
  stopifnot(r >= 0)
  f.inv <- isotropic(function(x) ifelse(x < r, x, x + (x-r)^p))
  dee.f.inv <- function(x) ifelse(x < r, 1, 1 + p * (x-r)^(p-1))
  if (p == 3) {
    g <- function(x) {
      n <- sqrt((27*r-27*x)^2 + 108) + 27 * (r - x)
      r + (2/n)^(1/3) - (n/2)^(1/3)/3
    }
    f <- isotropic(function(x) ifelse(x < r, x, g(x)))
  } else {
    # sadly, no general closed form solutions exist.  However, sine p > 2
    # and R >= 0, x + (x-R)^p is strictly increasing.  A simple binary
    # search will work.
    f <- NOTWRITTENYET
  }
  morph.isotropic(f, f.inv, dee.f.inv)
}

morph.isotropic <- function(f, f.inv, dee.f.inv) {
  logdetdeeh <- function(x) {
    x.norm <- sqrt(sum(x * x))
    if (x.norm == 0)
      length(x) * log(dee.f.inv(0))
    else
      log(dee.f.inv(x.norm)) + (k - 1)*(log(f(x.norm)) - log(x.norm))
  }
  morph(isotropic(f), isotropic(f.inv), logdetdeeh)
}

morph <- function(f, f.inv, logdetdeeh) {
  # is forcing f, f.inv, logdeeh necessary?
  out <- list()
  out$outfun <- function(outfun) {
    if (is.null(outfun)) {
      return(function(x, ...) f.inv(x))
    } else if (is.function(outfun)) {
      return(function(x, ...) outfun(f.inv(x), ...))
    } else {
      return(function(x, ...) f.inv(x)[outfun])
    }
  }

  out$transform <- function(x) f(x)
  out$inverse <- function(x) f.inv(x)
  out$lud <- function(lud) {
    function(lud, ...) lud(f.inv(x), ...) + logdetdeeh(x)
  }

  return(out)
}


morph.set.object <- function(out, morph=NULL) {
  if (is.null(morph) || is.null(out)) {
    return(out)
  }
  out$scale   <- morph$inverse(out$scale)
  out$initial <- morph$inverse(out$initial)
  out$final   <- morph$inverse(out$final)
  # what about out${z,proposal,current}?
  return(out)
}
