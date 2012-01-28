euclid.norm <- function(x) {
  sqrt(sum(x * x))
}

isotropic <- function(f) {
  f(1); # force evaluation of f
  function(x) {
    x.norm <- euclid.norm(x)
    if (x.norm == 0)
      rep(0, length(x))
    else
      f(x.norm) * x / x.norm
  }
}

isotropic.logjacobian <- function(f, d.f) {
  f(1); # force evaluation of components
  d.f(1);
  function(x) {
    x.norm <- euclid.norm(x)
    k <- length(x)
    if (x.norm == 0) {
      k * log(d.f(x.norm))
    } else {
      log(d.f(x.norm)) + (k - 1) * (log(f(x.norm)) - log(x.norm))
    }
  }
}

newton.raphson <- function(f, df, x, r) {
  x0 <- 2 * r
  while(abs(f(x0) - x) > sqrt(.Machine$double.eps)) {
    x0 <- x0 - (f(x0) - x) / df(x0)
  }
  return(x0)
}

subexponential <- function(b=1) {
  if (is.null(b)) b <- 1
  stopifnot(b > 0)
  b;
  f.inv <- function(x) ifelse(x > 1/b,
                              exp(b * x) - exp(1)/3,
                              (x * b)^3 * exp(1) / 6 + x * b * exp(1) / 2)
  d.f.inv <- function(x) ifelse(x > 1/b,
                                b * exp(b * x),
                                b * (x * b)^2 * exp(1) / 2 + b * exp(1) / 2)
  f <- function(x) {
    # x > exp(b * 1 / b) - exp(1) / 3
    if (x > 2 * exp(1) / 3) {
      log(x + exp(1)/3) / b
    } else {
      poly.inv <- exp(1/3) *
        (sqrt(b^12 * (9 * x^2 + exp(2))) - 3 * b^6 * x)^(-1/3)
      poly.inv * b - 1 / (poly.inv * b^3)
    }
  }
  return(list(f=f, f.inv=f.inv, d.f.inv=d.f.inv))
}

exponential <- function(r=1, p=3) {
  if (is.null(p)) p <- 3
  if (is.null(r)) r <- 0
  stopifnot(p > 2)
  stopifnot(r >= 0)
  f.inv <- function(x) x + (x-r)^p * (x > r)
  d.f.inv <- function(x) 1 + p * (x-r)^(p-1) * (x > r)
  if (p == 3) {
    g <- function(x) {
      n <- sqrt((27*r-27*x)^2 + 108) + 27 * (r - x)
      r + (2/n)^(1/3) - (n/2)^(1/3)/3
    }
    f <- function(x) ifelse(x < r, x, g(x))
  } else {
    # No general closed form solution exists.  However, since the
    # transformation has polynomial form, using the Newton-Raphson method
    # should work well.
    f <- function(x) ifelse(x < r,
                            x,
                            newton.raphson(f.inv, d.f.inv, x, r))
  }
  return(list(f=f, f.inv=f.inv, d.f.inv=d.f.inv))
}

morph.identity <- function() {
  morph(f=function(x) x,
        f.inv=function(x) x,
        logjacobian=function(x) 0)
}

morph <- function(f=NULL, f.inv=NULL, logjacobian=NULL,
                  scale.fun=NULL, scale.inv=NULL,
                  r=NULL, p=NULL, b=NULL,
                  center=0) {
  first.set <- c(is.null(f), is.null(f.inv), is.null(logjacobian))
  second.set <- c(is.null(r), is.null(p), is.null(b))
  if (!xor(all(first.set), all(second.set))) {
    stop(paste("Exactly one of the sets of arguments",
               "(f, f.inv, logjacobian) and (r, p, b)",
               "should be non NULL."))
  }
  if (any(first.set) & !all(first.set)) {
    stop("f, f.inv and logjacobian must all be specified.")
  }
  if (any(second.set)) {
    exp.f <- NULL
    sub.f <- NULL
    if (!is.null(r) | !is.null(p)) {
      exp.f <- exponential(r, p)
    }
    if (!is.null(b)) {
      sub.f <- subexponential(b)
    }
    if (!is.null(exp.f) & !is.null(sub.f)) {
      f <- function(x) isotropic(sub.f$f)(isotropic(exp.f$f)(x))
      f.inv <-
        function(x) isotropic(sub.f$f.inv)(isotropic(exp.f$f.inv)(x))
      logjacobian <-
        function(x)
          isotropic.logjacobian(sub.f$f.inv, sub.f$d.f.inv)(exp.f$f.inv(x)) +
            isotropic.logjacobian(exp.f$f.inv, exp.f$d.f.inv)(x)
    } else if (!is.null(exp.f)) {
      f <- isotropic(exp.f$f)
      f.inv <- isotropic(exp.f$f.inv)
      logjacobian <- isotropic.logjacobian(exp.f$f.inv, exp.f$d.f.inv)
    } else if (!is.null(sub.f)) {
      f <- isotropic(sub.f$f)
      f.inv <- isotropic(sub.f$f.inv)
      logjacobian <- isotropic.logjacobian(sub.f$f.inv, sub.f$d.f.inv)
    }
  }
  # force evaluation of transformation functions.
  #f; f.inv; logjacobian; center;
  if (missing(scale.fun)) scale.fun <- f
  if (missing(scale.inv)) scale.inv <- f.inv
  out <- list(f=f, f.inv=f.inv, logjacobian=logjacobian, center=center,
              scale.fun=scale.fun, scale.inv=scale.inv)
  transform <- function(state) out$f(state - out$center)
  inverse <- function(state) out$f.inv(state) + out$center
  transform; inverse;
  out$transform <- transform
  out$inverse <- inverse

  # pay attention to how '...' arguments are handled.  If care is not
  # taken, weird bugs will be introduced that will break handling
  # done by metrop.* functions.
  out$outfun <- function(outfun, ...) {
    list(...);
    if (is.null(outfun)) {
      return(function(state) out$inverse(state))
    } else if (is.function(outfun)) {
      outfun;
      return(function(state) outfun(out$inverse(state), ...))
    } else {
      return(function(state) out$inverse(state)[outfun])
    }
  }

  out$lud <- function(lud, ...) {
    lud; list(...);
    function(state)
      lud(out$inverse(state), ...) + out$logjacobian(state)
  }

  return(out)
}

morph.set.object <- function(out) {
  if (is.null(out) || is.null(out$morph)) {
    return(out)
  }
  out$scale   <- out$morph$scale.inv(out$scale)
  out$initial <- out$morph$inverse(out$initial)
  out$final   <- out$morph$inverse(out$final)
  # what about out${z,proposal,current}?
  return(out)
}
