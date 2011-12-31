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

isotropic.logjacobian <- function(f.inv, d.f.inv) {
  f.inv(1); # force evaluation of components
  d.f.inv(1);
  function(x) {
    x.norm <- euclid.norm(x)
    k <- length(x)
    if (x.norm == 0) {
      k * log(d.f.inv(x.norm))
    } else {
      log(d.f.inv(x.norm)) + (k - 1) * (log(f.inv(x.norm)) - log(x.norm))
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
  b
  f.inv <- function(x) ifelse(x > 1/b,
                              exp(b * x) - exp(1)/3,
                              (x * b)^3 * exp(1) / 6 + x * b * exp(1) / 2)
  d.f.inv <- function(x) ifelse(x > 1/b,
                                n * exp(b * x),
                                (x * b)^2 * exp(1) / 2 + b * exp(1) / 2)
  f <- function(x) {
    if (x > 2 * exp(1) / 3) {
      log(x + exp(1)/3) / b
    } else {
      g <- (sqrt(b^12 * (9*x^2 + exp(2)) - 3 * b^6 * x))^(1/3)
      exp(1/3) * b / g - g / (exp(1/3) * b^3)
    }
  }
  return(list(f=f, f.inv=f.inv, d.f.inv=d.f.inv))
}

exponential <- function(r=1, p=3) {
  if (is.null(p)) p <- 3
  if (is.null(r)) r <- 1
  stopifnot(p > 2)
  stopifnot(r >= 0)
  f.inv <- function(x) x + (x-r)^p * (x > r)
  dee.f.inv <- function(x) 1 + p * (x-r)^(p-1) * (x > r)
  if (p == 3) {
    g <- function(x) {
      n <- sqrt((27*r-27*x)^2 + 108) + 27 * (r - x)
      r + (2/n)^(1/3) - (n/2)^(1/3)/3
    }
    f <- function(x) ifelse(x < r, x, g(x))
  } else {
    # No general closed form solution exists.  However, since the transformation
    # has polynomial form, using the Newton-Raphson method should work well.
    f <- function(x) ifelse(x < r, x, newton.raphson(f.inv, dee.f.inv, x, r))
  }
  return(list(f=f, f.inv=f.inv, d.f.inv=d.f.inv))
}

morph <- function(f=NULL, f.inv=NULL, logjacobian=NULL,
                  r=NULL, p=NULL, b=NULL,
                  center=0) {
  first.set <- c(is.null(f), is.null(f.inv), is.null(logjacobian))
  second.set <- c(is.null(r), is.null(p), is.null(b))
  if (!xor(any(first.set), any(second.set))) {
    stop("Exactly one of the sets of arguments (f, f.inv, logjacobian) and (r, p, b) should be non NULL.")
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
  # is forcing f, f.inv, logdeeh necessary?
  out <- list()
  out$outfun <- function(outfun) {
    if (is.null(outfun)) {
      return(function(x, ...) f.inv(x) + center)
    } else if (is.function(outfun)) {
      return(function(x, ...) outfun(f.inv(x) + center, ...))
    } else {
      return(function(x, ...) (f.inv(x) + center)[outfun])
    }
  }

  out$transform <- function(x) f(x - center)
  out$inverse <- function(x) f.inv(x) + center
  out$lud <- function(lud) {
    function(lud, ...) lud(f.inv(x) + center, ...) + logjacobian(x)
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
