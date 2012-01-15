library(mcmc)

# make sure morph.identity works properly
ident.func <- function(x) x
m.id <- morph.identity()
all.equal(m.id$transform(1:10), 1:10)
all.equal(m.id$inverse(1:10), 1:10)
x <- seq(-1,1, length.out=15)
all.equal(sapply(x, m.id$lud(function(x) dnorm(x, log=TRUE))),
          dnorm(x, log=TRUE))
all.equal(m.id$outfun(ident.func)(x), x)

# test center parameter, univariate version
zero.func <- function(x) 0
center <- 2
morph.center <- morph(f=ident.func,
                      f.inv=ident.func,
                      logjacobian=zero.func,
                      center=center)
all.equal(morph.center$transform(x), x-center)
all.equal(morph.center$inverse(x), x+center)
all.equal(sapply(x, morph.center$lud(function(y) dnorm(y, log=TRUE))),
          dnorm(x, log=TRUE))
# test center parameter, multivariate version
center <- 1:4
x <- rep(0, 4)
morph.center <- morph(f=ident.func,
                      f.inv=ident.func,
                      logjacobian=zero.func,
                      center=center)
lud.mult.dnorm <- function(x) prod(dnorm(x, log=TRUE))
all.equal(morph.center$transform(x), x-center)
all.equal(morph.center$inverse(x), x+center)
all.equal(morph.center$lud(lud.mult.dnorm)(x),
          lud.mult.dnorm(x - center))
# test 'r'.
r <- 1
morph.r <- morph(r=r)
x <- seq(-1, 1, length.out=20)
all.equal(sapply(x, morph.r$lud(function(x) dnorm(x, log=TRUE))),
          dnorm(x, log=TRUE))
x <- seq(1.1, 2, length.out=10)
all(sapply(x, morph.r$lud(function(x) dnorm(x, log=TRUE)))
    !=
    dnorm(x, log=TRUE))

# make sure '...' arguments work properly with lazy-loading
mean <- 2
ident.morph <- morph.identity()
dnorm.morph <- ident.morph$lud(function(x, mean=0)
                                 dnorm(x, mean=mean, log=TRUE),
                               mean=mean)
mean <- 3
all.equal(dnorm.morph(2), dnorm(2, mean=2, log=TRUE))

# check lazy-loading of '...' arguments by outfun
outfun.orig <- function(x, mean) x + mean
ident.morph <- morph.identity()
mean <- 1
outfun.morph <- ident.morph$outfun(outfun.orig, mean=mean)
mean <- 2
all.equal(outfun.morph(1:10), 1:10+1)

# make sure morph$outfun passes '...' arguments.
outfun.orig <- function(x, mean) x + mean
ident.morph <- morph.identity()
mean <- 1
outfun.morph <- ident.morph$outfun(outfun.orig, mean=mean)
all.equal(outfun.morph(1:10), 1:10+mean)

# check lazy-loading of f, f.inv, logjacobian and center
f <- function(x) x^2
f.inv <- function(x) sqrt(x)
d.f.inv <- function(x) 1/2 * x^(-1/2)
logjacobian <- isotropic.logjacobian(f.inv, d.f.inv)
center <- 3
morph.test <- morph(f=isotropic(f),
                    f.inv=isotropic(f.inv),
                    logjacobian=logjacobian,
                    center=center)
f <- function(x) x
f.inv <- function(x) x
logjacobian <- function(x) 0
center <- 2
x <- seq(1, 3, length.out=20)
# are f and f.inv the passed version?
all.equal(sapply(x, morph.test$f), x^2)
all.equal(sapply(x, morph.test$f.inv), sqrt(x))
# is centered still the passed value?
x.centered <- x - 3
norm.x.centered <- sqrt(x.centered * x.centered)
all.equal(0, morph.test$transform(3))
all.equal(3, morph.test$inverse(0))
# make sure that transformation/inverse weren't harmed.
all.equal(sapply(x, morph.test$transform),
          norm.x.centered * x.centered)
all.equal(sapply(x, morph.test$inverse),
          sqrt(x) + 3)
# does lud use the passed values (this uses logjacobian and f.inv).
all.equal(sapply(x, morph.test$lud(dnorm, log=TRUE)),
          dnorm(sqrt(x) + 3, log=TRUE) +
          sapply(x, isotropic.logjacobian(sqrt,
                                          function(x) 1/2 * x^(-1/2))))

          
