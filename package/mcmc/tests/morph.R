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

# make sure morph$outfun passes '...' arguments.
outfun.orig <- function(x, mean) x + mean
ident.morph <- morph.identity()
mean <- 1
outfun.morph <- ident.morph$outfun(outfun.orig, mean=mean)
all.equal(outfun.morph(1:10), 1:10+mean)

# check lazy-loading of '...' arguments by outfun
outfun.orig <- function(x, mean) x + mean
ident.morph <- morph.identity()
mean <- 1
outfun.morph <- ident.morph$outfun(outfun.orig, mean=mean)
mean <- 2
all.equal(outfun.morph(1:10), 1:10+1)
