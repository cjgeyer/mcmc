library(mcmc)

# make sure morph.identity works properly
id <- function(x) x
m.id <- morph.identity()
all.equal(m.id$transform(1:10), 1:10)
all.equal(m.id$inverse(1:10), 1:10)
x <- seq(-1,1, length.out=15)
all.equal(sapply(x, m.id$lud(function(x) dnorm(x, log=TRUE))),
          dnorm(x, log=TRUE))
all.equal(m.id$outfun(id)(x), x)
