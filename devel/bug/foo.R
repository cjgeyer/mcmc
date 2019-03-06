
# use example from metrop help page

library("mcmc", lib.loc = "../../package/mcmc.Rcheck")
packageVersion("mcmc")
set.seed(42)

h <- function(x) if (all(x >= 0) && sum(x) <= 1) return(0) else return(-Inf)
out <- metrop(h, rep(0, 5), 1000)
out$accept
# acceptance rate too low
out <- metrop(out, scale = 0.1)
t.test(out$accept.batch)$conf.int
# acceptance rate o. k. (about 25 percent)

# Now check that random numbers OK

out <- metrop(out, debug = TRUE)
names(out)

u <- out$u
u <- u[! is.na(u)]
length(u) == length(unique(u))

z <- out$z
dim(z)
length(z) == length(unique(z))

# try longer run?
out <- metrop(out, blen = 1000)

u <- out$u
u <- u[! is.na(u)]
length(u) == length(unique(u))
u[duplicated(u)]
seq(along = u)[duplicated(u)]

z <- out$z
dim(z)
length(z) == length(unique(z))

# Looks OK
# Try example from vignette

rm(list = ls())

data(logit)
out <- glm(y ~ x1 + x2 + x3 + x4, data = logit, family = binomial, x = TRUE)

# use function factory pattern to avoid ... and global variables

lupost_factory <- function(x, y) function(beta) {
    eta <- as.numeric(x %*% beta)
    logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
    logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
    logl <- sum(logp[y == 1]) + sum(logq[y == 0])
    return(logl - sum(beta^2) / 8)
}

lupost <- lupost_factory(out$x, out$y)

beta.init <- as.numeric(coefficients(out))

out <- metrop(lupost, beta.init, 1e3)
out$accept
# acceptance rate too low
out <- metrop(out, scale = 0.1)
out$accept
out <- metrop(out, scale = 0.3)
out$accept
out <- metrop(out, scale = 0.5)
out$accept
out <- metrop(out, scale = 0.4)
t.test(out$accept.batch)$conf.int
# acceptance rate o. k. (about 25 percent)

# Now check that random numbers OK

out <- metrop(out, debug = TRUE)
names(out)

u <- out$u
u <- u[! is.na(u)]
length(u) == length(unique(u))

z <- out$z
dim(z)
length(z) == length(unique(z))

# try longer run?
out <- metrop(out, blen = 1000)

u <- out$u
u <- u[! is.na(u)]
length(u) == length(unique(u))
u[duplicated(u)]
seq(along = u)[duplicated(u)]

z <- out$z
dim(z)
length(z) == length(unique(z))

