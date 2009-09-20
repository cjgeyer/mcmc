
 library(mcmc)

 set.seed(42)

 d <- 9
 witch.which <- c(0.1, 0.3, 0.5, 0.7, 1.0)
 ncomp <- length(witch.which)

 neighbors <- matrix(FALSE, ncomp, ncomp)
 neighbors[row(neighbors) == col(neighbors) + 1] <- TRUE
 neighbors[row(neighbors) == col(neighbors) - 1] <- TRUE

 ludfun <- function(state) {
     stopifnot(is.numeric(state))
     stopifnot(length(state) == d + 1)
     icomp <- state[1]
     stopifnot(icomp == as.integer(icomp))
     stopifnot(1 <= icomp && icomp <= ncomp)
     theta <- state[-1]
     if (any(theta > 1.0)) return(-Inf)
     bnd <- witch.which[icomp]
     if (any(theta > bnd)) return(0)
     return(- d * log(bnd))
 }

 thetas <- matrix(0.5, ncomp, d)
 out <- temper(ludfun, initial = thetas, neighbors = neighbors, nbatch = 20,
     blen = 10, nspac = 5, scale = 0.56789, parallel = TRUE, debug = TRUE)

 names(out)

 out$acceptx

 out$accepti

 ### seems to be o. k.

