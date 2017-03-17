
 library(mcmc, lib.loc = "mcmc.Rcheck")
 library(MASS)

 set.seed(42)

 d <- 9
 witch.which <- c(0.1, 0.3, 0.5, 0.7, 1.0)
 ncomp <- length(witch.which)

 neighbors <- matrix(FALSE, ncomp, ncomp)
 neighbors[row(neighbors) == col(neighbors) + 1] <- TRUE
 neighbors[row(neighbors) == col(neighbors) - 1] <- TRUE

 ludfun <- function(state, log.pseudo.prior) {
     stopifnot(is.numeric(state))
     stopifnot(length(state) == d + 1)
     icomp <- state[1]
     stopifnot(icomp == as.integer(icomp))
     stopifnot(1 <= icomp && icomp <= ncomp)
     stopifnot(is.numeric(log.pseudo.prior))
     stopifnot(length(log.pseudo.prior) == ncomp)
     theta <- state[-1]
     if (any(theta > 1.0)) return(-Inf)
     bnd <- witch.which[icomp]
     lpp <- log.pseudo.prior[icomp]
     if (any(theta > bnd)) return(lpp)
     return(- d * log(bnd) + lpp)
 }

 theta.initial <- c(1, rep(0.5, d))
 qux <- rep(0, ncomp)

 out <- temper(ludfun, initial = theta.initial, neighbors = neighbors,
     nbatch = 50, blen = 30, nspac = 2, scale = 0.56789,
     parallel = FALSE, debug = FALSE, log.pseudo.prior = qux)

 repeat {

     out <- temper(out, nbatch = 200, blen = 100, nspac = 1,
         log.pseudo.prior = qux)
     foo <- apply(out$ibatch, 2, mean)
     bar <- qux - log(foo)
     bar <- bar - min(bar)
     baz <- bar > qux + log(100)
     bar[baz] <- (qux + log(100))[baz]
     print(bar)

     qux <- bar

     if(all(baz == FALSE)) break
 }

 out <- temper(out, nbatch = 1000, log.pseudo.prior = qux)

 foo <- apply(out$ibatch, 2, mean)
 foo
 sqrt(apply(out$ibatch, 2, var) / out$nbatch)
 sqrt(apply(out$ibatch, 2, function(x) initseq(x)$var.con) / out$nbatch)

 bar <- qux - log(foo)
 bar <- bar - min(bar)
 baz <- bar > qux + log(100)
 bar[baz] <- (qux + log(100))[baz]
 print(bar)

