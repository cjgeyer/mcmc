
 library(mcmc, lib.loc = "mcmc.Rcheck")
 library(MASS)

 set.seed(42)

 n <- 100
 rho <- 0.5
 beta0 <- 0.25
 beta1 <- 0.5
 beta2 <- 1
 beta3 <- 1.5

 Sigma <- matrix(rho, 3, 3)
 diag(Sigma) <- 1
 Sigma <- 0.75 * Sigma
 Mu <- rep(0, 3)

 foo <- mvrnorm(n, Mu, Sigma)

 x1 <- foo[ , 1]
 x2 <- foo[ , 2]
 x3 <- foo[ , 3]

 modmat <- cbind(1, foo)

 eta <- beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3
 p <- 1 / (1 + exp(- eta))
 y <- as.numeric(runif(n) < p)

 out <- glm(y ~ x1 + x2 + x3, family = binomial())
 summary(out)

 models <- cbind(rep(0:1, each = 4), rep(rep(0:1, times = 2), each = 2),
               rep(0:1, times = 4))

 exes <- paste("x", 1:3, sep = "")
 models[nrow(models), ]
 beta.initial <- c(nrow(models), out$coefficients)

 neighbors <- matrix(FALSE, nrow(models), nrow(models))
 for (i in 1:nrow(neighbors)) {
     for (j in 1:ncol(neighbors)) {
         foo <- models[i, ]
         bar <- models[j, ]
         if (sum(foo != bar) == 1) neighbors[i, j] <- TRUE
     }
 }
 neighbors

 ludfun <- function(state, log.pseudo.prior) {
     stopifnot(is.numeric(state))
     stopifnot(length(state) == ncol(models) + 2)
     stopifnot(length(state) == ncol(models) + 2)
     icomp <- state[1]
     stopifnot(icomp == as.integer(icomp))
     stopifnot(1 <= icomp && icomp <= nrow(models))
     stopifnot(is.numeric(log.pseudo.prior))
     stopifnot(length(log.pseudo.prior) == nrow(models))
     beta <- state[-1]
     inies <- c(TRUE, as.logical(models[icomp, ]))
     beta.logl <- beta
     beta.logl[! inies] <- 0
     eta <- as.numeric(modmat %*% beta.logl)
     logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
     logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
     logl <- sum(logp[y == 1]) + sum(logq[y == 0])
     val <- logl - sum(beta^2) / 2 + log.pseudo.prior[icomp]
     return(val)
 }

 qux <- rep(0, nrow(models))

 out <- temper(ludfun, initial = beta.initial, neighbors = neighbors,
     nbatch = 20, blen = 10, nspac = 5, scale = 0.56789, debug = TRUE,
     log.pseudo.prior = qux)

 repeat {

     out <- temper(out, nbatch = 200, blen = 100, nspac = 1, debug = FALSE,
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

