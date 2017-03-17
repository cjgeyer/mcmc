
# new feature batching acceptance rates

 set.seed(42)

 library(mcmc)

 h <- function(x) if (all(x >= 0) && sum(x) <= 1) return(1) else return(-Inf)
 out <- metrop(h, rep(0, 5), nbatch = 100, blen = 100, scale = 0.1,
     debug = TRUE)

 all.equal(out$accept, mean(out$accept.batch))

 foo <- matrix(out$debug.accept, nrow = out$blen)
 bar <- colMeans(foo)
 all.equal(out$accept.batch, bar)

 options(digits = 4) # try to keep insanity of computer arithmetic under control

 out$accept
 t.test(out$accept.batch)$conf.int

