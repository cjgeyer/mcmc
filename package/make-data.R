
 library(MASS)

 RNGversion("2.9.1")
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

 foo <- data.frame(x1 = x1, x2 = x2, x3 = x3, y)

 write.table(foo, file = "foo.txt", row.names = FALSE, quote = FALSE)

