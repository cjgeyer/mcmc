
 library(mcmc)

# should give intelligible error (unlike before ver 0.9-8)

 suppressMessages(try(metrop(function(x) x, double(0), nbatch = 10)))
