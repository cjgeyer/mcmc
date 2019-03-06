
set.seed(42)
RNGkind()

u <- runif(1e7)
length(u) == length(unique(u))
mean(duplicated(u))
