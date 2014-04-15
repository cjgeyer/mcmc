
library(mcmc)

set.seed(42)

d <- 10
r <- floor((d - 1) / 4)

i <- 1:d
a2 <- rbind(rep(1, d), i, as.numeric(i <= (d + 1) / 2) -
    as.numeric(i >= (d + 1) / 2))
b2 <- c(1, (d + 1) / 2, 0)
a1 <- - diag(d)
b1 <- rep(0, d)
a1 <- rbind(a1, - as.numeric(i <= (d + 1) / 2 - r) -
    as.numeric(i >= (d + 1) / 2 + r))
b1 <- c(b1, - 1 / 2)
dimnames(a2) <- NULL

hout <- hitrun(function(x) return(0), nbatch = 100,
    a1 = a1, b1 = b1, a2 = a2, b2 = b2, debug = TRUE)
names(hout)

# check random numbers
set.seed(42)
identical(.Random.seed, hout$initial.seed)
myz <- matrix(NA_real_, nrow(hout$current), ncol(hout$current))
myu1 <- rep(NA_real_, nrow(hout$current))
myu2 <- rep(NA_real_, nrow(hout$current))
for (i in seq(along = myu1)) {
    myz[i, ] <- rnorm(ncol(hout$current))
    myu1[i] <- runif(1)
}
identical(.Random.seed, hout$final.seed)
identical(myz, hout$z)
identical(myu1, hout$u1)
identical(myu2, hout$u2)

# check green ratios (trivial for this problem)
all(hout$log.green == 0)

# check construction of basis, origin, amat, bvec, rip
hrep4 <- makeH(a1 = a1, b1 = b1)
hrep3 <- makeH(a2 = a2, b2 = b2)
hrep4 <- d2q(hrep4)
hrep3 <- d2q(hrep3)
dim(hrep3)
dim(hrep4)
vrep3 <- scdd(hrep3)$output
is.line <- vrep3[ , 1] == "1" & vrep3[ , 2] == "0"
is.point <- vrep3[ , 1] == "0" & vrep3[ , 2] == "1"
all(is.point | is.line)
sum(is.point) == 1
sum(is.line) == d - nrow(hrep3)
foo <- vrep3[ , - c(1, 2)]
origin <- foo[is.point, ]
basis <- foo[is.line, ]
basis <- t(basis)
a1q <- d2q(a1)
b1q <- d2q(b1)
amat <- qmatmult(a1q, basis)
bvec <- qmq(b1q, qmatmult(a1q, cbind(origin)))
hrep5 <- makeH(a1 = amat, b1 = bvec)
vrep5 <- scdd(hrep5)$output
is.point <- vrep5[ , 1] == "0" & vrep5[ , 2] == "1"
all(is.point)
v <- vrep5[ , - c(1, 2)]
x <- apply(v, 2, qsum)
x <- qdq(x, rep(as.character(nrow(v)), length(x)))

origin <- q2d(origin)
basis <- q2d(basis)
amat <- q2d(amat)
bvec <- q2d(bvec)
x.initial <- q2d(x)
identical(origin, hout$origin)
identical(basis, hout$basis)
identical(amat, hout$amat)
identical(bvec, hout$bvec)
identical(x.initial, hout$initial)
identical(x.initial, hout$current[1, ])

# check bounds
mys1 <- rep(Inf, nrow(hout$current))
mys2 <- rep(-Inf, nrow(hout$current))
for (i in seq(along = myu1)) {
   z <- hout$z[i, ]
   x <- hout$current[i, ]
   ax <- as.numeric(amat %*% x)
   az <- as.numeric(amat %*% z)
   bnd <- (bvec - ax) / az
   mys1[i] <- max(bnd[az < 0])
   mys2[i] <- min(bnd[az > 0])
}
all.equal(mys1, hout$s1)
all.equal(mys2, hout$s2)

# check proposal
myproposal <- matrix(NA, nrow(hout$current), ncol(hout$current))
for (i in seq(along = myu1)) {
   x <- hout$current[i, ]
   z <- hout$z[i, ]
   smin <- hout$s1[i]
   smax <- hout$s2[i]
   u <- hout$u1[i]
   myproposal[i, ] <- x + z * (u * smin + (1 - u) * smax)
}
all.equal(myproposal, hout$proposal)
identical(hout$current[- 1, ], hout$proposal[- nrow(hout$proposal), ])
identical(hout$current[1, ], hout$initial)
identical(hout$proposal[nrow(hout$proposal), ], hout$final)

# check path is feasible (reduced coordinates)
foo <- hout$current %*% t(amat)
foo <- sweep(foo, 2, bvec)
all(foo <= 0)

# check transformation
foo <- hout$proposal %*% t(basis)
foo <- sweep(foo, 2, origin, "+")
identical(foo, hout$batch)

# check path is feasible (original coordinates)
foo <- hout$batch %*% t(a1)
foo <- sweep(foo, 2, b1)
bar <- hout$batch %*% t(a2)
bar <- sweep(bar, 2, b2)
all(foo <= sqrt(.Machine$double.eps))
all(abs(bar) <= sqrt(.Machine$double.eps))

# everything checks when ludfun is flat, now for non-flat
ludfun <- function(x) {
    stopifnot(is.numeric(x))
    stopifnot(is.finite(x))
    stopifnot(length(x) == d)
    1.3 * sum(log(x))
}
hout <- hitrun(ludfun, nbatch = 100,
    a1 = a1, b1 = b1, a2 = a2, b2 = b2, debug = TRUE)

foo <- hout$current %*% t(basis)
foo <- sweep(foo, 2, origin, "+")
bar <- hout$proposal %*% t(basis)
bar <- sweep(bar, 2, origin, "+")
my.log.green <- apply(bar, 1, ludfun) - apply(foo, 1, ludfun)
identical(my.log.green, hout$log.green)
identical(is.na(hout$u2), hout$log.green >= 0)
my.accept.1 <- is.na(hout$u2) | hout$u2 < exp(hout$log.green)
foo <- hout$proposal
bar <- hout$current[- 1, ]
bar <- rbind(bar, hout$final)
baz <- foo - bar
my.accept.2 <- apply(baz == 0, 1, all)
identical(my.accept.1, my.accept.2)

# now check restart property

.Random.seed <- hout$initial.seed
hout1 <- hitrun(ludfun, nbatch = 50, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
hout2 <- hitrun(hout1)
identical(hout$batch, rbind(hout1$batch, hout2$batch))
identical(hout$final, hout2$final)
identical(hout$final.seed, hout2$final.seed)

# now check batching and spacing

hout3 <- hitrun(hout, nbatch = 17, nspac = 3, blen = 11)
hout4 <- hitrun(hout, nspac = 1, blen = 1,
    nbatch = hout3$nbatch * hout3$blen * hout3$nspac)
mybatch <- hout4$batch[seq(1, hout4$nbatch) %% hout3$nspac == 0, ]
dim(mybatch)
mybatch <- array(as.vector(mybatch),
    c(hout3$blen, hout3$nbatch, ncol(mybatch)))
dim(mybatch)
mybatch <- apply(mybatch, c(2, 3), mean)
dim(mybatch)
all.equal(mybatch, hout3$batch)

