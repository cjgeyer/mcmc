
library(mcmc, lib.loc = "../package/mcmc.Rcheck")

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

hout <- hitrun(function(x) return(0), rep(0, d), nbatch = 100,
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

all(hout$log.green == 0)

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
identical(x.initial, hout$rip)
identical(x.initial, hout$current[1, ])

mys1 <- rep(Inf, nrow(hout$current))
mys2 <- rep(-Inf, nrow(hout$current))
for (i in seq(along = myu1)) {
   z <- myz[i, ]
   x <- hout$current[i, ]
   ax <- as.numeric(amat %*% x)
   az <- as.numeric(amat %*% z)
   bmax <- bvec - ax
   mys1[i] <- min(ifelse(az > 0, bmax / az, Inf))
   mys2[i] <- max(ifelse(az < 0, bmax / az, -Inf))
}
identical(mys1, hout$s1)
identical(mys2, hout$s2)
head(mys1)
head(hout$s1)
head(mys2)
head(hout$s2)

dim(hout$batch)
dim(a1)
length(b1)
foo <- hout$batch %*% t(a1)
foo <- sweep(foo, 2, b1)
foo <- apply(foo, 1, max)
max(foo)

