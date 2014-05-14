
 library(rcdd)

 # problem as specified by the user
 # median = (d + 1) / 2
 # mean  >= (d + 1) / 2
 d <- 10
 stopifnot(d %% 2 == 0)
 x <- 1:d
 a1 <- rbind(- x)
 b1 <- (- (d + 1) / 2)
 a2 <- matrix(0, 1, d)
 a2[ , x < (d + 1) / 2] <- 1
 b2 <- 1 / 2
 # maximum change in adjacent probabilities is 1 / 16
 foo <- matrix(0, 2 * (d - 1), d)
 for (i in 1:(d - 1)) {
     foo[2 * i - 1, i] <- 1
     foo[2 * i - 1, i + 1] <- - 1
     foo[2 * i, i] <- - 1
     foo[2 * i, i + 1] <- 1
 }
 bar <- rep(1 / 16, 2 * (d - 1))
 a1 <- rbind(a1, foo)
 b1 <- c(b1, bar)

 # problem as specified by the computer
 a1 <- rbind(- diag(d), a1)
 b1 <- c(rep(0, d), b1)
 a2 <- rbind(1, a2)
 b2 <- c(1, b2)

 # H-representation
 # to keep track of nonnegativity constraints, don't use makeH
 hrep1 <- rbind(cbind(0, b1, - a1), cbind(1, b2, - a2))
 dimnames(hrep1) <- NULL
 hrep1 <- d2q(hrep1)
 hrep1

 # V-representation
 vrep1 <- scdd(hrep1, representation = "H")$output
 # check that there are at least 2 vertices
 nrow(vrep1) >= 2
 apply(vrep1[ , - c(1, 2)], 2, qmin)
 # round(q2d(vrep1[ , - c(1, 2)]), 4)

 # find point in the relative interior of the constraint set
 v1 <- vrep1[ , - c(1, 2)]
 x1 <- apply(v1, 2, function(x) qdq(qsum(x), d2q(length(x))))
 x1

 # non-redundant H-representation
 rout <- redundant(hrep1)
 names(rout)
 hrep2 <- rout$output
 # keep track of where nonnegativity constraints are now
 nonneg.position <- rout$new.position[1:d]
 nonneg.position
 # only variables with nonneg.position != 0 can be zero somewhere in C

 # affine hull
 is.equality <- hrep2[ , 1] == "1"
 hrep3 <- hrep2[is.equality, , drop = FALSE]
 hrep3
 vrep3 <- scdd(hrep3, representation = "H")$output
 vrep3
 # 0 1 in first column is point
 # 0 0 in first column is ray
 # 1 0 in first column is line
 # 1 1 in first column is affine generator
 is.point <- vrep3[ , 1] == "0" & vrep3[ , 2] == "1"
 is.line <- vrep3[ , 1] == "1" & vrep3[ , 2] == "0"
 all(is.point | is.line)
 basis <- vrep3[is.line, , drop = FALSE]
 basis <- basis[ , - c(1, 2), drop = FALSE]
 basis <- t(basis)
 basis
 origin <- x1
 origin
 p <- ncol(basis)
 p
 # map from new coordinates (NC) to original coordinates (OC) is
 # function(w) as.numeric(origin + basis %*% w)

 # constraints in new coordinates
 hrep4 <- hrep2[! is.equality, , drop = FALSE]
 b4 <- hrep4[ , 2]
 a4 <- qneg(hrep4[ , - c(1, 2), drop = FALSE])
 a4
 b4
 amat <- qmatmult(a4, basis)
 bvec <- qmq(b4, qmatmult(a4, cbind(origin)))
 amat
 bvec
 # constraints in NC are amat %*% w <= bvec
 hrep5 <- makeH(amat, bvec)
 vrep5 <- scdd(hrep5)$output
 v5 <- vrep5[ , - c(1, 2)]

 # check, not done in condir function
 # do all vertices in v5 match up with vertices in v1?
 foo <- t(apply(v5, 1, function(w) qpq(origin, qmatmult(basis, cbind(w)))))
 nrow(foo) == nrow(v1)
 foo.too <- apply(foo, 1, paste, collapse = " ")
 bar.too <- apply(v1, 1, paste, collapse = " ")
 length(unique(foo.too)) == length(foo.too)
 length(unique(bar.too)) == length(bar.too)
 identical(sort(foo.too), sort(bar.too))

 # now where are the nonnegativity constraints?  and which are which?
 # we know where the are in hrep2
 nonneg.position
 nonneg.position.varb <- seq(along = nonneg.position)[nonneg.position != 0]
 nonneg.position.cons <- nonneg.position[nonneg.position != 0]
 nonneg.position.varb
 nonneg.position.cons
 hrep2[nonneg.position.cons, ]
 # where are they in hrep4?
 idx.ineq <- (1:nrow(hrep2))[! is.equality]
 idx.ineq
 i5 <- match(nonneg.position.cons, idx.ineq)
 i5
 amat[i5, ]
 # this says that variables
 nonneg.position.varb
 # correspond to NC constraints numbered
 i5
 # do we buy this ?????
 # have this many vertices
 nrow(foo)
 nrow(v5)
 # and this many non-redundant inequality constraints in NC
 nrow(amat)
 length(bvec)
 wtf <- matrix(NaN, nrow(foo), nrow(amat))
 for (i in 1:nrow(foo))
     for (j in 1:nrow(amat))
         wtf[i, j] <- qsign(qmq(bvec[j],
             qmatmult(amat[j, , drop = FALSE], t(v5[i, , drop = FALSE]))))
 all(wtf >= 0)
 # wtf
 wtf.too <- list()
 for (j in 1:nrow(amat)) {
     fred <- foo[wtf[ , j] == 0, , drop = FALSE]
     sally <- apply(fred == "0", 2, all)
     wtf.too[[j]] <- seq(along = sally)[sally]
 }
 alfred.e.newman <- wtf.too[i5]
 all(sapply(alfred.e.newman, length) == 1)
 all(unlist(alfred.e.newman) == nonneg.position.varb)


