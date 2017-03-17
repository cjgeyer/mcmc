\name{olbm}
\alias{olbm}
\title{Overlapping Batch Means}
\description{
    Variance of sample mean of time series calculated using overlapping
    batch means.
}
\usage{
olbm(x, batch.length, demean = TRUE)
}
\arguments{
  \item{x}{a matrix or time series object.  Each column of \code{x} is
      treated as a scalar time series.}
  \item{batch.length}{length of batches.}
  \item{demean}{when \code{demean = TRUE} (the default) the sample mean
      is subtracted from each batch mean when estimating the variance.
      Using \code{demean = FALSE} would essentially assume the true mean
      is known to be zero, which might be useful in a toy problem where
      the answer is known.}
}
\value{
  The estimated variance of the sample mean.
}
\seealso{
  \code{\link{ts}}
}
\examples{
h <- function(x) if (all(x >= 0) && sum(x) <= 1) return(1) else return(-Inf)
out <- metrop(h, rep(0, 5), 1000)
out <- metrop(out, scale = 0.1)
out <- metrop(out, nbatch = 1e4)
foo <- olbm(out$batch, 150)
# monte carlo estimates (true means are same by symmetry)
apply(out$batch, 2, mean)
# monte carlo standard errors (true s. d. are same by symmetry)
sqrt(diag(foo))
# check that batch length is reasonable
acf(out$batch, lag.max = 200)
}
\keyword{ts}
