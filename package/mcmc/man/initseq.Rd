\name{initseq}
\alias{initseq}
\title{Initial Sequence Estimators}
\description{
    Variance of sample mean of functional of reversible Markov chain
    using methods of Geyer (1992).
}
\usage{
initseq(x)
}
\arguments{
  \item{x}{a numeric vector that is a scalar-valued functional of a reversible
      Markov chain.}
}
\details{
Let
\deqn{\gamma_k = \textrm{cov}(X_i, X_{i + k})}{gamma[k] = cov(x[i], x[i + k])}
considered as a function of the lag \eqn{k} be
the autocovariance function of the input time series.
Define
\deqn{\Gamma_k = \gamma_{2 k} + \gamma_{2 k + 1}}{Gamma[k] = gamma[2 k] + gamma[2 k + 1]}
the sum of consecutive pairs of autocovariances.  Then Theorem 3.1 in
Geyer (1992) says that \eqn{\Gamma_k}{Gamma[k]} considered as a function of
\eqn{k} is strictly positive, strictly decreasing, and strictly convex,
assuming the input time series is a scalar-valued functional of a reversible Markov
chain.  All of the MCMC done by this package is reversible.
This \R function estimates the \dQuote{big gamma} function,
\eqn{\Gamma_k}{Gamma[k]} considered as a function of
\eqn{k}, subject to three different constraints, (1) nonnegative,
(2) nonnegative and nonincreasing, and (3) nonnegative, nonincreasing,
and convex.  It also estimates the variance in the Markov chain central
limit theorem (CLT)
\deqn{\gamma_0 + 2 \sum_{k = 1}^\infty \gamma_k = - \gamma_0 + 2 \sum_{k = 0}^\infty \Gamma_k}{- gamma0 + 2 * sum(gamma) = - gamma0 + 2 * sum(Gamma)}

\strong{Note:} The batch means provided by \code{\link{metrop}} are also
scalar functionals of a reversible Markov chain.  Thus these initial sequence
estimators applied to the batch means give valid standard errors for the
mean of the match means even when the batch length is too short to provide
a valid estimate of asymptotic variance.  One does, of course, have to
multiply the asymptotic variance of the batch means by the batch length
to get the asymptotic variance for the unbatched chain.
}
\value{
a list containing the following components:

  \item{gamma0}{the scalar \eqn{\gamma_0}{gamma[0]}, the marginal variance
  of \code{x}.}
  \item{Gamma.pos}{the vector \eqn{\Gamma}{Gamma}, estimated so as to be nonnegative,
  where, as always, \R uses one-origin indexing so \code{Gamma.pos[1]} is 
  \eqn{\Gamma_0}{Gamma[0]}.}
  \item{Gamma.dec}{the vector \eqn{\Gamma}{Gamma}, estimated so as to be nonnegative
  and nonincreasing, where, as always,
  \R uses one-origin indexing so \code{Gamma.dec[1]} is 
  \eqn{\Gamma_0}{Gamma[0]}.}
  \item{Gamma.con}{the vector \eqn{\Gamma}{Gamma}, estimated so as to be nonnegative
  and nonincreasing and convex, where, as always,
  \R uses one-origin indexing so \code{Gamma.con[1]} is 
  \eqn{\Gamma_0}{Gamma[0]}.}
  \item{var.pos}{the scalar \code{- gamma0 + 2 * sum(Gamma.pos)}, which is
  the asymptotic variance in the Markov chain CLT.  Divide by \code{length(x)}
  to get the approximate variance of the sample mean of \code{x}.}
  \item{var.dec}{the scalar \code{- gamma0 + 2 * sum(Gamma.dec)}, which is
  the asymptotic variance in the Markov chain CLT.  Divide by \code{length(x)}
  to get the approximate variance of the sample mean of \code{x}.}
  \item{var.con}{the scalar \code{- gamma0 + 2 * sum(Gamma.con)}, which is
  the asymptotic variance in the Markov chain CLT.  Divide by \code{length(x)}
  to get the approximate variance of the sample mean of \code{x}.}
}
\section{Bugs}{
Not precisely a bug, but \code{var.pos}, \code{var.dec}, and \code{var.con}
can be negative.  This happens only when the chain is way too short to estimate
the variance, and even then rarely.  But it does happen.
}
\references{
Geyer, C. J. (1992)
Practical Markov Chain Monte Carlo.
\emph{Statistical Science} \bold{7} 473--483.
}
\seealso{
\code{\link{metrop}}
}
\examples{
n <- 2e4
rho <- 0.99
x <- arima.sim(model = list(ar = rho), n = n)
out <- initseq(x)
\dontrun{
plot(seq(along = out$Gamma.pos) - 1, out$Gamma.pos,
   xlab = "k", ylab = expression(Gamma[k]), type = "l")
lines(seq(along = out$Gamma.dec) - 1, out$Gamma.dec, col = "red")
lines(seq(along = out$Gamma.con) - 1, out$Gamma.con, col = "blue")
}
# asymptotic 95\% confidence interval for mean of x
mean(x) + c(-1, 1) * qnorm(0.975) * sqrt(out$var.con / length(x))
# estimated asymptotic variance
out$var.con
# theoretical asymptotic variance
(1 + rho) / (1 - rho) * 1 / (1 - rho^2)
# illustrating use with batch means
bm <- apply(matrix(x, nrow = 5), 2, mean)
initseq(bm)$var.con * 5
}
\keyword{ts}
