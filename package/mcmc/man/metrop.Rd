\name{metrop}
\alias{metrop}
\alias{metrop.function}
\alias{metrop.metropolis}
\title{Metropolis Algorithm}
\description{
    Markov chain Monte Carlo for continuous random vector using a Metropolis
    algorithm.
}
\usage{
metrop(obj, initial, nbatch, blen = 1, nspac = 1, scale = 1, outfun,
    debug = FALSE, ...)
}
\arguments{
  \item{obj}{an R function that evaluates the log unnormalized probability
      density of the desired equilibrium distribution of the Markov chain.
      First argument is the state vector of the Markov chain.  Other arguments
      arbitrary and taken from the \code{...} arguments of this function.
      Should return \code{- Inf} for points of the state space having
      probability zero under the desired equilibrium distribution.
      Alternatively, an object of class \code{"metropolis"} from a
      previous run can be supplied, in which case any missing arguments
      (including the log unnormalized density function) are taken from
      this object.}
  \item{initial}{a real vector, the initial state of the Markov chain.}
  \item{nbatch}{the number of batches.}
  \item{blen}{the length of batches.}
  \item{nspac}{the spacing of iterations that contribute to batches.}
  \item{scale}{controls the proposal step size.  If scalar or
          vector, the proposal is \code{x + scale * z} where \code{x} is
          the current state and \code{z} is a standard normal random vector.
          If matrix, the proposal is \code{x + scale \%*\% z}.}
  \item{outfun}{controls the output.  If a function, then the batch means
          of \code{outfun(state, ...)} are returned.  If a numeric
          or logical vector, then the batch means of \code{state[outfun]}
          (if this makes sense).  If missing, the the batch means
          of \code{state} are returned.}
  \item{debug}{if \code{TRUE} extra output useful for testing.}
  \item{...}{additional arguments for \code{obj} or \code{outfun}.}
}
\details{
Runs a \dQuote{random-walk} Metropolis algorithm with multivariate normal
proposal
producing a Markov chain with equilibrium distribution having a specified
unnormalized density.  Distribution must be continuous.  Support of the
distribution is the support of the density specified by argument \code{obj}.
The initial state must satisfy \code{obj(state, ...) > 0}.
}
\value{
  an object of class \code{"mcmc"}, subclass \code{"metropolis"},
  which is a list containing at least the following components:
  \item{accept}{fraction of Metropolis proposals accepted.}
  \item{batch}{\code{nbatch} by \code{p} matrix, the batch means, where
      \code{p} is the dimension of the result of \code{outfun}
      if \code{outfun} is a function, otherwise the dimension of
      \code{state[outfun]} if that makes sense, and the dimension
      of \code{state} when \code{outfun} is missing.}
  \item{initial}{value of argument \code{initial}.}
  \item{final}{final state of Markov chain.}
  \item{initial.seed}{value of \code{.Random.seed} before the run.}
  \item{final.seed}{value of \code{.Random.seed} after the run.}
  \item{time}{running time of Markov chain from \code{system.time()}.}
  \item{lud}{the function used to calculate log unnormalized density,
  either \code{obj} or \code{obj$lud} from a previous run.}
  \item{nbatch}{the argument \code{nbatch} or \code{obj$nbatch}.}
  \item{blen}{the argument \code{blen} or \code{obj$blen}.}
  \item{nspac}{the argument \code{nspac} or \code{obj$nspac}.}
  \item{outfun}{the argument \code{outfun} or \code{obj$outfun}.}
}
\examples{
h <- function(x) if (all(x >= 0) && sum(x) <= 1) return(1) else return(-Inf)
out <- metrop(h, rep(0, 5), 1000)
out$accept
# acceptance rate too low
out <- metrop(out, scale = 0.1)
out$accept
# acceptance rate o. k. (about 25 percent)
plot(out$batch[ , 1])
# but run length too short (few excursions from end to end of range)
out <- metrop(out, nbatch = 1e4)
out$accept
plot(out$batch[ , 1])
hist(out$batch[ , 1])
}
\keyword{misc}
