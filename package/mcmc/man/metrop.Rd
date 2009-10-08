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
      this object (up until version 0.7-2 this was incorrect with respect
      to the \code{debug} argument, now it applies to it too).}
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
Description of a complete MCMC analysis (Bayesian logistic regression)
using this function can be found in the vignette \code{demo}
(\url{../doc/demo.pdf}).
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
  Description of additional output when \code{debug = TRUE} can be
  found in the vignette \code{debug} (\url{../doc/debug.pdf}).
}
\section{Warning}{
If \code{outfun} is missing or not a function, then the log unnormalized
density can be defined without a \ldots argument and that works fine.
One can define it starting \code{ludfun <- function(state)} and that works
or \code{ludfun <- function(state, foo, bar)}, where \code{foo} and \code{bar}
are supplied as additional arguments to \code{metrop}.

If \code{outfun} is a function, then both it and the log unnormalized
density function can be defined without \ldots arguments \emph{if they
have exactly the same arguments list} and that works fine.  Otherwise it
doesn't work.  Start the definitions \code{ludfun <- function(state, foo)}
and \code{outfun <- function(state, bar)} and you get an error about
unused arguments.  Instead start the definitions
\code{ludfun <- function(state, foo, \ldots)}
and \code{outfun <- function(state, bar, \ldots)}, supply
\code{foo} and \code{bar} as additional arguments to \code{metrop},
and that works fine.

In short, the log unnormalized density function and \code{outfun} need
to have \ldots in their arguments list to be safe.  Sometimes it works
when \ldots is left out and sometimes it doesn't.

Of course, one can avoid this whole issue by always defining the log
unnormalized density function and \code{outfun} to have only one argument
\code{state} and use global variables (objects in the R global environment) to
specify any other information these functions need to use.  That too
follows the R way.  But some people consider that bad programming practice.
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
