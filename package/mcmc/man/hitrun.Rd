\name{hitrun}
\alias{hitrun}
\alias{hitrun.function}
\alias{hitrun.hitrun}
\title{Hit and Run Algorithm on Convex Polytope}
\description{
    Markov chain Monte Carlo for continuous random vector on a convex polytope
    using a hit and run algorithm.
}
\usage{
hitrun(obj, initial, nbatch, blen = 1, nspac = 1, a1, b1, a2, b2, outfun,
    debug = FALSE, ...)
}
\arguments{
  \item{obj}{an R function that evaluates the log unnormalized probability
      density of the desired equilibrium distribution of the Markov chain.
      First argument is the state vector of the Markov chain.  Other arguments
      arbitrary and taken from the \code{...} arguments of this function.
      Should return \code{- Inf} for points of the state space having
      probability zero under the desired equilibrium distribution.
      Alternatively, an object of class \code{"hitrun"} from a
      previous run can be supplied, in which case any missing arguments
      (including the log unnormalized density function) are taken from
      this object.}
  \item{initial}{a real vector, the initial state of the Markov chain.}
  \item{nbatch}{the number of batches.}
  \item{blen}{the length of batches.}
  \item{nspac}{the spacing of iterations that contribute to batches.}
  \item{a1}{a numeric matrix. See details.}
  \item{b1}{a numeric vector. See details.}
  \item{a2}{a numeric matrix. See details.}
  \item{b2}{a numeric vector. See details.}
  \item{outfun}{controls the output.  If a function, then the batch means
          of \code{outfun(state, ...)} are returned.  If a numeric
          or logical vector, then the batch means of \code{state[outfun]}
          (if this makes sense).  If missing, the the batch means
          of \code{state} are returned.}
  \item{debug}{if \code{TRUE} extra output useful for testing.}
  \item{...}{additional arguments for \code{obj} or \code{outfun}.}
}
\details{
Runs a hit and run algorithm (for which see the references)
producing a Markov chain with equilibrium distribution having a specified
unnormalized density.  Distribution must be continuous.  Support of the
distribution is the support of the density specified by argument \code{obj}.
And by a polyhedral constraint set specified by the arguments
\code{a1}, \code{b1}, \code{a2}, and \code{b2}.  It is the set of all
\code{x} satisfying
\preformatted{
    a1 %*% x <= b1
    a2 %*% x == b2
}
Since the constraint set must be bounded \code{a1} and \code{b2} are not
allowed to be missing.  But \code{a2} and \code{b2} are allowed to be missing.

The initial state must be in the support of the equilibrium distribution.
That is, it must be in the constraint set and also must have a finite
value of the log unnormalized density function, that is,
\code{obj(initial, ...)} is finite.  Normally, the log unnormalized density
function is written that it never returns infinite values, so the constraints
define the support.

Suppose the function coded by the log unnormalized function is continuous
and nonzero on the constraint set.  Then the sampler is uniformly ergodic.
}
\value{
  an object of class \code{"mcmc"}, subclass \code{"hitrun"},
  which is a list containing at least the following components:
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
are supplied as additional arguments to \code{hitrun}.

If \code{outfun} is a function, then both it and the log unnormalized
density function can be defined without \ldots arguments \emph{if they
have exactly the same arguments list} and that works fine.  Otherwise it
doesn't work.  Start the definitions \code{ludfun <- function(state, foo)}
and \code{outfun <- function(state, bar)} and you get an error about
unused arguments.  Instead start the definitions
\code{ludfun <- function(state, foo, \ldots)}
and \code{outfun <- function(state, bar, \ldots)}, supply
\code{foo} and \code{bar} as additional arguments to \code{hitrun},
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
\references{
Belisle, C. J. P., Romeijn, H. E. and Smith, R. L. (1993)
Hit-and-run algorithms for generating multivariate distributions.
\emph{Mathematics of Operations Research}, \bold{18}, 255--266.

Chen, M. H. and Schmeiser, B. (1993)
Performance of the Gibbs, hit-and-run, and Metropolis samplers.
\emph{Journal of Computational and Graphical Statistics}, \bold{2}, 251--272.

}
\examples{
2 + 2
}
\keyword{misc}
