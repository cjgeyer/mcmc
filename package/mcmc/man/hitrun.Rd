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
hitrun(obj, nbatch, blen = 1, nspac = 1, a1, b1, a2, b2, outfun,
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
    a1 \%*\% x <= b1
    a2 \%*\% x == b2
}
Since the constraint set must be bounded \code{a1} and \code{b1} are not
allowed to be missing.  But \code{a2} and \code{b2} are allowed to be missing.
\code{a1}, \code{b1}, \code{a2}, and \code{b2} are ignored if \code{obj}
is an object of class \code{"hitrun"} (the constraint set and
log unnormalized density function are the same as for the previous run).

Suppose the function coded by the log unnormalized function is continuous
and nonzero on the constraint set.  Then the sampler is uniformly ergodic.
}
\value{
  an object of class \code{"mcmc"}, subclass \code{"hitrun"},
  which is a list containing at least the following components:
  \item{batch}{\code{nbatch} by \code{p} matrix, the batch means, where
      \code{p} is the dimension of the result of \code{outfun}
      if \code{outfun} is a function and the dimension
      of \code{state} when \code{outfun} is missing.}
  \item{initial}{initial state of Markov chain.}
  \item{final}{final state of Markov chain.}
  \item{initial.seed}{value of \code{.Random.seed} before the run.}
  \item{final.seed}{value of \code{.Random.seed} after the run.}
  \item{time}{running time of Markov chain from \code{system.time()}.}
  \item{ludfun}{the function used to calculate log unnormalized density,
  either \code{obj} or \code{obj$ludfun} from a previous run.}
  \item{nbatch}{the argument \code{nbatch} or \code{obj$nbatch}.}
  \item{blen}{the argument \code{blen} or \code{obj$blen}.}
  \item{nspac}{the argument \code{nspac} or \code{obj$nspac}.}
  \item{outfun}{the argument \code{outfun} or \code{obj$outfun}.}
%   Description of additional output when \code{debug = TRUE} can be
%   found in the vignette \code{debug} (\url{../doc/debug.pdf}).
}
\section{Warning}{
If \code{outfun} is missing, then the log unnormalized
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
# Bayesian inference for discrete probability distribution on {1, ..., d}
# state is probability vector p of length d
d <- 10
x <- 1:d
# equality constraints
#     mean equal to (d + 1) / 2, that is, sum(x * p) = (d + 1) / 2
#     probabilities sum to one, that is, sum(p) = 1
# inequality constraints
#     median less than or equal to (d + 1) / 2, that is,
#         sum(p[x <= (d + 1) / 2]) <= 1 / 2
#     probabilities nonnegative, that is, all(p >= 0)
a2 <- rbind(rep(1, d), x)
b2 <- c(1, (d + 1) / 2)
a1 <- - diag(d)
b1 <- rep(0, d)
a1 <- rbind(a1, as.numeric(x <= (d + 1) / 2))
b1 <- c(b1, 1 / 2)
# simulate prior, which Dirichlet(alpha)
# posterior would be another Dirichlet with n + alpha - 1,
#    where n is count of IID data for each value
alpha <- rep(2.3, d)
ludfun <- function(x) {
    stopifnot(is.numeric(x))
    stopifnot(is.finite(x))
    stopifnot(length(x) == d)
    stopifnot(is.numeric(alpha))
    stopifnot(is.finite(alpha))
    stopifnot(length(alpha) == d)
    stopifnot(alpha > 0)
    sum((alpha - 1) * log(x))
}
out <- hitrun(ludfun, nbatch = 30, blen = 250,
    a1 = a1, b1 = b1, a2 = a2, b2 = b2)
# prior means
round(colMeans(out$batch), 3)
# Monte Carlo standard errors
round(apply(out$batch, 2, sd) / sqrt(out$nbatch), 3)
}
\keyword{misc}
