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
\method{metrop}{function}(obj, initial, nbatch, blen = 1, nspac = 1,
    scale = 1, outfun, debug = FALSE, ...)
\method{metrop}{metropolis}(obj, initial, nbatch, blen = 1, nspac = 1,
    scale = 1, outfun, debug = FALSE, ...)
}
\arguments{
  \item{obj}{Either an \R function or an object of class \code{"metropolis"}
      from a previous invocation of this function.

      If a function, it evaluates the log unnormalized probability
      density of the desired equilibrium distribution of the Markov chain.
      Its first argument is the state vector of the Markov chain.  Other
      arguments arbitrary and taken from the \code{...} arguments of this
      function.
      It should return \code{-Inf} for points of the state space having
      probability zero under the desired equilibrium distribution.
      See also Details and Warning.

      If an object of class \code{"metropolis"}, any missing arguments
      (including the log unnormalized density function) are taken from
      this object.  Also \code{initial} is ignored and the initial state
      of the Markov chain is the final state from the run recorded in
      \code{obj}.
  }
  \item{initial}{a real vector, the initial state of the Markov chain.
      Must be feasible, see Details.  Ignored if \code{obj} is of
      class \code{"metropolis"}.}
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
Runs a \dQuote{random-walk} Metropolis algorithm, terminology introduced
by Tierney (1994), with multivariate normal proposal
producing a Markov chain with equilibrium distribution having a specified
unnormalized density.  Distribution must be continuous.  Support of the
distribution is the support of the density specified by argument \code{obj}.
The initial state must satisfy \code{obj(state, ...) > -Inf}.
Description of a complete MCMC analysis (Bayesian logistic regression)
using this function can be found in the vignette
\code{vignette("demo", "mcmc")}.

Suppose the function coded by the log unnormalized function (either
\code{obj} or \code{obj$lud}) is actually a log unnormalized density,
that is, if \eqn{w} denotes that function, then \eqn{e^w}{exp(w)} integrates
to some value strictly between zero and infinity.  Then the \code{metrop}
function always simulates a reversible, Harris ergodic Markov chain having
the equilibrium distribution with this log unnormalized density.
The chain is not guaranteed to be geometrically ergodic.  In fact it cannot
be geometrically ergodic if the tails of the log unnormalized density are
sufficiently heavy.  The \code{\link{morph.metrop}} function deals with this
situation.
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
  \item{accept.batch}{a vector of length \code{nbatch}, the batch means
      of the acceptances.}
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
doesn't work.  Define these functions by
\preformatted{
ludfun <- function(state, foo)
outfun <- function(state, bar)
}
and you get an error about unused arguments.  Instead define these functions by
\preformatted{
ludfun <- function(state, foo, \ldots)
outfun <- function(state, bar, \ldots)
}
and supply \code{foo} and \code{bar} as additional arguments to \code{metrop},
and that works fine.

In short, the log unnormalized density function and \code{outfun} need
to have \ldots in their arguments list to be safe.  Sometimes it works
when \ldots is left out and sometimes it doesn't.

Of course, one can avoid this whole issue by always defining the log
unnormalized density function and \code{outfun} to have only one argument
\code{state} and use global variables (objects in the \R global environment) to
specify any other information these functions need to use.  That too
follows the \R way.  But some people consider that bad programming practice.
}
\section{Philosophy of MCMC}{
This function follows the philosophy of MCMC explained
the introductory chapter of the
\emph{Handbook of Markov Chain Monte Carlo} (Geyer, 2011).

This function automatically does batch means in order to reduce
the size of output and to enable easy calculation of Monte Carlo standard
errors (MCSE), which measure error due to the Monte Carlo sampling (not
error due to statistical sampling --- MCSE gets smaller when you run the
computer longer, but statistical sampling variability only gets smaller
when you get a larger data set).  All of this is explained in the package
vignette \code{vignette("demo", "mcmc")} and in Section 1.10 of Geyer (2011).

This function does not apparently
do \dQuote{burn-in} because this concept does not actually help with MCMC
(Geyer, 2011, Section 1.11.4) but the re-entrant property of this
function does allow one to do \dQuote{burn-in} if one wants.
Assuming \code{ludfun}, \code{start.value}, \code{scale}
have been already defined
and are, respectively, an \R function coding the log unnormalized density
of the target distribution, a valid state of the Markov chain,
and a useful scale factor,
\preformatted{
out <- metrop(ludfun, start.value, nbatch = 1, blen = 1e5, scale = scale)
out <- metrop(out, nbatch = 100, blen = 1000)
}
throws away a run of 100 thousand iterations before doing another run of
100 thousand iterations that is actually useful for analysis, for example,
\preformatted{
apply(out$batch, 2, mean)
apply(out$batch, 2, sd) / sqrt(out$nbatch)
}
give estimates of posterior means and their MCSE assuming the batch length
(here 1000) was long enough to contain almost all of the signifcant
autocorrelation (see Geyer, 2011, Section 1.10, for more on MCSE).
The re-entrant property of this function (the second run starts
where the first one stops) assures that this is really \dQuote{burn-in}.

The re-entrant property allows one to do very long runs without having to
do them in one invocation of this function.
\preformatted{
out2 <- metrop(out)
out3 <- metrop(out2)
batch <- rbind(out$batch, out2$batch, out3$batch)
}
produces a result as if the first run had been three times as long.
}
\section{Tuning}{
The \code{scale} argument must be adjusted so that the acceptance rate
is not too low or too high to get reasonable performance.  The rule of
thumb is that the acceptance rate should be about 25\%.
But this recommendation (Gelman, et al., 1996) is justified by analysis
of a toy problem (simulating a spherical multivariate normal distribution)
for which MCMC is unnecessary.  There is no reason to believe this is optimal
for all problems (if it were optimal, a stronger theorem could be proved).
Nevertheless, it is clear that at very low acceptance rates the chain makes
little progress (because in most iterations it does not move) and that at
very high acceptance rates the chain also makes little progress (because
unless the log unnormalized density is nearly constant, very high acceptance
rates can only be achieved by very small values of \code{scale} so the
steps the chain takes are also very small.

Even in the Gelman, et al. (1996) result, the optimal rate for spherical
multivariate normal depends on dimension.  It is 44\% for \eqn{d = 1}
and 23\% for \eqn{d = \infty}{d = infinity}.
Geyer and Thompson (1995) have an example, admittedly for
simulated tempering (see \code{\link{temper}}) rather than random-walk
Metropolis, in which no acceptance rate less than 70\% produces an ergodic
Markov chain.  Thus 25\% is merely a rule of thumb.  We only know we don't
want too high or too low.  Probably 1\% or 99\% is very inefficient.
}
\references{
Gelman, A., Roberts, G. O., and Gilks, W. R. (1996)
Efficient Metropolis jumping rules.
In \emph{Bayesian Statistics 5: Proceedings of the Fifth Valencia
    International Meeting}.  Edited by J. M. Bernardo,
    J. O. Berger, A. P. Dawid, and A. F. M. Smith.
Oxford University Press, Oxford, pp. 599--607. 

Geyer, C. J. (2011)
Introduction to MCMC.
In \emph{Handbook of Markov Chain Monte Carlo}. Edited by S. P. Brooks,
A. E. Gelman, G. L. Jones, and X. L. Meng.
Chapman & Hall/CRC, Boca Raton, FL, pp. 3--48.

Geyer, C. J. and Thompson, E. A. (1995).
Annealing Markov chain Monte Carlo with applications to ancestral inference.
\emph{Journal of the American Statistical Association} \bold{90} 909--920.

Tierney, L. (1994)
Markov chains for exploring posterior distributions (with discussion).
\emph{Annals of Statistics} \bold{22} 1701--1762.
}
\seealso{
\code{\link{morph.metrop}} and \code{\link{temper}}
}
\examples{
h <- function(x) if (all(x >= 0) && sum(x) <= 1) return(1) else return(-Inf)
out <- metrop(h, rep(0, 5), 1000)
out$accept
# acceptance rate too low
out <- metrop(out, scale = 0.1)
out$accept
t.test(out$accept.batch)$conf.int
# acceptance rate o. k. (about 25 percent)
plot(out$batch[ , 1])
# but run length too short (few excursions from end to end of range)
out <- metrop(out, nbatch = 1e4)
out$accept
plot(out$batch[ , 1])
hist(out$batch[ , 1])
acf(out$batch[ , 1], lag.max = 250)
# looks like batch length of 250 is perhaps OK
out <- metrop(out, blen = 250, nbatch = 100)
apply(out$batch, 2, mean) # Monte Carlo estimates of means
apply(out$batch, 2, sd) / sqrt(out$nbatch) # Monte Carlo standard errors
t.test(out$accept.batch)$conf.int
acf(out$batch[ , 1]) # appears that blen is long enough
}
\keyword{misc}
