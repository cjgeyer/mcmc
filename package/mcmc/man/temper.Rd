\name{temper}
\alias{temper}
\alias{temper.function}
\alias{temper.tempering}
\title{Simulated Tempering and Umbrella Sampling}
\description{
    Markov chain Monte Carlo for continuous random vectors using parallel
    or serial simulated tempering, also called umbrella sampling.  For
    serial tempering the state of the Markov chain is a pair \eqn{(i, x)},
    where \eqn{i} is an integer between 1 and \eqn{k} and \eqn{x} is a vector
    of length \eqn{p}.  This pair is represented as a single real vector
    \code{c(i, x)}.  For parallel tempering the state of the Markov chain
    is vector of vectors \eqn{(x_1, \ldots, x_k)}{(x[1], \ldots, x[k])},
    where each \code{x} is of length \eqn{p}.  This vector of vectors is
    represented as a \eqn{k \times p}{k by p} matrix.
}
\usage{
temper(obj, initial, neighbors, nbatch, blen = 1, nspac = 1, scale = 1,
    outfun, debug = FALSE, parallel = FALSE, \dots)
}
\arguments{
  \item{obj}{either an \R function or an object of class \code{"tempering"} from
    a previous run.  If a function, should evalutate the log unnormalized
    density \eqn{\log h(i, x)}{log h(i, x)} of the desired equilibrium
    distribution of the Markov chain for serial tempering (the same function
    is used for both serial and parallel tempering, see details below for
    further explanation).  If an object, the log unnormalized density function
    is \code{obj$lud}, and missing arguments of \code{temper} are
    obtained from the corresponding elements of \code{obj}.
    The first argument of the log unnormalized density function is the
    state for simulated tempering \eqn{(i, x)} is supplied as an \R vector
    \code{c(i, x)}; other arguments are arbitrary and taken from
    the \code{\dots} arguments of \code{temper}.  The log unnormalized density
    functon should return \code{-Inf} for points of the state space having
    probability zero.}
  \item{initial}{for serial tempering, a real vector \code{c(i, x)} as
    described above.  For parallel tempering, a real
    \eqn{k \times p}{k by p} matrix as described above.  In either case,
    the initial state of the Markov chain.}
  \item{neighbors}{a logical symmetric matrix of dimension \code{k}
    by \code{k}.  Elements that are \code{TRUE} indicate jumps
    or swaps attempted by the Markov chain.}
  \item{nbatch}{the number of batches.}
  \item{blen}{the length of batches.}
  \item{nspac}{the spacing of iterations that contribute to batches.}
  \item{scale}{controls the proposal step size for real elements of the state
    vector.  For serial tempering, proposing a new value for the \eqn{x}
    part of the state \eqn{(i, x)}.  For parallel tempering, proposing
    a new value for the \eqn{x_i}{x[i]} part of the state
    \eqn{(x_1, \ldots, x_k)}{(x[1], \ldots, x[k])}.  In either case, the proposal
    is a real vector of length \eqn{p}.  If scalar or vector, the proposal
    is \code{x + scale * z} where \code{x} is the part \eqn{x} or
    \eqn{x_i}{x[i]} of the state the proposal may replace.
    If matrix, the proposal is
    \code{x + scale \%*\% z}.  If list, the length must be \code{k},
    and each element must be scalar, vector, or matrix, and operate as
    described above.  The \eqn{i}-th component of the list is used to update
    \eqn{x} when the state is \eqn{(i, x)} or \eqn{x_i}{x[i]} otherwise.}
  \item{outfun}{controls the output.  If a function, then the batch means
      of \code{outfun(state, \dots)} are returned.  The argument \code{state}
      is like the argument \code{initial} of this function.  If missing, the
      batch means of the real part of the state vector or matrix are returned,
      and for serial tempering the batch means of a multivariate Bernoulli
      indicating the current component are returned.}
  \item{debug}{if \code{TRUE} extra output useful for testing.}
  \item{parallel}{if \code{TRUE} does parallel tempering, if \code{FALSE} does
      serial tempering.}
  \item{...}{additional arguments for \code{obj} or \code{outfun}.}
}
\details{
Serial tempering simulates a mixture of distributions of a continuous random
vector.  The number of components of the mixture is \code{k}, and the dimension
of the random vector is \code{p}.  Denote the state \eqn{(i, x)}, where \eqn{i}
is a positive integer between 1 and \eqn{k}, and let \eqn{h(i, x)} denote
the unnormalized joint density of their equilibrium distribution.
The logarithm of this function is what \code{obj} or \code{obj$lud} calculates.
The mixture distribution is the marginal for \eqn{x} derived from
the equilibrium distribution \eqn{h(i, x)}, that is,
\deqn{h(x) = \sum_{i = 1}^k h(i, x)}{h(x) = sum[i = 1 to k] h(i, x)}

Parallel tempering simulates a product of distributions of a continuous random
vector.  Denote the state \eqn{(x_1, \ldots, x_k)}{(x[1], \ldots, x[k])},
then the unnormalized joint density of the equilibrium distribution is
\deqn{h(x_1, \ldots, x_k) = \prod_{i = 1}^k h(i, x_i)}{h(x[1], \dots, x[k]) = prod[i = 1 to k] h(i, x[i])}

The update mechanism of the Markov chain combines two kinds of elementary
updates: jump/swap updates (jump for serial tempering, swap for parallel
tempering) and within-component updates.  Each iteration of the Markov chain
one of these elementary updates is done.  With probability 1/2 a jump/swap
update is done, and with probability 1/2 a with-component update is done.

Within-component updates are the same for both serial and parallel tempering.
They are \dQuote{random-walk} Metropolis updates with multivariate normal
proposal, the proposal distribution being determined by the argument
\code{scale}.  In serial tempering, the \eqn{x} part of the current state
\eqn{(i, x)} is updated preserving \eqn{h(i, x)}.
In parallel tempering, an index \eqn{i} is chosen at random and the part
of the state \eqn{x_i}{x[i]} representing that component is updated,
again preserving \eqn{h(i, x)}.

Jump updates choose uniformly at random a neighbor of the current component:
if \eqn{i} indexes the current component, then it chooses uniformly at random
a \eqn{j} such that \code{neighbors[i, j] == TRUE}.  It then does does a
Metropolis-Hastings update for changing the current state from \eqn{(i, x)}
to \eqn{(j, x)}.

Swap updates choose a component uniformly at random and a neighbor of that
component uniformly at random: first an index \eqn{i} is chosen uniformly
at random between 1 and \eqn{k}, then an index \eqn{j} is chosen
uniformly at random such that \code{neighbors[i, j] == TRUE}.  It then does
does a Metropolis-Hastings update for swapping the states of the
two components: interchanging \eqn{x_i}{x[i, ]} and \eqn{x_j}{x[j, ]}
while perserving \eqn{h(x_1, \ldots, x_k)}{h(x[1], \dots, x[k])}.

The initial state must satisfy \code{lud(initial, ...) > - Inf} for serial
tempering or must satisfy \code{lud(initial[i, ], ...) > - Inf} for each
\code{i} for parallel tempering, where \code{lud} is either \code{obj}
or \code{obj$lud}.
That is, the initial state must have positive probability.
}
\value{
  an object of class \code{"mcmc"}, subclass \code{"tempering"},
  which is a list containing at least the following components:
  \item{batch}{the batch means of the continuous part of the state.
    If \code{outfun} is missing, an \code{nbatch} by \code{k} by \code{p}
    array.  Otherwise, an \code{nbatch} by \code{m} matrix, where \code{m}
    is the length of the result of \code{outfun}.}
  \item{ibatch}{(returned for serial tempering only) an \code{nbatch}
    by \code{k} matrix giving batch means for the multivariate Bernoulli
    random vector that is all zeros except for a 1 in the \code{i}-th place
    when the current state is \eqn{(i, x)}.}
  \item{acceptx}{fraction of Metropolis within-component proposals accepted.
    A vector of length \code{k} giving the acceptance rate for each component.}
  \item{accepti}{fraction of Metropolis jump/swap proposals accepted.
    A \code{k} by \code{k} matrix giving the acceptance rate for each allowed
    jump or swap component.  \code{NA} for elements such that the corresponding
    elements of \code{neighbors} is \code{FALSE}.}
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
If \code{outfun} is missing, then the log unnormalized
density function can be defined without a \ldots argument and that works fine.
One can define it starting \code{ludfun <- function(state)} and that works
or \code{ludfun <- function(state, foo, bar)}, where \code{foo} and \code{bar}
are supplied as additional arguments to \code{temper} and that works too.

If \code{outfun} is a function, then both it and the log unnormalized
density function can be defined without \ldots arguments \emph{if they
have exactly the same arguments list} and that works fine.  Otherwise it
doesn't work.  Start the definitions \code{ludfun <- function(state, foo)}
and \code{outfun <- function(state, bar)} and you get an error about
unused arguments.  Instead start the definitions
\code{ludfun <- function(state, foo, \ldots)}
and \code{outfun <- function(state, bar, \ldots)}, supply
\code{foo} and \code{bar} as additional arguments to \code{temper},
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
d <- 9
witch.which <- c(0.1, 0.3, 0.5, 0.7, 1.0)
ncomp <- length(witch.which)

neighbors <- matrix(FALSE, ncomp, ncomp)
neighbors[row(neighbors) == col(neighbors) + 1] <- TRUE
neighbors[row(neighbors) == col(neighbors) - 1] <- TRUE

ludfun <- function(state, log.pseudo.prior = rep(0, ncomp)) {
    stopifnot(is.numeric(state))
    stopifnot(length(state) == d + 1)
    icomp <- state[1]
    stopifnot(icomp == as.integer(icomp))
    stopifnot(1 <= icomp && icomp <= ncomp)
    stopifnot(is.numeric(log.pseudo.prior))
    stopifnot(length(log.pseudo.prior) == ncomp)
    theta <- state[-1]
    if (any(theta > 1.0)) return(-Inf)
    bnd <- witch.which[icomp]
    lpp <- log.pseudo.prior[icomp]
    if (any(theta > bnd)) return(lpp)
    return(- d * log(bnd) + lpp)
}

# parallel tempering
thetas <- matrix(0.5, ncomp, d)
out <- temper(ludfun, initial = thetas, neighbors = neighbors, nbatch = 20,
    blen = 10, nspac = 5, scale = 0.56789, parallel = TRUE, debug = TRUE)

# serial tempering
theta.initial <- c(1, rep(0.5, d))
# log pseudo prior found by trial and error
qux <- c(0, 9.179, 13.73, 16.71, 20.56)

out <- temper(ludfun, initial = theta.initial, neighbors = neighbors,
    nbatch = 50, blen = 30, nspac = 2, scale = 0.56789,
    parallel = FALSE, debug = FALSE, log.pseudo.prior = qux)
}
\keyword{misc}
