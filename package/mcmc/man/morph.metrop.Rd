\name{morph.metrop}
\alias{morph.metrop}
\alias{morph.metrop.function}
\alias{morph.metrop.morph.metropolis}
\title{Morphemetric Metropolis Algorithm}
\description{
  Markov chain Monte Carlo for continuous random vector using a
  Metropolis algorithm for an induced density.
}
\usage{
morph.metrop(obj, initial, nbatch, blen = 1, nspac = 1, scale = 1,
  outfun, debug = FALSE, b, r, p, center = 0, ...)
}
\arguments{
  \item{obj}{see \code{\link{metrop}}.}
  \item{initial}{see \code{\link{metrop}}.}
  \item{nbatch}{see \code{\link{metrop}}.}
  \item{blen}{see \code{\link{metrop}}.}
  \item{nspac}{see \code{\link{metrop}}.}
  \item{scale}{see \code{\link{metrop}}.}
  \item{outfun}{see \code{\link{metrop}}.}
  \item{debug}{see \code{\link{metrop}}.}
  \item{r}{Non-negative real number.  May be missing.  If \code{p} is
    specified, \code{r} defaults to 0.}
  \item{p}{Real number strictly greater than 2.  May me missing.  If
    \code{r} is specified, \code{p} defaults to 3.}
  \item{b}{Positive real number.  May be missing.}
  \item{center}{Real scalar or vector.  If \code{center} is a vector it
    should be the same length of the state of the Markov chain.}
  \item{...}{see \code{\link{metrop}}.}
}
\details{
  \code{morph.metrop} implements morphemetric methods for Markov
  chains.  The caller specifies a log unnormalized probability density
  and a transformation.  The specified transformation is used to induce
  a new log unnormalized probability density, a Metropolis algorithm is
  run for the induced density.  The Markov chain is transformed back to
  the original scale.  Running the Metropolis algorithm for the induced
  density, instead of the original density, can result in a Markov chain
  with better convergence properties.  For more details see Johnson and Geyer
  (submitted).

  If \eqn{X} is a real vector valued continuous random variable, and
  \eqn{Y = f(X - c)} where \eqn{f} is a diffeomorphism, then the pdf of
  \eqn{Y} is given by \deqn{f_Y(y) = f_X(f^{-1}(y)+c) | \nabla f^{-1}(y)
  |}{ f_Y(y) = f_X(f^{-1}(y)+c) | del f^{-1}(y) |} where \eqn{f_X} is
  the pdf of \eqn{X}.  Because \eqn{f} is a diffeomorphism, a Markov
  chain for \eqn{f_Y} may be transformed into a Markov chain for
  \eqn{f_X}.  Furthermore, these Markov chains are isomorphic (Johnson
  and Geyer(submitted)) and have the same convergence rate.  The
  parameters \code{r}, \code{p}, \code{b}, and \code{center} specify a
  diffemorphism (details below).  \code{morph.metrop} uses this
  diffeomorphism to induce the log unnormalized density, \eqn{\log
  f_Y}{log f_Y} based on the user supplied log unnormalized density,
  \eqn{\log f_X}{log f_X}.  \code{morph.metrop} runs a Metropolis
  algorithm for \eqn{\log f_Y}{log f_Y} and transforms the resulting
  Markov chain into a Markov chain for \eqn{f_X}.  The user useable
  output components are the same as those that come from
  \code{\link{metrop}}, see the documentation for \code{\link{metrop}}
  for details.

  The parameters \code{r}, \code{p}, \code{b} and \code{center} specify the
  transformation function.  In all cases, \code{center} gives the center
  of the transformation, which is the value \eqn{c} in the equation
  \deqn{Y = f(X - c).}
  The parameters \code{r}, \code{p} and \code{b} specify a function
  \eqn{g}, which is a monotonically increasing bijection from the
  non-negative reals to the non-negative reals.  Then
  \deqn{f(X) = g(|X|) \frac{X}{|X|}}{f(X) = g(|X|) * X / |X|}
  where \eqn{|X|} represents the Euclidean norm of the vector \eqn{X}.
}
\value{
  an object of class \code{mcmc}, subclass \code{morph.metropolis}.
  This object is a list containing all of the elements from an object
  returned by \code{\link{metrop}}, plus at least the following
  components:
  \item{morph}{an internal object used for the transformations.}
  \item{morph.final}{the final state of the markov chain on the
    transformed scale.}
}
\examples{
}
\references{
Johnson, L. T. and Geyer, C. J. (submitted)
Variable Transformation to Obtain Geometric Ergodicity
    in the Random-walk Metropolis Algorithm.
}
\seealso{
  \code{\link{metrop}}
}
\keyword{misc}
