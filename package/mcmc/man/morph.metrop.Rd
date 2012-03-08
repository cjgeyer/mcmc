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
  outfun, debug = FALSE, morph, ...)
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
  \item{morph}{morph object used for transformations.  See \code{\link{morph}}.}
  \item{...}{see \code{\link{metrop}}.}
}
\details{
  \code{morph.metrop} implements morphemetric methods for Markov
  chains.  The caller specifies a log unnormalized probability density
  and a transformation.  The transformation specified by the
  \code{morph} parameter is used to induce a new log unnormalized
  probability   density, a Metropolis algorithm is
  run for the induced density.  The Markov chain is transformed back to
  the original scale.  Running the Metropolis algorithm for the induced
  density, instead of the original density, can result in a Markov chain
  with better convergence properties.  For more details see Johnson and Geyer
  (submitted).  Except for \code{morph}, all parameters are
  passed to \code{\link{metrop}}, transformed when necessary.  The
  \code{scale} parameter is \emph{not} transformed.

  If \eqn{X} is a real vector valued continuous random variable, and
  \eqn{Y = f(X)} where \eqn{f} is a diffeomorphism, then the pdf of
  \eqn{Y} is given by \deqn{f_Y(y) = f_X(f^{-1}(y)) | \nabla f^{-1}(y)
  |}{ fY(y) = fX(f^{-1}(y)) | del f^{-1}(y) |} where \eqn{f_X}{fX} is
  the pdf of \eqn{X} and \eqn{\nabla f^{-1}}{del f^{-1}} is the Jacobian
  of \eqn{f^{-1}}.  Because \eqn{f} is a diffeomorphism, a Markov chain
  for \eqn{f_Y}{fY} may be transformed into a Markov chain for
  \eqn{f_X}{fX}.  Furthermore, these Markov chains are isomorphic
  (Johnson and Geyer(submitted)) and have the same convergence rate.
  The \code{\link{morph}} variable provides a diffeomorphism,
  \code{morph.metrop} uses this diffeomorphism to induce the log
  unnormalized density, \eqn{\log f_Y}{log fY} based on the user
  supplied log unnormalized density, \eqn{\log f_X}{log fX}.
  \code{morph.metrop} runs a Metropolis algorithm for \eqn{\log f_Y}{log
  fY} and transforms the resulting Markov chain into a Markov chain for
  \eqn{f_X}{fX}.  The user accesible output components are the same as
  those that come from \code{\link{metrop}}, see the documentation for
  \code{\link{metrop}} for details.

  Subsequent calls of \code{morph.metrop} may change to the
  transformation by specifying a new value for \code{morph}.

  Any of the other parameters to \code{morph.metrop} may also be
  modified in subsequent calls.  See \code{\link{metrop}} for more details.
}
\value{
  an object of class \code{mcmc}, subclass \code{morph.metropolis}.
  This object is a list containing all of the elements from an object
  returned by \code{\link{metrop}}, plus at least the following
  components:
  \item{morph}{the morph object used for the transformations.}
  \item{morph.final}{the final state of the markov chain on the
    transformed scale.}
}
\examples{
out <- morph.metrop(function(x) dt(x, df=3, log=TRUE), 0, blen=100,
  nbatch=100, morph=morph(b=1))
# change the transformation.
out <- morph.metrop(out, morph=morph(b=2))
out$accept
# accept rate is high, increase the scale.
out <- morph.metrop(out, scale=4)
# close to 0.20 is about right.
out$accept
}
\references{
Johnson, L. T. and Geyer, C. J. (submitted)
Variable Transformation to Obtain Geometric Ergodicity
    in the Random-walk Metropolis Algorithm.
}
\seealso{
  \code{\link{metrop}}, \code{\link{morph}}.
}
\keyword{misc}
