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
  \item{p}{Real number strictly greater than 2.  May be missing.  If
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
  \eqn{Y} is given by
  \deqn{f_Y(y) = f_X(f^{-1}(y)+c) | \nabla f^{-1}(y)
    |}{ fY(y) = fX(f^{-1}(y)+c) | del f^{-1}(y) |}
  where \eqn{f_X}{fX} is the pdf of \eqn{X} and \eqn{\nabla f^{-1}}{del
  f^{-1}} is the Jacobian of \eqn{f^{-1}}.  Because \eqn{f} is a
  diffeomorphism, a Markov chain for \eqn{f_Y}{fY} may be transformed into a
  Markov chain for \eqn{f_X}{fX}.  Furthermore, these Markov chains are
  isomorphic (Johnson and Geyer(submitted)) and have the same
  convergence rate.  The parameters \code{r}, \code{p}, \code{b}, and
  \code{center} specify a diffemorphism (details below).
  \code{morph.metrop} uses this diffeomorphism to induce the log
  unnormalized density, \eqn{\log  f_Y}{log fY} based on the user
  supplied log unnormalized density, \eqn{\log f_X}{log fX}.
  \code{morph.metrop} runs a Metropolis algorithm for \eqn{\log f_Y}{log
  fY} and transforms the resulting Markov chain into a Markov chain for
  \eqn{f_X}{fX}.  The user useable output components are the same as those
  that come from \code{\link{metrop}}, see the documentation for
  \code{\link{metrop}} for details.

  The parameters \code{r}, \code{p}, \code{b} and \code{center} specify the
  transformation function.  In all cases, \code{center} gives the center
  of the transformation, which is the value \eqn{c} in the equation
  \deqn{Y = f(X - c).}
  The parameters \code{r}, \code{p} and \code{b} specify a function
  \eqn{g}, which is a monotonically increasing bijection from the
  non-negative reals to the non-negative reals.  Then
  \deqn{f(X) = g\bigl(|X|\bigr) \frac{X}{|X|}}{f(X) = g(|X|) * X / |X|}
  where \eqn{|X|} represents the Euclidean norm of the vector \eqn{X}.
  The inverse function is given by
  \deqn{f^{-1}(Y) = g^{-1}\bigl(|Y|\bigr) \frac{Y}{|Y|}.}{f^{-1}(Y) = g^{-1}(|Y|) * Y / |Y|.}

  The parameters \code{r} and \code{p} are used to define the function
  \deqn{g_1(x) = x + (x-r)^p I(x > r)}{g1(x) = x + (x-r)^p * I(x > r)}
  where \eqn{I( \cdot )}{I(GET UNICODE DOT)}  is the indicator
  function.  We require that \code{r} is non-negative and \code{p} is
  strictly greater than 2.  The parameter \code{b} is used to define the
  function
  \deqn{g_2(x) = \bigl(e^{bx} - e / 3\bigr) I(x > \frac{1}{b}) + 
    \bigl(x^3 b^3 e / 6 + x b e / 2\bigr) I(x \leq
    \frac{1}{b})}{
    g2(x) = (exp(b * x) - exp(1) / 3) * I(x > 1 / b) +
            (x^3 * b^3 exp(1) / 6 + x * b * exp(1) / 2) * I(x <= 1 / b).}
  We require that \eqn{b} is non-negative.

  \code{morph.metrop} uses the transformation specified by the inverse
  relationship, \deqn{X = f^{-1}(Y) + c}, where \eqn{c} is given by the
  value of \code{center}.  The parameters \code{r}, \code{p} and
  \code{b} specify \eqn{f^{-1}} in the following manner.
  
  If one of \code{r} and \code{p} is specified, and \code{b} is not
  specified, then 
  \deqn{f^{-1}(X) = g_1(|X|) \frac{X}{|X|}.}{f^{-1}(X) = g1(|X|) * X / |X|.}
  In the case of only specifying \code{r}, \code{p=3}.  If only \code{p}
  is specified, \code{r = 0} is used.
  
  If only \code{b} is specified, then 
  \deqn{f^{-1}(X) = g_2(|X|) \frac{X}{|X|}.}{f^{-1}(X) = g2(|X|) * X / |X|.}
  
  If one of \code{r} and \code{p} is specified, and \code{b} is also
  specified, then
  \deqn{f^{-1}(X) = g_2(g_1(|X|)) \frac{X}{|X|}.}{
    f^{-1}(X) = g2(g1(|X|)) * X / |X|.}

  Subsequent calls of \code{morph.metrop} may change to the
  transformation by specifying new values for \code{r}, \code{p},
  \code{b} or \code{center}.

  Any of the other parameters to \code{morph.metrop} may also be
  modified in subsequent calls.  See \code{\link{metrop}} for more details.
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
out <- morph.metrop(function(x) dt(x, df=3, log=TRUE), 0, blen=100,
  nbatch=100, b=1)
# change the transformation.
out <- morph.metrop(out, b=2)
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
  \code{\link{metrop}}
}
\keyword{misc}
