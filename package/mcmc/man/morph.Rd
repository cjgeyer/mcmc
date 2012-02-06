\name{morph}
\alias{morph}
\alias{morph.identity}
\title{Change of Variable}
\description{
  Specify change-of-variable and inverse change-of-variable used
  to make metropolis algorithm geometrically ergodic.
  Used to provide argument for \code{\link{metrop}} function.
}
\usage{
morph.identity()
morph(f=NULL, f.inv=NULL, logjacobian=NULL,
      scale.fun=NULL, scale.inv=NULL,
      r=NULL, p=NULL, b=NULL,
      center=0)
}
\arguments{
  \item{f}{Transformation funtion.}
  \item{f.inv}{Inverse transformation function.}
  \item{logjacobian}{Function to compute the log-jacobian of the
    transformation.}
  \item{scale.fun}{Function to transform the \eqn{X} proposal scale to the
    \eqn{Y} proposal scale.  If \code{scale.fun} is missing, {f} is used.}
  \item{scale.inv}{Inverse function of \code{scale.inv}.  If
    \code{scale.inv} is missing, \code{f.inv} is used.}
  \item{r}{Non-negative real number.  May be missing.  If \code{p} is
    specified, \code{r} defaults to 0.}
  \item{p}{Real number strictly greater than 2.  May me missing.  If
    \code{r} is specified, \code{p} defaults to 3.}
  \item{b}{Positive real number.  May be missing.}
  \item{center}{Real scalar or vector.  If \code{center} is a vector it
    should be the same length of the state of the Markov chain.}
}
\details{
  The \code{morph} functions implement variable transformation for
  continuous random variables.  If \eqn{X} is a real vector valued
  continuous random variable, and \eqn{Y = f(X - c)} then the pdf
  of \eqn{Y} is given by
  \deqn{f_Y(y) = f_X(f^{-1}(y)+c) | \nabla f^{-1}(y) |}{
    f_Y(y) = f_X(f^{-1}(y)+c) | del f^{-1}(y) |}
  where \eqn{f_X} is the pdf of \eqn{X}.

  \code{morph.identity} implements the identity transformation,
  \eqn{Y=X}.

  In \code{morph}, the arguments \code{f}, \code{f.inv} and
  \code{logjacobian} are the functions \eqn{f}, \eqn{f^{-1}} and
  \eqn{y \mapsto \log | \nabla f^{-1}(y) |}{function(y) \log | del
    f^{-1}(y) |} respectively.  \code{center} gives the value \eqn{c}.

  \code{morph} is used by \code{\link{metrop}} to run isomorphic Markov
  chains as in Johnson and Geyer (Submitted).  Let \code{lud.x} be the
  log-unnormalized density of the random variable \code{X}, and \code{Y
  = f(X)} where \code{f} is a diffeomorphism with inverse \code{f.inv}. 
  Hence, \code{f.inv} (or \code{f}) will transform a Markov chain for
  the distribution of \code{X} (or \code{Y}) into an isomorphic Markov
  chain for the distribution of \code{Y} (or \code{X}).  By default,
  \code{\link{metrop}} runs a random-walk Metropolis algorithm for the
  distribution of \code{X}.  If \code{morph} is used, the
  \code{\link{metrop}} runs a random-walk Metropolis algorithm for the
  distribution of \code{Y}, but presents the isomoprhic Markov chain for
  the distribution of \code{X} to the caller.  As demonstrated in
  Johnson and Geyer (Submitted), a random-walk Metropolis algorithm for
  the distribution of \code{Y} may be geometrically ergodic when one for
  the corresponding distribution of \code{X} can not be, e.g. when
  \code{X} has a \eqn{t}-distribution.

  When using the \code{morph} parameter of \code{\link{metrop}}, the
  other parameters are kept on the '\code{X}' scale.  Specifically: the
  Markov chain will be started at \code{X[0] = initial} (\code{Y[0] =
  f(initial - center)}); the batch means of \code{outfun(X[i]) =
  outfun(f.inv(Y[i]) + center)} will be returned; and the proposal
  density for the random-walk Metropolis algorithm for the distribution
  of \code{Y} will have scale \code{scale.fun(scale)}.

  The functions \code{scale.fun} and \code{scale.inv} are used to
  transform and invert the \code{scale} (see \code{\link{metrop}})
  respectively.  If \code{scale.fun} is missing, \code{f} is used; if
  \code{scale.inv} is missing \code{f.inv} is used.  This is appropriate
  in the case of spherical transformation functions centered at zero
  (\code{\link{isotropic}} functions are spherical), but possibly
  inappropriate if the transformation function is non-spherical, and
  certainly inappropriate if a non-zero center is used.  In these cases,
  the user should provide \code{scale.fun} and \code{scale.inv}.  Failing to
  do so may cause a negative scale to be used, resulting in an error, or
  other undesired results.
  
}
\value{
  a list of functions;
  \item{outfun}{A function that operates on functions.  For a given
  function, \code{outX}, \code{outfun(outX)(y) = outX(f.inv(y) +
    center)}.}
  \item{transform}{A function that transforms \code{X} into \code{Y}.
    \code{transform(x) = f(x-center) = y}.}
  \item{inverse}{A function that transforms \code{Y} into \code{X}.
    \code{inverse{y} = f.inv(y) + center = x}.}
  \item{lud}{A function that calculates the log-unnormalized density of
  \code{Y}, \code{lud.y} from the log-unnormalized density of \code{X}, \code{lud.x}.
    \code{lud(lud.x)(y) = lud.y(y)}.}
  \item{f}{The uncentered transformation function.}
  \item{f.inv}{The inverse of \code{f}.}
  \item{logjacobian}{The log-jacobian function of \code{f.inv}.}
  \item{scale.fun}{Function to transform the \eqn{X} proposal scale to
    the \eqn{Y} proposal scale.}
  \item{scale.inv}{The inverse of \code{scale}.}
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
  \code{\link{isotropic}}
}
\keyword{misc}
