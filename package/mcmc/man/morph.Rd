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
      r=NULL, p=NULL, b=NULL,
      center=0)
}
\arguments{
  \item{f}{Transformation funtion.}
  \item{f.inv}{Inverse transformation function.}
  \item{logjacobian}{Function to compute the log-jacobian of the
    transformation.}
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
}
\value{
  a list of functions used by \code{metrop}.
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
