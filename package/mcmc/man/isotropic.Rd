\name{isotropic}
\alias{isotropic}
\alias{isotropic.logjacobian}
\title{Isotropic Functions}
\description{
    Implements isotropic functions and the log-Jacobian of isotropic
    functions.
}
\usage{
isotropic(f)
isotropic.logjacobian(f, d.f)
}
\arguments{
  \item{f}{an R function that maps a non-negative real number to a
    non-negative real number.}
  \item{d.f}{an R function that computes the first derivative of
    \code{f} at the given point.}
}
\details{
  A function \eqn{h} is isotropic if it has the form
  \deqn{h(x) = \biggl\{ \begin{array}{l l}f(|x|) \frac{x}{|x|}, & x \neq 0
    \\ 0, & x = 0 \end{array}}{
    h(x) = if (|x| == 0)
             rep(0, length(x))
           else
	     f(|x|) * x / |x|}
  where \eqn{x} is a real vector, \eqn{|x|} is the usual Euclidean norm,
  and \eqn{f} maps \eqn{(0,\infty)}{(0, infinity)} to
  \eqn{(0,\infty)}{(0, infinity)}. The log-Jacobian of an isotropic
  function, \eqn{h} is given by
  \deqn{\log |\nabla h(x)|  = \biggl\{ \begin{array}{l
      l}\log\bigl(f'(|x|)\bigr) + (k-1) \log\bigl( f(|x|) / |x| \bigr),
         & x \neq 0 \\
    k \log\bigl(f'(|x|)\bigr), & x = 0\end{array}}{
    log | del h(x) | = if (|x| == 0)
                         log( f'(|x|) ) + (k-1)*log( f(|x|) / |x| )
                       else
 		         k * log( f'(|x|) )
  }
  where \eqn{k = length(x)}.

  A brief discussion of isotropic functions is contained in \cite{Johnson and
  Geyer (submitted)}.		 
}
\value{
  a function to calculate the corresponding isotropic function, or
  log-Jacobian of the isotropic function.}
\examples{
  identity <- function(x) x
  d.identity <- function(x) 1
  i.id <- isotropic(identity)
  lj.id <- isotropic.logjacobian(identity, d.identity)
  i.id(1:10) # This is an identity function.
  lj.id(1:10) # This is equal to 0.

  sq <- function(x) x^2
  d.sq <- function(x) 2 * x
  i.sq <- isotropic(sq)
  lj.sq <- isotropic.logjacobian(sq, d.sq)
  i.sq(c(1,1))
  lj.sq(c(1,1))
}
\references{
Johnson, L. T. and Geyer, C. J. (submitted)
Variable Transformation to Obtain Geometric Ergodicity
    in the Random-walk Metropolis Algorithm.
}
\seealso{
  \code{\link{metrop}}
  \code{\link{morph}}
}
\keyword{misc}