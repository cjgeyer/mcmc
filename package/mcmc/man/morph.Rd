\name{morph}
\encoding{UTF-8}
\alias{morph}
\alias{morph.identity}
\title{Variable Transformation}
\description{
  Utility functions for variable transformation.
}
\usage{
morph(b, r, p, center)
morph.identity()
}
\arguments{
  \item{b}{Positive real number.  May be missing.}
  \item{r}{Non-negative real number.  May be missing.  If \code{p} is
    specified, \code{r} defaults to 0.}
  \item{p}{Real number strictly greater than 2.  May be missing.  If
    \code{r} is specified, \code{p} defaults to 3.}
  \item{center}{Real scalar or vector.  May be missing.  If
    \code{center} is a vector it should be the same length of the state
    of the Markov chain, \code{center} defaults to 0}
}
\section{Warning}{
  The equations for the returned \code{transform} function (see below)
  do not have a general analytical solution when \code{p} is not equal
  to 3.  This implementation uses numerical approximation to calculate
  \code{transform} when \code{p} is not equal to 3.  If computation
  speed is a factor, it is advisable to use \code{p=3}.  This is not a
  factor when using \code{\link{morph.metrop}}, as \code{transform} is
  only called once during setup, and not at all while running the Markov chain.
}
\details{
  The \code{morph} function facilitates using variable transformations
  by providing functions to (using \eqn{X} for the original random
  variable with the pdf \eqn{f_X}{f.X}, and \eqn{Y} for the transformed
  random variable with the pdf \eqn{f_Y}{f.Y}):
  \itemize{
    \item Calculate the log unnormalized probability density for \eqn{Y}
    induced by the transformation.
    \item Transform an arbitrary function of \eqn{X} to a function of
    \eqn{Y}.
    \item Transform values of \eqn{X} to values of \eqn{Y}.
    \item Transform values of \eqn{Y} to values of \eqn{X}
       (the inverse transformation).
  }
  for a select few transformations.

  \code{morph.identity} implements the identity transformation,
  \eqn{Y=X}.
  
  The parameters \code{r}, \code{p}, \code{b} and \code{center} specify the
  transformation function.  In all cases, \code{center} gives the center
  of the transformation, which is the value \eqn{c} in the equation
  \deqn{Y = f(X - c).}  If no parameters are specified, the identity
  transformation, \eqn{Y=X}, is used.
  
  The parameters \code{r}, \code{p} and \code{b} specify a function
  \eqn{g}, which is a monotonically increasing bijection from the
  non-negative reals to the non-negative reals.  Then
  \deqn{f(X) = g\bigl(|X|\bigr) \frac{X}{|X|}}{f(X) = g(|X|) * X / |X|}
  where \eqn{|X|} represents the Euclidean norm of the vector \eqn{X}.
  The inverse function is given by
  \deqn{f^{-1}(Y) = g^{-1}\bigl(|Y|\bigr) \frac{Y}{|Y|}.}{f^{-1}(Y) = g^{-1}(|Y|) * Y / |Y|.}

  The parameters \code{r} and \code{p} are used to define the function
  \deqn{g_1(x) = x + (x-r)^p I(x > r)}{g1(x) = x + (x-r)^p * I(x > r)}
  where \eqn{I( \cdot )}{I(•)}  is the indicator
  function.  We require that \code{r} is non-negative and \code{p} is
  strictly greater than 2.  The parameter \code{b} is used to define the
  function
  \deqn{g_2(x) = \bigl(e^{bx} - e / 3\bigr) I(x > \frac{1}{b}) + 
    \bigl(x^3 b^3 e / 6 + x b e / 2\bigr) I(x \leq
    \frac{1}{b})}{
    g2(x) = (exp(b * x) - exp(1) / 3) * I(x > 1 / b) +
            (x^3 * b^3 exp(1) / 6 + x * b * exp(1) / 2) * I(x <= 1 / b).}
  We require that \eqn{b} is positive.

  The parameters \code{r}, \code{p} and \code{b} specify \eqn{f^{-1}} in
  the following manner:
  \itemize{
    \item  If one or both of \code{r} and \code{p} is specified, and \code{b}
    is not specified, then \deqn{f^{-1}(X) = g_1(|X|)
      \frac{X}{|X|}.}{f^{-1}(X) = g1(|X|) * X / |X|.}  If only
    \code{r} is specified, \code{p = 3} is used.  If only \code{p} is specified,
    \code{r = 0} is used.

    \item If only \code{b} is specified, then \deqn{f^{-1}(X) = g_2(|X|)
      \frac{X}{|X|}.}{f^{-1}(X) = g2(|X|) * X / |X|.}

    \item If one or both of \code{r} and \code{p} is specified, and \code{b} is
    also specified, then \deqn{f^{-1}(X) = g_2(g_1(|X|))
      \frac{X}{|X|}.}{f^{-1}(X) = g2(g1(|X|)) * X / |X|.}
  }
}
\value{
  a list containing the functions
  \itemize{
  \item \code{outfun(f)}, a function that operates on functions.
  \code{outfun(f)} returns the function \code{function(state, ...)
    f(inverse(state), ...)}.
  \item \code{inverse}, the inverse transformation function.
  \item \code{transform}, the transformation function.
  \item \code{lud}, a function that operates on functions.  As input,
  \code{lud} takes a function that calculates a log unnormalized
  probability density, and returns a function that calculates the
  log unnormalized density by transforming a random variable using the
  \code{transform} function.  \code{lud(f) = function(state, ...)
  f(inverse(state), ...) + log.jacobian(state)}, where
  \code{log.jacobian} represents the function that calculate the log
  Jacobian of the transformation.  \code{log.jacobian} is not returned.
}
}
\examples{

# use an exponential transformation, centered at 100.
b1 <- morph(b=1, center=100)
# original log unnormalized density is from a t distribution with 3
# degrees of freedom, centered at 100.
lud.transformed <- b1$lud(function(x) dt(x - 100, df=3, log=TRUE))
d.transformed <- Vectorize(function(x) exp(lud.transformed(x)))
\dontrun{
curve(d.transformed, from=-3, to=3, ylab="Induced Density")
}
}
\seealso{
  \code{\link{morph.metrop}}
}
\keyword{misc}
