\name{morph}
\alias{morph}
\title{Change of Variable}
\description{
    Specify change-of-variable and inverse change-of-variable used
    to make metropolis algorithm geometrically ergodic.
    Used to provide argument for \code{\link{metrop}} function.
}
\usage{
morph(type = c("none", "exponential", "subexponential"), R, p, b)
}
\arguments{
  \item{type}{If "none", no change-of-variable is used.  See details section
    for other types.}
  \item{R}{Positive real number.  Ignored and may be missing
    if \code{type == "none"}.}
  \item{p}{Real number strictly greater than 2.  Ignored and may be missing
    if \code{type == "none"}.}
  \item{b}{Real number strictly greater than 0.  Ignored and may be missing
    if \code{type == "none"} or \code{type == "exponential"}.}
}
\details{
Define
\deqn{f(x) = \left\{\begin{array}{ll} x, & x < R \\ x + (x - R)^p, & x \geq R \end{array}\right.}{f(x) = x + (x - R)^p I(x > R)}
where \eqn{R \ge 0}{R >= 0} and \eqn{p > 2}\ifelse{latex}{}{and where I is the indicator function}.
Also define
\deqn{h(\gamma) = \left\{\begin{array}{ll} f(|\gamma|) \frac{\gamma}{|\gamma|}, & \gamma \neq 0 \\ 0, & \gamma = 0 \end{array}\right.}{h(\gamma) = f(|\gamma|) \gamma / |\gamma| I(\gamma != 0)}
whose inverse function is defined the same way except with \eqn{f}
replaced by \eqn{f^{-1}}{f^(-1)}.

In case \code{type == "exponential"} if the unnormalized density specified
by the user is for the random vector \eqn{\beta}, then the random-walk
Metropolis algorithm is actually done on the corresponding unnormalized density
for the random vector \eqn{\gamma} related to \eqn{\beta} by
\eqn{\beta = h(\gamma)}.  But the output is adjusted to be for the original
random vector \eqn{\beta} rather than for the transformed random vector.

Also define
\deqn{g(x) = \left\{\begin{array}{ll} e^{b x} - e / 3, & x > 1 / b \\ x^3 b^3 e / 6 + x b e / 2 & x \leq 1 / b \end{array}\right.}{g(x) = [e^{b x} - e / 3] I(x > 1 / b) + [x^3 b^3 e / 6 + x b e / 2] I(x <= 1 / b)}
where \eqn{b > 0}.

In case \code{type == "subexponential"} if the unnormalized density specified
by the user is for the random vector \eqn{\beta}, then the random-walk
Metropolis algorithm is actually done on the corresponding unnormalized density
for the random vector \eqn{\gamma} related to \eqn{\beta} by
\eqn{\beta = h(\gamma)}, where \eqn{h} is given by the equation above with
\eqn{f} replaced by the composition \eqn{g} and \eqn{f}.  As in the other
case the output is adjusted to be for the original
random vector \eqn{\beta} rather than for the transformed random vector.

If the log unnormalized density has super-exponential tail behavior, as when
normal priors are used for all parameters, then, typically, no transformation
need be used to obtain geometric ergodicity.

If the log unnormalized density has exponential tail behavior, as when
conjugate priors are used for an exponential family model, such as logistic
or Poisson regression or log-linear models for categorical data analysis,
then, typically (always for exponential families having natural statistic
that is bounded in any direction), the transformation of
\code{type == "exponential"} must be used to obtain geometric ergodicity
(with any valid \eqn{R} and \eqn{p}).

If the log unnormalized density has sub-exponential tail behavior, then,
typically, the transformation of
\code{type == "subexponential"} must be used to obtain geometric ergodicity
(with any valid \eqn{R} and \eqn{p} and \eqn{b}).
}
\value{
a list containing the non-missing arguments of this function.
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
