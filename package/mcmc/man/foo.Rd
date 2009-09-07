\name{foo}
\docType{data}
\alias{foo}
\title{Simulated logistic regression data.}
\description{
  Like it says
}
\usage{data(foo)}
\format{
  A data frame with variables
   \describe{
    \item{x1}{quantitative predictor.}
    \item{x2}{quantitative predictor.}
    \item{x3}{quantitative predictor.}
    \item{y}{Bernoulli response.}
  }
}
\examples{
library(mcmc)
data(foo)
out <- glm(y ~ x1 + x2 + x3, family = binomial, data = foo)
summary(out)
}
\keyword{datasets}
