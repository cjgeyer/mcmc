\name{logit}
\docType{data}
\alias{logit}
\title{Simulated logistic regression data.}
\description{
  Like it says
}
\usage{data(logit)}
\format{
  A data frame with variables
   \describe{
    \item{x1}{quantitative predictor.}
    \item{x2}{quantitative predictor.}
    \item{x3}{quantitative predictor.}
    \item{x4}{quantitative predictor.}
    \item{y}{Bernoulli response.}
  }
}
\examples{
library(mcmc)
data(logit)
out <- glm(y ~ x1 + x2 + x3 + x4, family = binomial, data = logit)
summary(out)
}
\keyword{datasets}
