\name{logLik.SATcurve}
\alias{logLik.SATcurve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-Likelihood of a 'SATcurve' Object
}
\description{
Extract the log-likelihood of the fitted SAT curve.
}
\usage{
\method{logLik}{SATcurve}(object, useMSE = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class \code{SATCurve}
}
  \item{useMSE}{
logical. If TRUE, MSE instead of logLik is retured.
}
  \item{\dots}{
further arguments (currently unused)
}
}
\value{
an object of class \code{logLik}, with the value corresponding to the log-likelihood (or MSE if \code{useMSE=TURE} of the fitted SAT curve.
}
\keyword{models}
