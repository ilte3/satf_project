\name{extractAIC}
\alias{extractAIC.SATcurve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
extractAIC method for class 'SATcurve'
}
\description{
Extract Akaike's An Information Criteria from a 'SATcurve' object
}
\usage{
\method{extractAIC}{SATcurve}(fit, scale, k = 2, useMSE = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{fit}{object of class \code{glc}, \code{gqc}, \code{gcjc}, or \code{grg}}
 \item{scale}{unused argument}
 \item{k}{numeric specifying the penalty per parameter to be used in calculating AIC. Default to 2.}
 \item{useMSE}{logical. If \code{TRUE}, MSE is used to calculate Burnham and Anderson's (1998) AIC analog. Default to \code{FALSE}}
 \item{\dots}{further arguments (currently not used).}
}

\details{
As with the default method, the criterion used is
  \deqn{AIC = - 2\log L +  k \times \mbox{df},}{AIC = - 2*log L +  k * df,}
  where \eqn{L} is the likelihood and \eqn{df} is the degrees
  of freedom (i.e., the number of free parameters) of \code{fit}.
}
\value{
  A numeric vector of length 2 including:
  \item{df}{the degrees of freedom for the fitted model \code{fit}.}
  \item{AIC}{the Akaike's Information Criterion for \code{fit}.}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
