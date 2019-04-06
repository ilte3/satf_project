\name{summary.SATcurve}
\alias{summary.SATcurve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize SATcurve object
}
\description{
Summary method for an "SATcurve" object
}
\usage{
\method{summary}{SATcurve}(object, ...)
}

\arguments{
  \item{object}{
  object of class \code{SATcurve}
}
  \item{...}{other arguemnts to be passed. Not currently used.}
}
\value{

THe value returned is essentially the same as that returned by \code{\link{fit.SATcurve}}, but is rearrenged into a data farme with following columns:
    \item{asym}{asymptote of the fitted curves. If there is more than one conditions, more than one columns of asym are returned, named as asym1, asym2, ... asymN, where N = number of conditions. }
    \item{rate}{rate of the fitted curves. If there is more than one conditions, more than one columns of rate are returned, named as rate1, rate2, ... rateN, where N = number of conditions. }
    \item{incp}{intercept of the fitted curves. If there is more than one conditions, more than one columns of incp are returned, named as incp1, incp2, ... incpN, where N = number of conditions. }
    \item{R2}{R-squared of the model}
    \item{adjR2}{Adjucted R-squared of the model}
    \item{RMSEfit}{Root Mean Sqaured Error of the model}
    \item{method}{The optimazation method used}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fit.SATcurve}}, \code{\link{SATsummary.list}}
}

