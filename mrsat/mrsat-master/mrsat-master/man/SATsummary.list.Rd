\name{SATsummary.list}
\alias{SATsummary.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summarize a list of SATcurve objects
}
\description{
A convenience function to summarize a list of SATcurve objects in to a data frame.
}
\usage{
SATsummary.list(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class \code{SATcurve}}
}
\details{
This is what this function does:
\code{do.call("rbind", lapply(x, summary.SATcurve))}
}
\value{
A data farme with following columns
    \item{asym}{asymptote of the fitted curves. If there is more than one conditions, more than one columns of asym are returned, named as asym1, asym2, ... asymN, where N = number of conditions. }
    \item{rate}{rate of the fitted curves. If there is more than one conditions, more than one columns of rate are returned, named as rate1, rate2, ... rateN, where N = number of conditions. }
    \item{incp}{intercept of the fitted curves. If there is more than one conditions, more than one columns of incp are returned, named as incp1, incp2, ... incpN, where N = number of conditions. }
    \item{R2}{R-squared of the model}
    \item{adjR2}{Adjucted R-squared of the model}
    \item{RMSEfit}{Root Mean Sqaured Error of the model}
    \item{method}{The optimazation method used}
    
    The rows are named according to the names of the supplied list \code{x}.
}

\seealso{
\code{\link{fit.SATcurve}}, \code{\link{summary.SATcurve}}
}


