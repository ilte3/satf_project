\name{get.param}
\alias{get.param}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get a data frame of starting values and bounds for fitting SAT curves.
}
\description{
Construct a data frame containing starting values and upper/lower bounds for fitting SAT curves having certain constraints.
}
\usage{
get.param(par.cond, 
  asym = c(1, 0.5, 4.7), 
  rate = c(3, 0.1, 12), 
  incp = c(0.25, 0.1, 2), 
  random.start = FALSE, 
  auto.asym = FALSE, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par.cond}{
a list specifying the configuraiton of parameters per condition. Should contain three named elements: \code{asym}, \code{rate}, and \code{incp}. See 'Details'.
}
  \item{asym}{
a numeric vector specifying the starting value, minimum and maximum values for asymptote
}
  \item{rate}{
 a numeric vector specifying the starting value, minimum and maximum values for rate.
}
  \item{incp}{
 a numeric vector specifying the starting value, minimum and maximum values for intercept.
}
  \item{random.start}{
  logical. If \code{TRUE}, starting values will be based on random numbers that fall between lower and upper bounds.
  }
  \item{auto.asym}{
  logical. If \code{TRUE}, values for the asymptote parameters are automatically calculated as the average of the last 4 dprime values specified in the \code{data}.
  }
  \item{data}{
  data.frame of class \code{mrsat.data}
  }
}
\details{
This is a convinience function to generate a data frame of inital starting parameter values as well as lower and upper bounds that need to be passed when using the function \code{fit.SATcurve}. The default values are set based on prior studies, but are somewhat arbitrarily chosen.

}
\value{
A data frame containing:
  \item{start}{the starting values for each parameters to be optimized}
  \item{lower}{the lower bound for each parameters to be optimized}
  \item{upper}{the upper bound for each parameters to be optimized}
%% ...
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fit.SATcurve}}
}
\examples{
pc1 <- list(asym=c(1,2,3), rate=c(1,2,2), incp=c(1,1,1))
get.param(pc1)

pc2 <- list(asym=c(1,2,3,4,5), rate=c(1,2,2,3,3), incp=c(1,1,1,2,2))
get.param(pc2)
}

