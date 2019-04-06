\name{vareq}
\alias{vareq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to return expected values of an SAT curve
}
\description{
A function to return expected values of an SAT curve given the speficified parameter at given time points
}
\usage{
vareq(par, lags, condition, 
  par.cond = list(asym = condition, rate = condition, incp = condition))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
a numeric vector of parameters that determine the shape of SAT curve. Should be ordered as follows: parameter(s) for asymptote, rate and then intercept. If there is more than one parameter in any of the three gruop, the values should be ordered within the group according to the condition to which the values belong (See Details).
}
  \item{lags}{
a numeric vector of time points at which expected values need be calculated.
}
  \item{condition}{
vector specifying condition to which the curve(s) belong. Must be the same length as \code{length(lags)}.
}
  \item{par.cond}{
list or matrix, containing following named elements (in case of list) or columns (in case of matrix): \code{asym}, \code{rate}, and \code{incp}. The supplied list/matrix should specify the configuration of each parameter for different conditions speficied in \code{condition} (See 'Details'). Default to have all the parameters different for each of the conditions specified in \code{condition|}.
}
}
\details{
The SAT curve in this function is given by the following equation:
\deqn{d'(t) = \lambda (1 - exp(-\beta (t - \delta)), for t < \delta else 0}
where \eqn{\lambda} serves to estimate the asymptotic level of performance, \eqn{\delta}
represent the intercept of the curve on x axes, and \eqn{\beta} indexes the rate at which curve turns to asymptote.

If there are more than one condition, then both \code{condition} and \code{par.cond} need to be specified so that values in \code{par} are applied to the equation appropriately. 

The argument \code{par.cond} specifies whether or not and how certain parameter(s) are shared across conditions. The argument \code{par.cond} can take either list or matrix, as long as it meets following specification. If it is a list, then it should be a named list containing following elements: \code{asym}, \code{rate}, and \code{incp}. Each element of the list must be of same length, being equal either to \code{length(condition)} or to \code{length(levels(condition))} (or \code{length(unique(condition))} if \code{condition} is not a factor). If it is a matrix, then it should have \code{asym}, \code{rate}, and \code{incp} as its columns, and rows should represent different condiitons.

As an example, there are three elements in \code{par} (i.e., \code{length(par) = 3}), and there are three conditions so that the three conditions share all of \eqn{\lambda} (asymptote), \eqn{\beta} (rate) and \eqn{\delta} (intercept). In this case, \code{par.cond} should be specified as \code{list(asym=c(1,1,1), rate=c(1,1,1), incp=c(1,1,1))}.

As another example, there are six elements in \code{par}, consisting of three \eqn{\lambda} (asymptote; one for each condition) , two \eqn{\beta} (rate; one for condition 1 and anothor for condition 2 and 3), and one \eqn{\delta} (intercept; for all the conditions). Here, \code{par.cond = list(asym=c(1,2,3), rate=c(1,2,2), incp=c(1,1,1))}

}
\value{
a numeric vector of length equal to \code{length(lags)}. The values represent the expected values at time points specified in \code{lags}, given the curves speficied with the parameter in \code{par}.
}

\seealso{
\code{\link{fit.SATcurve}}
}
\examples{

#3 conditions, 10 time points from 0.1 to 5
my.cond <- rep(1:3, each=10)
my.lags <- rep(seq(0.1,5, length.out=10), 3)

#some random parameters with 6 elewments for the purpose of demo
my.par <- c(4, 3, 2, 2, 2, 1)

#par.cond configuration according to the second example in the 'Detail'
my.pc <- list(asym=c(1,2,3), rate=c(1,2,2), incp=c(1,1,1))

#then get the expected value

vareq(par=my.par, lags = my.lags, condition = my.cond, par.cond=my.pc)
}

