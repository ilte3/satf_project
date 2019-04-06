\name{mrsat-package}
\alias{mrsat-package}
\alias{mrsat}
\docType{package}
\title{
Multiple Response Speed-Accuracy Tradeoff (MR-SAT).
}
\description{
This package contains a set of functions useful for analyzing data for psychology experiments based on Multiple Response Speed-Accuracy Tradeoff (MR-SAT) Method (Reed, 1973, 1976; McElree, 1993).

The much of functions in this package is based on the script written by Matthew Wagers called "mrsatfxns.R", which were written specifically to analyze the data from MR-SAT experiments collected via the McElree Lab E-prime scripts. 
}
\details{
\tabular{ll}{
Package: \tab mrsat\cr
Type: \tab Package\cr
Version: \tab 0.1.1\cr
Date: \tab 2015-08-21\cr
License: \tab GPL(>=2)\cr
}


~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Julie Van Dyke, Matt Wagers, Pyeongwhan Cho, Kazunaga Matsuki
Maintainer: Kazunaga Matsuki <matsukk@mcmaster.ca>
}
\references{
Reed, A. V. (1973). \emph{Speed-accuracy trade-off in recognition memory.} Science, 181, 574--576.
Reed, A. V. (1976). \emph{The time course of recognition in human memory.} Memory & Cognition, 4, 16--30.
}
\keyword{ package }
\examples{

#loading the demo data
data(Auditory_demo)

## define how conditions are grouped
my.signal <- list(noint = c(1,3), 
                obrel = c(5,8), 
                obrelsub = c(11, 14))
my.noise <- list(noint = c(2,4), 
               obrel = c(6, 7, 9, 10), 
               obrelsub = c(12, 13, 15, 16))

## get bins
d.bins <- get.bins(Auditory_demo, auditory=TRUE)

#check how the RT are binned.
plot(d.bins$opt.bins)

## obtain dprime values
d.dprime <- get.dprime(d.bins$bins, signal.list=my.signal, noise.list=my.noise, 
    is.index=TRUE, binmax=14)


## define structure of parameters
## in this case, different parameters for each condition
pc333 <- list(asym=c(1, 2, 3), rate=c(1, 2, 3), incp=c(1, 2, 3))

## fit the curve assuming, and plot
fit.333 <- fit.SATcurve(d.dprime, par.cond = pc333)
plot(fit.333, main="333")

## compare that to the curves with different asymptotes but the same rate and intercept

pc311 <- list(asym=c(1, 2, 3), rate=c(1, 1, 1), incp=c(1, 1, 1))
fit.311 <- fit.SATcurve(d.dprime, par.cond = pc311)
plot(fit.311, main="311")

#compare outputs of the two models side by side
SATsummary.list(list(fit.333, fit.311))

#or just compare AIC of the two models
AIC(fit.333, fit.311)

#fitting a 311 model with fixed asymptote
fit.311fa <- fit.SATcurve(d.dprime, fix.asym=TRUE, par.cond = pc311)
#and compare with the other two models
SATsummary.list(list(fit.333, fit.311, fit.311fa))
}
