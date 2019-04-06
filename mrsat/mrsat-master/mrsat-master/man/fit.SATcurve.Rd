\name{fit.SATcurve}
\alias{fit.SATcurve}
\alias{print.SATcurve}
\title{
fit SAT curves
}
\description{
A function to fit curve to d primes obtained from MR-SAT experiment.
}
\usage{
fit.SATcurve(data, par.cond, 
             fix.asym = FALSE,
             params = get.param(par.cond, auto.asym=fix.asym, data=data),
             opt = c("acp" ,"nlminb", "nlm", "optim"), 
             multi.opt = FALSE,
             rep=1, quiet=TRUE, maxit=500, strict=TRUE, trace=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
data frame contining columns named \code{dprime}, \code{lags}, and \code{condition} (as in the output of \code{\link{get.dprime}}), or character specifying the location of dat file (see \code{\link{read.dat}} for detail of the file format)
}
  \item{par.cond}{
a list specifying the configuraiton of parameters per condition. Should contain three named elements: \code{asym}, \code{rate}, and \code{incp}. See 'Details'.
}
  \item{fix.asym}{
a logical specifying whether asymptote parameters of the model should be fixed to the starting values supplied in the  \code{params} argument.
  }
  \item{params}{
a data frame specifying parameters, containing columns named \code{start}, \code{lower}, and \code{upper}. Defalt to \code{get.param(par.cond)}, see \code{\link{get.param}}. It can also be a character vector specifying the location of 'par' file (see \code{\link{read.param}} for detail of the file format).
}
  \item{opt}{
charaster or a character vector indication the optimizer to be used. If \code{length(opt) > 1} and \code{multi.opt=TRUE}, all of the specified methods will be used in the specified order. If \code{length(opt) > 1} and \code{multi.opt=FALSE}, the first element of \code{opt} will be used.
}
  \item{multi.opt}{
logical, specifying wheather to fit the data using multiple optimizer.
}
  \item{rep}{
numeric, representing number of times the optimazation process is repeated. If \code{multi.opt = TRUE}, each optimizer specified in \code{opt} is run as many times as \code{rep}.
}
  \item{quiet}{
logical. Default to \code{TRUE}. If \code{FALSE}, a progress bar is printed.
}
  \item{maxit}{
The maximun number of iteration. Defaults to \code{500}.
}
  \item{strict}{
logical. Only relevant if \code{method="acp"}. See 'Details' of \code{\link{acp}}
}
  \item{trace}{
logical. Only relevant if \code{method="acp"}. See 'Details' of \code{\link{acp}}
}
}
\details{

The argument \code{par.cond} specifies whether or not and how certain parameter(s) are shared across conditions (which are specified in \code{condition} column of \code{data}). The argument \code{par.cond} can take either list or matrix, as long as it meets following specification. If it is a list, then it should be a named list containing following elements: \code{asym}, \code{rate}, and \code{incp}. Each element of the list must be of same length, being equal either to \code{length(condition)} or to \code{length(levels(condition))} (or \code{length(unique(condition))} if \code{condition} is not a factor). If it is a matrix, then it should have \code{asym}, \code{rate}, and \code{incp} as its columns, and rows should represent different condiitons.

As an example, there are three elements in \code{par} (i.e., \code{length(par) = 3}), and there are three conditions so that the three conditions share all of \eqn{\lambda} (asymptote), \eqn{\beta} (rate) and \eqn{\delta} (intercept). In this case, \code{par.cond} should be specified as \code{list(asym=c(1,1,1), rate=c(1,1,1), incp=c(1,1,1))}.

As another example, there are six elements in \code{par}, consisting of three \eqn{\lambda} (asymptote; one for each condition) , two \eqn{\beta} (rate; one for condition 1 and anothor for condition 2 and 3), and one \eqn{\delta} (intercept; for all the conditions). Here, \code{par.cond = list(asym=c(1,2,3), rate=c(1,2,2), incp=c(1,1,1))}



}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{fit}{a list of object returned from the optimizer (if multiple optimizers were used, the one that returned the smallest RMSE), with a few additinal information}
\item{SSE}{Sum of Squared Error of the model.}
\item{MSE}{Mean Squared Error of the model.  A object of class \code{logLik}.}
\item{RMSE}{Root Mean Squared Error of the model.}
\item{logLik}{Log likelihood of the model. A object of class \code{logLik}.}
\item{R2}{R squared of the model}
\item{adjR2}{adjusted R squared of the model}
\item{data}{same as the argument \code{data}}
\item{par.cond}{same as the argument \code{par.cond}}
\item{rep.out}{If \code{rep>1}, the data.frame containing the parameters and fit statistics obtained from each repeated fitting of the data. If \code{rep=1}, \code{NULL}}
   
}

\seealso{
\code{\link{acp}}, \code{\link{optim}}, \code{\link{nlminb}}, \code{\link{nlm}}
}
\examples{

#load data
data(Auditory_demo)

#need to define the scale.list
exp1.signal <- list(noint = c(1,3), 
                obrel = c(5,8), 
                obrelsub = c(11, 14))

exp1.noise <- list(noint = c(2,4), 
               obrel = c(6, 7, 9, 10), 
               obrelsub = c(12, 13, 15, 16))


# tag bins
s01.bins <- get.bins(Auditory_demo, auditory=TRUE, window = "fixed")

#calculate dprimes
s01.dp<- get.dprime(s01.bins$bins, 
  signal.list = exp1.signal, noise.list = exp1.noise, is.index = TRUE,
  binmax=14)

# specify par.cond
my.pc <- list(asym=c(1,2,3), rate=c(1,2,2), incp=c(1,1,1))

# fit the data
s01.fit <- fit.SATcurve(s01.dp, par.cond = my.pc)

#plot the data and fitter curves
plot(s01.fit)

#to calculate AIC and BIC, do
AIC(s01.fit)
BIC(s01.fit)

#you may also calculate AIC analogue using MSE
# as suggested by Burnham and Anderson (1998)
AIC(s01.fit$MSE)
}
