\name{average.dprime}
\alias{average.dprime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Average dprime and lag of each bins
}
\description{
A convenience function to aggregate the  dprime and lag values from a list of multiple data frames obtained through \code{get.dprime} function.
}
\usage{
average.dprime(bin.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bin.list}{
a list containing multiple data frames, each of which are separate output of the function \code{\link{get.dprime}}
}
}

\value{
An object of class \code{mrsat.data}, which is a data frame containing following columns:
  \item{bin}{timing bins}
  \item{lags}{means of time points at which the responses occured within each timing bin}
  \item{dprimes}{calculated as \code{dprime = qnorm(hit) - qnomr(fa)}}
  \item{condition}{values based on experimental conditions defined in the argument \code{scalelist} to which the dprime values belong to. }

}


\seealso{
\code{\link{get.dprime}}, \code{\link{plot.mrsat.data}}
}
\examples{
#need to define the scale.list
exp1.signal <- list(noint = c(1,3), 
                obrel = c(5,8), 
                obrelsub = c(11, 14))

exp1.noise <- list(noint = c(2,4), 
               obrel = c(6, 7, 9, 10), 
               obrelsub = c(12, 13, 15, 16))

#load data

data(Auditory_demo)

### this isn't really meaningful, but for the purpose of demo,
### tag bins in two diffent ways (with "fixed" vs. "free" window)
### and then obtain two data frames containing slightly different
### dprime and lag values.

s01.bins.fixed <- get.bins(Auditory_demo, auditory=TRUE, window = "fixed")
s01.bins.free <- get.bins(Auditory_demo, auditory=TRUE,  window = "free")


s01.dp.fixed <- get.dprime(s01.bins.fixed$bins, 
  signal.list = exp1.signal, noise.list = exp1.noise, is.index=TRUE,
  binmax=14)

s01.dp.free <- get.dprime(s01.bins.free$bins, 
  signal.list = exp1.signal, noise.list = exp1.noise, is.index=TRUE,
  binmax=14)
  
### average the two data frames

mean.dp <- average.dprime(list(s01.dp.fixed, s01.dp.free))

plot(mean.dp)
}
