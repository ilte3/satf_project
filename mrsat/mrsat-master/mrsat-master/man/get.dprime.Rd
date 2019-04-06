\name{get.dprime}
\alias{get.dprime}
\title{
get dprime
}
\description{
A function that calculates dprimes for each timing bin 
}
\usage{
get.dprime(x, signal.list, noise.list, is.index = FALSE,
  subset = (x$SeqFlag == 0 & x$Experiment == 1), 
  binmax = 16, 
  dprime.correction=c("extreme", "loglinear", "none"),
  negative.dprime = FALSE,
  bintail = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
data frame usually obtained as an output from \code{\link{get.bins}}
}
  \item{signal.list}{
a list specifying for each condition, which "Condition" in \code{x} correspond to signal. 
}
  \item{noise.list}{
a list specifying for each condition, which "Condition" in \code{x} correspond to noise. 
}
  \item{is.index}{
  a logical specifying whether the values in \code{signal.list} and \code{noise.list} should be treated as numeric indices of \code{x$Condition}. Default to \code{FALSE}.
  }
  \item{subset}{
an optional vector specifying a subset of data to be used in the calculating dprimes.
}
  \item{binmax}{
an integer specifying maximum number of bins to be retained. All the data exceeding the limit will be discarded.
}
  \item{dprime.correction}{
the method to be used for calculating the dprime when hit rate or false alarm rate is zero or one. If \code{dprime.correction="extreme"} which is the default, zeros are replaced with \code{0.5/n} and ones are replaced with \code{(n-0.5)/n}, where n is the number of signal or noise trials. If \code{dprime.correction="loglinear"}, hit and false alarm rates are calculated by adding .5 to the numerator (the number of hit/false alarm trials) and 1 to the denominator (the number of signal/noise trials). If \code{dprime.correction="none"}, no adjustment will be applied.  
}
  \item{negative.dprime}{
  logical. If \code{FALSE} (default), the values of zero is used in place of negative dprime values.
  }
  \item{bintail}{
an optional integer specifying the number of bins to be concatenated at the tail.
}
}

\value{
An object of class \code{mrsat.data}, which is a data frame containing following columns:
  \item{bin}{timing bins}
  \item{lags}{means of time points at which the responses occured within each timing bin}
  \item{dprimes}{calculated as \code{dprime = qnorm(hit) - qnomr(fa)}}
  \item{condition}{values based on experimental conditions defined in the argument \code{scalelist} to which the dprime values belong to. }
  
The output of \code{get.dprime} contains following additional columns:
  \item{hit}{hit rate (correct response to correct stimuli) in each bin}
  \item{hit.correction}{a character indicating the type of correction performed when calculating hit rate. }
  \item{hit.denom}{total number of response to correct stimuli in each bin}
  
  \item{fa}{false alarm rate (correct response to incorrect stimuli) in each bin}
  \item{fa.correction}{a character indicating the type of correction performed when calculating false alarm rate. }
  \item{fa.denom}{total number of responses to incorrect stimuli in each bin}

}
\references{
Stanislaw, H., & Todorov, N. (1999). Calculation of signal detection theory measures. Behavior Research Methods, Instruments, & Computers, 31, 137--149.
}

\examples{

#need to define the scale.list
exp1.signal <- list(noint = c(1,3), 
                obrel = c(5,8), 
                obrelsub = c(11, 14))

exp1.noise <- list(noint = c(2,4), 
               obrel = c(6, 7, 9, 10), 
               obrelsub = c(12, 13, 15, 16))

#then load data and tag bins

data(Auditory_demo)
s01.bins <- get.bins(Auditory_demo, auditory=TRUE)

# finally, obtain dprime
s01.dp <- get.dprime(s01.bins$bins, 
  signal.list = exp1.signal, noise.list = exp1.noise, 
  is.index = TRUE, binmax=14)


plot(s01.dp)
}

