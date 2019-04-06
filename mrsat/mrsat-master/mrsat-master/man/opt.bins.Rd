\name{opt.bins}
\alias{opt.bins}
\alias{print.optbin}
\alias{plot.optbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the optimal way to bin the response times.
}
\description{
A function to divide the response time data into a fixed number of small consecutive intervals or bins, in a way that is optimal given the distribution of the reponse time data. This involves finding an optimal offset (or start) of bining, and possily an optimal bin size. 
}
\usage{
opt.bins(RT, maxRT = 6300, breaksize = 7, 
  startbin = 1, endbin = 16, binsize = 350, 
  type = c("fixed", "free"),
  add.two=FALSE, shift.edge = FALSE)
}

\arguments{
  \item{RT}{
a numeric vector containing the raw response time information
}
  \item{maxRT}{
numeric. Maximum duration of response time above which will be trimmed for the estimation.
}
  \item{breaksize}{
a single number representing the size of bins (in milisecond) used in the intermediate step.
}
  \item{startbin}{
integer. 
}
  \item{endbin}{
integer.
}
  \item{binsize}{
a single number corresponding to the initial (and final if \code{type="fixed"}) window size to be used for binning
}
  \item{type}{
a character string indicating whether to use the window of prespecified width ("fixed") or the widow of data-specific width ("free").
}
  \item{add.two}{
logical. An option for compatibility. If \code{TRUE}, a magic number 2 is added to the window size when calculating the fit.
  }
  \item{shift.edge}{
  logical. An option forcompatibility with old codes. If \code{TRUE}, the timing bins are shifted 100 ms forward. Default to \code{FALSE}.
  }
}
\details{
This function is meant to be called internally from a wapper function \code{\link{get.bins}}.

The optimal timing bin size (if \code{window="free"}) and offset are estimated by fitting the half-rectified sine wave against the histgram of response data.
}
\value{
A list containing:
\item{maxRT, breaksize, endbin, binsize, type}{same as scpecified in the arguemnts}
\item{peaktime}{a vector of time points representing peak in each time bin}
\item{offset}{difference in time between the offset of the stimulus and the onset of the fastest response.}
\item{fit}{root mean squared distance between the actual distribution of the data and the half-rectified sine wave}
\item{phase}{phase (in ms) for the sign wave function.}
\item{omega}{fitted window size. If \code{type="fixed"}, the value should be the same as the argument \code{binsize}}
\item{rt.hist}{an object of class "histogram", with \code{x=RT}, and \code{breaks=seq(from=0, to=maxRT, by=breaksize)}}
\item{amp}{amplitude for the sign wave function.}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{get.bins}}
}
\examples{
data(Auditory_demo)
rt_col <- grep("rt[0-3][0-9]",colnames(Auditory_demo), value=TRUE)
RT <- Auditory_demo[,rt_col]
d.opt <- opt.bins(unlist(na.omit(RT)))
plot(d.opt)
}
