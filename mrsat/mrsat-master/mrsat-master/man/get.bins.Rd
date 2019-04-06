\name{get.bins}
\alias{get.bins}
\title{
Get bins
}
\description{
A convenience function to transform raw data from MR-SAT expeiment to list that contain information about bins used and a long-format data frame with bins specified. 
}
\usage{
get.bins(x, 
  rt_col = "rt[0-3][0-9]",
  res_col = "key[0-3][0-9]",
  window = c("fixed","free"),
  offset = TRUE, 
  sauce=FALSE, 
  endbin = 16, binsize = 350, 
  auditory = FALSE, audio.rt.base = 5006, 
  quiet = TRUE, 
  plot = FALSE, 
  add.two = FALSE,
  shift.edge = FALSE)
}
\arguments{
  \item{x}{
a data frame. See 'Details'.
}
  \item{rt_col}{
Either a numeric vector specifying the colum position of columns in \code{x} representing the response time information, or a charactor string containing a regular expression to be matched in the column names of \code{x}.
}
  \item{res_col}{
Either a numeric vector specifying the colum position of columns in \code{x} representing the reponse accuracy information, or a charactor string containing a regular expression to be matched in the column names of \code{x}.
}
  \item{window}{
character. If \code{fixed}, a fixed window with size specified in \code{binsize} will be used for binning. If \code{free}, the window size will be determined based on the pattern in the data.
}
  \item{offset}{
logical. If \code{TRUE}, a initial period of no response will be assumed.
}
  \item{sauce}{
logical. Currently unknown what this is for but it's there for histrical reasons.
}
  \item{endbin}{
Numeric. 
}
  \item{binsize}{
Numeric. The size of time bins in millisecond.
}
  \item{auditory}{
logical. If \code{TRUE}, reaction times in the \code{x} data frame will be adjusted based on the constrant values specified by \code{audio.rt.base} and the value speficied in the 'SentLength' column of \code{x}.
}
  \item{audio.rt.base}{
numeric.
}
  \item{quiet}{
logical. If FALSE and when the \code{offset} is not 'none', information regarding the parameters chosen for the offset will be printed.
}
  \item{plot}{
logical. If \code{TRUE}, a histram of the response data will be plotted along with the   
}
  \item{add.two}{
  logical. Setting this to \code{TRUE} make the results more compatible with the old codes. Only relevant when \code{window="fixed"}. Default to \code{FALSE}.
  }
  \item{shift.edge}{
  logical. Setting this to \code{TRUE} make the results more compatible with the old codes by shiting the timing bins 100 ms forward. Default to \code{FALSE}.
  }
}
\value{
%%  ~Describe the value returned
A list containing
\item{bins}{a data frame containing the information needed to calculate dprimes using \code{get.dprime}.}
\item{opt.bin}{An output form \code{opt.bins}}
%% ...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(Auditory_demo)
s01.bins <- get.bins(Auditory_demo, auditory = TRUE)

#this would provide a slightly different bins 
s01.bins.free <- get.bins(Auditory_demo, auditory = TRUE, window="free")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
