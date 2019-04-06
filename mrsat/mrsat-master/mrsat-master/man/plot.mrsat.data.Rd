\name{plot.mrsat.data}
\alias{plot.mrsat.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the dprimes obtained using \code{get.dprime}
}
\description{
Plot the dprimes obtained using \code{\link{get.dprime}} (and also \code{\link{average.dprime}})
}
\usage{
\method{plot}{mrsat.data}(x, condition, xlim, ylim, legend=TRUE, pch, col, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object of class \code{mrsat.data}
}
  \item{condition}{
a vector specifying the condition(s) to be included in the plot. The values should be either numeric index or part of condition labels used in the SAT curve object.
}
  \item{xlim}{
the x limits of the plot
}
  \item{ylim}{
the y limits of the plot
}
  \item{legend}{
  logical. If \code{TRUE}, legends will be printed.
}
  \item{pch}{
  a vector of plotting characters or symbols per condition.
  }
  \item{col}{
  The colors for lines and points per condition.
  }
  \item{...}{
other graphical parameters to be passed.
  }
}

\seealso{
\code{\link{get.dprime}}, \code{\link{average.dprime}}, \code{\link{fit.SATcurve}}, \code{\link{plot.SATcurve}}, 
}
\examples{
#see help(get.dprime) or help(average.dprime) for examples
}

