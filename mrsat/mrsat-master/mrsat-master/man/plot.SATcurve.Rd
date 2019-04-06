\name{plot.SATcurve}
\alias{plot.SATcurve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot methods for the object of class \code{SATcurve}
}
\description{
Plot the data points and curve(s) that were fit using \code{fit.SATcurve}
}
\usage{
\method{plot}{SATcurve}(x, condition, 
  xlim, ylim, npoints = 100, 
  xlab = "Time", ylab="d\'",
  legend=TRUE, labels = TRUE, lab.cex = 1,
  residual = FALSE, pch, col, lwd = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object of class \code{SATcurve}
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
  \item{npoints}{
number of data points to use for plotting a curve. Larger value means smoother curve. Default to \code{100}.
}
  \item{xlab}{a label for the x axis, default to "Time"}
  \item{ylab}{a label for the y axis, default to "d'""}
  \item{legend}{
  logical. If \code{TRUE}, legends (and fitted parameters of the models) will be printed.
}
  \item{labels}{
logical. If \code{TRUE}, Adjusted R-squared and RMSE of the fitted model will be printed.
}
\item{lab.cex}{
cex value for the labels. Default to 1.
}
  \item{residual}{
logical. If \code{TRUE}, a residual plot instead of the fitted curves are drawn. Default to \code{FALSE}
}
  \item{pch}{
  a vector of plotting characters or symbols per condition.
  }
  \item{col}{
  The colors for lines and points per condition.
  }
  \item{lwd}{
  The width of the line use to draw the curve. Defaul to \code{1}.
  }
  \item{...}{
other graphical parameters to be passed.
  }
}
\seealso{
\code{\link{fit.SATcurve}}, \code{\link{plot.mrsat.data}}, 
}
\examples{
#see help(fit.SATcurve) for examples
}