\name{acp}
\alias{acp}
\title{
Optimization using an iterative hill-climbing algorithm
}
\description{
Box-constrained optimization using an iterative hill-climbing algorithm
}
\usage{
acp(start, objective, ..., control = list(), lower, upper)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
Initial staring values for the parameters to be optimized
}
  \item{objective}{
A function to be minimized. The fist argument to the function should be the vector of parameters to be optimized. The function should return a scalar value.
}
  \item{\dots}{
Further arguments to be supplied to \code{objective}.
}
  \item{control}{
A list of control parameters. See 'Details' for more information.
}
  \item{lower, upper}{
Vectors, replicated to be as long as start. Lower and upper bounds of the parameter space. 
}

}
\details{
The \code{control} argument is a list that can supply any of the
  following components:
  \describe{
    \item{\code{trace}}{logical. If true, tracing information on the
      progress of the optimization is produced.}
    \item{\code{maxit}}{The maximum number of iterations. Defaults to \code{1000}}
    \item{\code{lmax}}{The number of stepping cycles per iteration. Defaults to \code{20}}
    \item{\code{strict}}{logical. If \code{TRUE}, iteration continues until all the values of stepsize are less than the value of \code{abs.tol}. Defaults to \code{FALSE}}
    \item{\code{abs.tol}}{The absolute convergence tolerance. Defaults to \code{.Machine$double.eps}. Relevant only when \code{strict=TURE}}
    \item{\code{auto.correct}}{logical indicating wheather to adjust the starting value with warning when value is equal to upper or lower bound (which may lead to convergence error). Defaults to \code{TRUE}}
    }
    
    
}
\value{
A list containing:
 \item{par}{The best set of parameters found through iterative searching.}
 \item{convergence}{An integer indicating success or possible error. 0 indicates successful completion.}
 \item{iterations}{The number of iterations took took to converge.}
 \item{value}{The value of \code{objective} with \code{par}.}
 \item{trace}{If argument \code{control(trace=TRUE)}, a dataframe containing the values and pars at each iteration is returened. Otherwise, \code{NULL}}
}
\references{
Chandler, J.P. (1969).  Subroutine STEPIT -- finds local minimum of a smooth function of several parameters.  Behavioral Science, 14, 81-82.

Judd, C. M., & McClelland, G. H. (1989). Data analysis: A model-comparison approach. San Diego: Harcourt Brace Jovanovich. 
}
\author{
\R port: Kazunaga Matsuki

original Fortran code by Unknown
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[stats]{optim}}, \code{\link[stats]{nlm}}, \code{\link[stats]{nlminb}}
}
\examples{
## example taken from nlmimb
x <- rnbinom(100, mu = 10, size = 10)
hdev <- function(par)
    -sum(dnbinom(x, mu = par[1], size = par[2], log = TRUE))
    
nlminb(c(20, 20), hdev, lower = 0.001, upper = Inf) 

## acp produce comparable results
acp(c(20, 20), hdev, lower = 0.001, upper = Inf) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}
