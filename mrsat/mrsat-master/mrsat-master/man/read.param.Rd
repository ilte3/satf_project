\name{read.param}
\alias{read.param}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read a parameter file
}
\description{
A function to read the parameters speficied in the 'par' file format specific to the original 'acp.f' code.
}
\usage{
read.param(file, quiet = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
the name of the file from which the parameters are to be read.
}
  \item{quiet}{
logical. If FALSE, it will print a message saying how many items have been read per line.
}
}
\details{
The file with the file name specified in the \code{file} should contain following information for each line, separated by comma where multiple values are present:
\enumerate{
\item (1st line): the number of parameters
\item (2nd line): the starting value for each parameter
\item (3rd line): the lower limit for each parameter
\item (4th line): the upper limit for each parameter
}


}
\value{
A data frame containing starting values and their corresponding upper and lower limits.
}

\author{
Pyeongwhan Cho, Kazunaga Matsuki
}


\seealso{
\code{\link{read.dat}}
}

\keyword{file}
