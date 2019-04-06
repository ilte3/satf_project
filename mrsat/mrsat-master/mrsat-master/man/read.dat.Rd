\name{read.dat}
\alias{read.dat}
\title{
Read .dat files for a MR-SAT experiment
}
\description{
A function to read a file with .dat extension which is originally used to store lags, dprimes, and condition obtained/calculated from an output of a MR-SAT experiment via C program (MRsat-datacrunch.c).
}
\usage{
read.dat(file, quiet = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
the name of the dat file from which the data are to be read. 
}
  \item{quiet}{
logical. If FALSE, internally used function \code{scan()} will print how many items have been read per line.
}
}
\details{
The .dat file should be in a specific format, as follows:
First line should contain a number representing data points to be scanned (a).
Then, (the number of data points per condition)-by-(the number of condition) matrix for "lags", "dprimes", and "condition" should follow.
}
\value{
A data frame containing "lags", "dprimes", and "condition".
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
