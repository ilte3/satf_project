\name{read.edat.txt}
\alias{read.edat.txt}
\title{
Read a text file exported from E-prime data file.
}
\description{
read a tab-delimited text file exported through E-DataAid, a part of E-prime application suite. The function is essentially a wrapper around the function \code{read.delim}
}
\usage{
read.edat.txt(file, collapse = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
The name of the file from which the data are to be read. Can be a vector of multile file names.
}
  \item{collapse}{
logical. Relevant only when the \code{file} contains more than one file names. If \code{TRUE}, the ouput will be a single data.frame. If \code{FALSE}, the list of data.frame will be returned.
}
}
\details{
This function is essentially a wrapper around the function function \code{read.delim}.


}
\value{
A data frame, or a list of data frames (if \code{collapse=TURE}), containing the data from the file.
}
\keyword{file}
