\name{scintec1.nonprofile}
\alias{scintec1.nonprofile}
\title{
  read the single-level variables from a file in Scintec Format-1 format
}
\description{
  get timeseries of non-profile variables from a 
  a Scintec atmospheric profiler data format "Format-1" (extension ".mnd")
}
\usage{
scintec1.nonprofile(file, header = list(), vars = list())
}
\arguments{
  \item{file}{
    name (and path if necessary) of the file to read}
  \item{header}{
    optionally, instead of reading the header from the file again, 
    the output of a previous call to \code{\link{scintec1.header}}
    on the same file can be supplied via this argument.}
  \item{vars}{
    optionally, instead of reading the header from the file again, 
    the output of a previous call to \code{\link{scintec1.variables}}
    on the same file can be supplied via this argument.}
}
\value{
  list(vector); length(vector) = n (number of times contained in the file)\cr
                each vevtor contains a non-profile variable (see variable list)
                special varaible name "time" contains the times as POSIXct
}
\references{
Scintec APRun software manual
}
\author{
  Clemens Druee
}
\examples{
dat <- scintec1.nonprofile(system.file("extdata", "example.mnd", package="readmet"))
head(dat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
\keyword{ IO }
\keyword{ Scintec }
