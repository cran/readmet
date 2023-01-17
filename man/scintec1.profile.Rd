\name{scintec1.profile}
\alias{scintec1.profile}
\title{
  read the multi-level (profile) variables from a Scintec Format-1 file
}
\description{
  get profile variables from
  a Scintec atmospheric profiler data format "Format-1" (extension ".mnd")
}
\usage{
scintec1.profile(file, header = list(), vars = list())
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
  list(matrix); dim(matrix)=c(n,m) \cr
  where n is number of times contained in the file and m is 
                the number of levels\cr
                i.e. rows represents timeseries of at a fixed level,
                columns represent profiles at a fixed time
}
\references{
Scintec APRun software manual
}
\author{
  Clemens Druee
}
\examples{
dat <- scintec1.profile(system.file("extdata", "example.mnd", package="readmet"))
dat
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
\keyword{ IO }
\keyword{ Scintec }
