\name{scintec1.header}
\alias{scintec1.header}
\title{
  read the header of a file in Scintec Format-1 format
}
\description{
  get timeseries of non-profile variables from a 
  a Scintec atmospheric profiler data format "Format-1" (extension ".mnd")
}
\usage{
scintec1.header(file)
}
\arguments{
  \item{file}{
    name (and path if necessary) of the file to read
  }
}
\value{
  named list; each entry corresponds to one header field:
  \item{starttime}{starting time of the measurement as POSIXct object}
  \item{filecount}{running number of files produced during the current measurement}
  \item{instrument}{model type of the device}
  \item{commentlines}{lines of text containing the coment list}
  \item{variables}{lines of text containing the list of the variables in encoded form}
  \item{heightlevels}{number of height of levels where data are produced}
}
\references{
Scintec APRun software manual
}
\author{
  Clemens Druee
}
\examples{
hdr <- read.scintec1(system.file("extdata", "example.mnd", package="readmet"))
hdr$instrument
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
\keyword{ IO }
\keyword{ Scintec }
