\name{scintec1.variables}
\alias{scintec1.variables}
\title{
  read the variable descriptions from a file in Scintec Format-1 format
}
\description{
  get timeseries of non-profile variables from a 
  a Scintec atmospheric profiler data format "Format-1" (extension ".mnd")
}
\usage{
scintec1.variables(file, header = list())
}
\arguments{
  \item{file}{
    name (and path if necessary) of the file to read}
  \item{header}{
    optionally, instead of reading the header from the file again, 
    the output of a previous call to \code{\link{scintec1.header}}
    on the same file can be supplied via this argument.}
}
\value{
  data.frame; each entry corresponds to one variable. The columns are:
  \item{label}{Name of variable}
  \item{symbol}{short name; corresponds to list name in 
               scintec1.profile,scintec1.nonprofile, and read.scintec1}
  \item{unit}{physical unit}
  \item{type}{code describing for example averaging, profile/nonprofile, 
              measured, derived or assimilated. See Scintec Software Manual}
  \item{error.mask}{See Scintec Software Manual}
  \item{gap.value}{value in files representing missing values of this variable}
}
\references{
Scintec APRun software manual
}
\author{
  Clemens Druee
}
\examples{
vars <- scintec1.variables(system.file("extdata", "example.mnd", package="readmet"))
head(vars)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
\keyword{ IO }
\keyword{ Scintec }
