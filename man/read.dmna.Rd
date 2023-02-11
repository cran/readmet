\name{read.dmna}
\alias{read.dmna}
\title{
  function to read contents of a file in DMNA-format
}
\description{
  reads gridded one to two dimensional data from a DMNA format
  that used for example
  by the German national dispersion model AUSTAL
}
\usage{
read.dmna(file, val=1, debug=FALSE)
}
\arguments{
  \item{file}{
    name (and path if necessary) of the file to read}
  \item{val}{
    name (number of name of the variable to be extracted. Only needed in case the
    file contains gridded data and more than one variable)}
  \item{debug}{
    if TRUE shwon debugging information}
}
\value{
  returns an object whos type depends on the data containt in the file.  
  
  If the file contains a timeseries, a data.frame is returned that contains
  each variable in the file as a column. Column names are taken from the
  format description embedded in the file. If a column for time (name "te")
  is present, times are returned as POSIXct.
  
  If the file contains gridded data an array of numeric is returned.  The number
  of dimensions of the array is the lowest needed to hold the data.  If the file
  contains mor than one variable, the function will return only one variable with
  each call.  Choose the varibale by providing paramter val (name or number of 
  variable). 
}
\references{
AUSTAL 3.1 model reference (by Ingenieurbuero Janicke)
}
\author{
  Clemens Druee
}
\note{
  Prior to version 1.2.4, this function did not arrays with more than
  two dimensions, or multiple variables in files of more than one
  dimension. Prior to version 1.5.0 timeseries were not implemented.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
so2 <- read.dmna(system.file("extdata", "example-grid.dmna", package="readmet"))
ax <- dmna.axes(system.file("extdata", "example-grid.dmna", package="readmet"))
filled.contour(ax$x, ax$y, so2)

monitor <- read.dmna(system.file("extdata", "example-monitor.dmna", package="readmet"))
plot(monitor$te, monitor$X2, type="l")

series <- read.dmna(system.file("extdata", "example-series.dmna", package="readmet"))
plot(series$te, series$ua, type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
\keyword{ file }
\keyword{ DMNA }
\keyword{ AUSTAL }
