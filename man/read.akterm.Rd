\name{read.akterm}
\alias{read.akterm}
\title{
  function to read contents of a file in AKTERM-format
}
\description{
  reads weather data time series from a AKTERM format file
  that is used for example
  by the German national dispersion model AUSTAL
}
\usage{
read.akterm(file)
}
\arguments{
  \item{file}{
    name (and path if necessary) of the file to read}
}
\value{
  returns a data.frame with the following columns:
  \enumerate{
    \item "Time" is the time of obeservation as POSIXct.
    \item "STA" is the WMO number of the observing statin
    \item "QDD" is the quality byte for DD
    \item "QFF" is the quality byte for FF
    \item "DD"  is the wind direction clockwise in degrees from north
    \item "FF"  is the wind speed in m/s
    \item "QQ1" is the quality byte for KM
    \item "KM"  is the Klug/Manier stability class
    \item "QQ2" is the quality byte for HM
    \item "HM"  is unsused
    \item "QQ3" is the overall quality byte.
  }
  Note that DD and FF are always in degrees and m/s, indepenent from
  the value of the quality byte. Values where the quality byte is 9
  are always NA. 
  I.e. the quality byte is given for informational purposes only.
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
series <- read.akterm(system.file("extdata", "example-series.akterm", package="readmet"))
plot(series$Time, series$FF, type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
\keyword{ file }
\keyword{ DMNA }
\keyword{ AUSTAL }
