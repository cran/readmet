\name{read.scintec1}
\alias{read.scintec1}
\title{
  function to read the contents of a file in Scintec Format-1
}
\description{
  reads output data (extension ".mnd") from 
  Scintec atmospheric measurement systems into a list that
  contains the variable names as names and the values as 
  either matrix or vector.
}
\usage{
read.scintec1(files)
}
\arguments{
  \item{files}{
    array containing names (and path if necessary) of the files to read
}
}
\value{
  list(different types)
  \item{vector}{lenght = n (number of times contained in the file)\cr
                each vevtor contains a non-profile variable (see variable list)\cr
                special varaible name "time" contains the times as POSIXct}
  \item{matrix}{dim(n,m) where n is number of times contained in the file and m is 
                the number of levels\cr
                i.e. rows represents timeseries of at a fixed level,
                columns represent profiles at a fixed time}
}
\references{
Scintec APRun software manual
}
\author{
  Clemens Druee
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
mnd <- read.scintec1(system.file("extdata", "example.mnd", package="readmet"))
mnd$CT2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
\keyword{ IO }
\keyword{ Scintec }
