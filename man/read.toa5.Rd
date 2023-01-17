\name{read.toa5}
\alias{read.toa5}
\title{
  function to read contents of a file in TOA5-format
}
\description{
  reads the contents of a Campbell Scientific table-orientad ASCII format 5 (TOA5) into a data frame, columns are named after the TOA5 vribale names
}
\usage{
read.toa5(file)
}
\arguments{
  \item{file}{
    name (and path if necessary) of the file to read
  }
}
\value{
  data.frame
}
\references{
Campbell Scientific CR3000 data logger manual
}
\author{
  Clemens Druee
}
\examples{
dat <- read.toa5(system.file("extdata", "example.dat", package="readmet"))
dat$time <- as.POSIXct(dat$TIMESTAMP)
plot (dat$time, dat$AirTC_Avg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ read }
\keyword{ format }
\keyword{ TAO5 }
