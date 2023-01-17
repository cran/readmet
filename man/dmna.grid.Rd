\name{dmna.grid}
\alias{dmna.grid}
\title{
  read (horizontal) grid definition from a file in DMNA-format
}
\description{
  Files in DMNA format contain gridded data.  The header
  provides information on the grid definition.  This function extracts
  the horizontal grid definition ans provides a list that may be
  supplied directly to function write.arcgrid.  DMNA is used for example
  by the German national dispersion model AUSTAL.
}
\usage{
dmna.grid(file)
}
\arguments{
  \item{file}{
    name (and path if necessary) of the file to read
}
}
\value{
  returns a list of character. List entries are:
    \item{xlen}{number of cells in x direction}
    \item{ylen}{number of cells in y direction}
    \item{xll}{easting of the lower left (i.e. southwest) corner}
    \item{yll}{northing of the lower left (i.e. southwest) corner}
    \item{delta}{grid node spacing}
}
\references{
AUSTAL 3.1 model reference (by Ingenieurbuero Janicke)
}
\author{
  Clemens Druee
}
\note{
  This function was added in version 1.3.0.
}
\examples{
dmna.grid(system.file("extdata", "example.dmna", package="readmet"))
}
\keyword{ IO }
\keyword{ file }
\keyword{ DMNA }
\keyword{ AUSTAL }
