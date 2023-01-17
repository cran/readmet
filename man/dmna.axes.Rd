\name{dmna.axes}
\alias{dmna.axes}
\title{
  read grid positions from a file in DMNA-format
}
\description{
  Files in DMNA format contain gridded data. The header
  provides information on how the data are organized.
  This routine expands the positions of the grid planes along
  each dimension used in the file. The values are put out in
  columns x, y, z for each respective axis.
}
\usage{
dmna.axes(file)
}
\arguments{
  \item{file}{
    name (and path if necessary) of the file to read
}
}
\value{
  returns a data.frame with the following columns:
  \item{x}{grid positions along the x axis}
  \item{y}{grid positions along the y axis}
  \item{z}{grid level heights}
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
  dimension.  The number of dimensions is currently 
  limited to three.
}

\examples{
ax <- dmna.axes(system.file("extdata", "example.dmna", package="readmet"))
ax$x
}

\keyword{ IO }
\keyword{ file }
\keyword{ DMNA }
\keyword{ AUSTAL }
