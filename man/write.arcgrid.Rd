\name{write.arcgrid}
\alias{write.arcgrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to write data (read from DMNA file) into an ESRI ArcInfo gridded ASCII file
}
\description{
This function writes a 2D matrix into an ESRI ArcInfo gridded ASCII file that
can be easily imported into most geographic information systems (GIS)
}
\usage{
write.arcgrid(z,file,xlen,ylen,xll,yll,delta,grid,naval=-9999)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
    2D matrix containing tha data
  }
  \item{file}{
    name (incuding path, if needed) of the file to write to
  }
  \item{xlen,ylen}{
    number of data alon x and y axis, resp.
  }
  \item{xll,yll}{
    position of the lower left (i.e. southwest) corner
  }
  \item{delta}{
    grid spacing
  }
  \item{grid}{a list containing grid parameters. 
  Instead of providing xlen,ylen,xll,yll, and delta individually, 
  containing these values way be provided. For example, such a list is returned
  when calling \code{\link{dmna.grid}} }
  \item{naval}{
    value written to file instead of NA and +/-Inf
  }
}
\details{
The standard plotting functions for R plot columns along the x axis and
rows along the y axis. Hence the matrix is rotated 90 degrees left (compared to 
\code{\link{write.table}} or \code{\link{write.csv}}) to yield the same
orientation of the data when plotted in R and in GIS.}
\value{
nothing.
}
\references{
  Wikipedia entry on Esri grid : \url{https://en.wikipedia.org/wiki/Esri_grid}
}
\author{
  Clemens Druee
}
\note{
  This function was added in version 1.3.0.
}
\seealso{
  \code{\link{read.dmna}}, \code{\link{dmna.grid}}
}
\examples{
# read data and grid info
infile <- system.file("extdata", "example.dmna", package="readmet")
so2 <- read.dmna(infile)
grid <- dmna.grid(infile)
# write file
write.arcgrid(so2, file = "myfile.grid", grid = grid)
# show head of file
writeLines(readLines("myfile.grid", n=7))
# delete file
file.remove("myfile.grid")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
\keyword{file}
\keyword{ESRI}
\keyword{ArcInfo}
