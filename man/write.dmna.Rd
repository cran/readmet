\name{write.dmna}
\alias{write.dmna}
\title{
  function to write data to a file in DMNA-format
}
\description{
  writes gridded on-e, two or three-dimensional data or timeseries to a DMNA format
  that is used for example by the German national dispersion model AUSTAL
}
\usage{
write.dmna(filename, values, axes=NULL, name=NULL, types=NULL, vldf="V", debug=FALSE)
}
\arguments{
  \item{filename}{
    name (and path if necessary) of the file to write}
  \item{values}{
    matrix or list(matrix) or data.frame. Matrix or list(matrix) implies writing gridded data.
    list(matrix) must be named using the variable names.
    data.frame implies writing a timeseries and must contain a column named "te" containing time
    as POXIXct.}
  \item{axes}{data.frame. Required for gridded data. Must contain columns "x" and column "y" if
    values are two-dimensional and additionally "sk" (or "z") if data are three-dimensional. 
    The spacing of all values in "x" and "y" must be identical.}
  \item{name}{character string. Variable name. Required if values is of class matrix.
    Is ignored else.}
  \item{types}{named list. Variable type for each variable. If `values` is a list, 
    names in `types` must match names in in `values`. If `values` is a matrix, a
    name in `types` must match `name`. Types are:
    \itemize{
      \item{`d`: integer number format}
      \item{`f`: floating point format (suitable for numbers 0.1 ... 99999.)}
      \item{`e`: exponential format}
      \item{`t`: timestamp}
    }}
  \item{vldf}{character string. specifies, where values are located in the model grid.
    "V" denotes volume average, "P" volumne center point values, or "" volume edge values.}
  \item{debug}{boolean. ignored. Is kept for compatibility}
    

}
\value{
  returns nothing
}
\references{
AUSTAL 3.1 model reference (by Ingenieurbuero Janicke)
}
\author{
  Clemens Druee
}
\note{
  Introduced in version 1.6.0.
}
\examples{
#create data and write file
len <- 25
val <- list(random = matrix(runif(len * len),nrow = len))
ax <- list(x = 1:len,
           y = 1:len)
write.dmna("myfile.dmna", val, ax)

# show first lines
writeLines(readLines('myfile.dmna', n=12))

# delete file
unlink("myfile.dmna")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
\keyword{ file }
\keyword{ DMNA }
\keyword{ AUSTAL }
