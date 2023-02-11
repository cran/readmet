\name{dmna.header}
\alias{dmna.header}
\title{
  read header information from a file in DMNA-format
}
\description{
  Files in DMNA format contain gridded data. The header
  provides information on how the data are organized
  and other user-specified meta data. DMNA is used for example
  by the German national dispersion model AUSTAL2000.
}
\usage{
dmna.header(file)
}
\arguments{
  \item{file}{
    name (and path if necessary) of the file to read
}
}
\value{
  list(character)
}
\references{
AUSTAL2000 model reference (by Ingenieurbuero Janicke)
}
\author{
  Clemens Druee
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
hdr <- dmna.header(system.file("extdata", "example-grid.dmna", package="readmet"))
hdr$idnt
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
\keyword{ file }
\keyword{ DMNA }
\keyword{ AUSTAL2000 }
