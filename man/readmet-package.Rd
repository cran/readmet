\name{readmet}
\alias{readmet-package}
\alias{readmet}
\docType{package}
\title{
Read some less Popular Formats Used in Meteorology
}
\description{
Contains tools for reading and writing data from or to files in dmna, Scintec Format-1, Campbell Scientific TOA5 formats.
}
\details{
\tabular{ll}{
Package: \tab readmet\cr
Type: \tab Package\cr
Version: \tab 1.6.9\cr
Date: \tab 2023-02-11\cr
License: \tab GPL-3\cr
}
Currently supported:\cr
\tabular{ll}{
DMNA \tab Self describing gridded ASCII by Ing.Buero Janicke/Umweltbundesamt\cr
     \tab format descrition in e.g. AUSTAL dispersion model description https://www.umweltbundesamt.de/sites/default/files/medien/2338/dokumente/austal_de.pdf (available only in German) \cr
Scintec Format-1 \tab Format used to store atmospheric data from wind profilers, SODAR, RASS and Scintillometer devices by Scintec AG, Rottenburg, Germany \cr
     \tab format descrition in the device manuals, available to customers upon registration from https://www.scintec.com \cr
TOA5 \tab Table oriented ACSII format #5 by Campbell Scientific \cr
     \tab format descrition in the device manuals, e.g. CR3000 data logger https://s.campbellsci.com/documents/us/manuals/cr3000.pdf) \cr
}

Available functions:\cr
\tabular{ll}{
Function \tab Format\cr
\cr
read.akterm \tab function to read contents of a file in AKTERM-format\cr
read.dmna \tab function to read contents of a file in DMNA-format\cr
dmna.axes \tab read grid positions from a file in DMNA-format \cr
dmna.header \tab read header information from a file in DMNA-format\cr
dmna.grid \tab read (horizontal) grid definition from a file in DMNA-format\cr
write.dmna \tab function to write data to a file in DMNA-format\cr
write.arcgrid \tab function to write data (read from DMNA file) into an ESRI ArcInfo gridded ASCII file\cr
\cr
read.toa5 \tab function to read contents of a file in TOA5-format\cr
\cr
read.scintec1 \tab function to read the contents of a file in Scintec Format-1\cr
scintec1.comments \tab read the comment fields of a file in Scintec Format-1 format\cr
scintec1.header \tab read the header of a file in Scintec Format-1 format\cr
scintec1.nonprofiles \tab read the single-level data from a file in Scintec Format-1 format\cr
scintec1.profile \tab read the multi-level (profile) variables from a Scintec Format-1 file\cr
scintec1.variables \tab read the variable descriptions from a file in Scintec Format-1 format\cr
}
}
\author{
Clemens Druee, Umweltmeteorologie, Universitaet Trier, Germany <druee@uni-trier.de>
}
\keyword{ package }
\keyword{ IO }
\keyword{ file }
\keyword{ DMNA }
\keyword{ TOA5 }
\keyword{ Scintec }
