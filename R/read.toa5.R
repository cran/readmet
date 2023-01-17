read.toa5 <- function( file ) {
  line=read.csv(file,header=F,nrows=1)
  if (line[,1] != "TOA5") {
    stop (paste(file,"does not have TOA5 format"))
  }
  data=read.csv(file,skip=4,header=F)
  names(data)=as.matrix(read.csv(file,skip=1,header=F,nrows=1))
  if ("TIMESTAMP" %in% names(data)) {
    data$TIMESTAMP = as.POSIXct(data$TIMESTAMP)
  }
  data
}
