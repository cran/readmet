write.arcgrid <- function(z,file,xlen,ylen,xll,yll,delta,grid,naval=-9999) {
  if ( ( missing(xlen) | missing(ylen) |
         missing(xll) | missing(yll) |
         missing(delta) ) & 
       missing(grid) ) {
    stop (paste('incomplete grid definition: ',
                'either give xlen,ylen,xll,yll, and delta',
                'or grid'))
  }
  if ( ( ! missing(xlen) | ! missing(ylen) |
         ! missing(xll) | ! missing(yll) |
         ! missing(delta) ) & 
       ! missing(grid) ) {
    stop (paste('duplitcate grid definition: ',
                'either give xlen,ylen,xll,yll, and delta',
                'or grid'))
  }
  # move values into "grid" list
  if(missing(grid)) {
    grid["xlen"] <- xlen
    grid["ylen"] <- ylen
    grid["xll"] <- xll
    grid["yll"] <- yll
    grid["delta"] <- delta
  } 
  grid["na"] <- naval
  # replaxe NA by naval
  z[is.na(z)]=naval
  z[is.infinite(z)]=naval
  # rotate 90deg left
  z <- apply(t(z), 2, rev)
  #
  # write file
  #
  f <- file(file,"w+t")
  write(sprintf("%-14s%-64d","ncols",       grid["xlen"]), file=f,append=TRUE)
  write(sprintf("%-14s%-64d","nrows",       grid["ylen"]), file=f,append=TRUE)
  write(sprintf("%-14s%-64f","xllcorner",   grid["xll"]),  file=f,append=TRUE)
  write(sprintf("%-14s%-64f","yllcorner",   grid["yll"]),  file=f,append=TRUE)
  write(sprintf("%-14s%-64f","cellsize",    grid["delta"]),file=f,append=TRUE)
  write(sprintf("%-14s%-64f","NODATA_value",grid["delta"]),file=f,append=TRUE)
  #
  write.table(z,file=f,append=TRUE,row.names = FALSE, col.names = FALSE)
  close(f)
}