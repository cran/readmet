dmna.grid <- function(file) {
  header <- dmna.header(file)
  #
  # get axes length
  #
  if ( ! ( "dims" %in% names(header))) {
    stop (paste(file,"does not contain number of dimensions"))
  }
  dims=header$dims
  if ( header$artp == "ZA" | dims<2 ) {
    stop (paste("this function does not apply to timeseries"))
  }
  if ( ! ( "hghb" %in% names(header) 
         &&  "lowb" %in% names(header))) {
    stop (paste(file,"does not contain index value ranges"))
  }
  lowb=read.table(header=F,text=header$lowb)
  hghb=read.table(header=F,text=header$hghb)
  xlen=hghb[,1]-lowb[,1]+1
  ylen=hghb[,2]-lowb[,2]+1
  #
  # get axes start
  #
  if ( ! (  "xmin" %in% names(header) 
         && "ymin" %in% names(header) 
         && "delta" %in% names(header))) {
    stop (paste(file,"does not contain all information on axes"))
  }
  xmin=as.numeric(header$xmin)
  ymin=as.numeric(header$ymin)
  delta=as.numeric(header$delta)
  #
  # reference position
  #
  if ( ! (  "refx" %in% names(header) 
            && "refy" %in% names(header))) {
    stop (paste(file,"does not contain all information on axes"))
  }
  refx=as.numeric(header$refx)
  refy=as.numeric(header$refy)
  #
  # lower left corner reference:
  #
  xll=refx+xmin
  yll=refy+ymin
  #
  # return list
  #
  out=c(xlen,ylen,xll,yll,delta)
  names(out)=c("xlen","ylen","xll","yll","delta")
  return(out)
}

