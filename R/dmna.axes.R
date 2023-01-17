dmna.axes <- function(file) {
  header <- dmna.header(file)
  #
  # get axes length
  #
  if ( ! ( "dims" %in% names(header))) {
    stop (paste(file,"does not contain number of dimensions"))
  }
  dims=header$dims
  if ( ! ( "axes" %in% names(header))) {
   if ( dims >= 1 & dims <= 3 ) {
     header$axes='xyz'[1:dims]
   }else{
     stop (paste(file,"does not contain names of axes"))
   } 
  }
  if ( ! ( "hghb" %in% names(header) &&  "lowb" %in% names(header))) {
    stop (paste(file,"does not contain index value ranges"))
  }
  lowb=read.table(header=F,text=header$lowb)
  hghb=read.table(header=F,text=header$hghb)

  xlen=hghb[,1]-lowb[,1]+1
  if (dims>1) {
    ylen=hghb[,2]-lowb[,2]+1
  } else {
    ylen=1
  }
  if (dims>2) {
    zlen=hghb[,3]-lowb[,3]+1
  } else {
    zlen=1
  }
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
  x=xmin+delta*((1:xlen)-1)
  if (dims>1) {
    y=ymin+delta*((1:ylen)-1)
  } else {
    y=0:0
  }
  if (zlen==1) {
    return(list(x=x,y=y))
  } else {
    if ( ! (  "sk" %in% names(header))) {
      stop (paste(file,"does not contain level heights"))
    }
    sk=as.matrix(read.table(header=F,text=header$sk))
    if (! ( hghb[,3] <= length(sk) && lowb[,3] >= 1 ) ) {
      stop (paste(file,"level indices outside givel level heights"))
    }
    z=sk[lowb[,3]:hghb[,3]]
    return(list(x=x,y=y,z=z))
  }
}
