read.dmna <- function( file, val=1, debug=FALSE ) {
  #
  # read header
  #
  header <- dmna.header(file)
  #
  # get number and kind of dimensions
  #
  if ( ! ( "dims" %in% names(header))) {
    stop (paste(file,"does not contain number of dimensions"))
  }
  dims=as.integer(header$dims)
  if ( dims > 3 ) {
    stop ("this function only supports up to three dimensions")
  } else {
    if (debug) print (paste("dims:",dims))
  }
  if ("artp" %in% names(header)) {
    artp=header$artp
  } else {
    artp=""
  }
  # get index oder and orientation
  #
  # index sequence gives order (slowest counting to fastest counting)
  # of numbers in file e.g. "k+,j-,i+"
  # index position is position of axis in list seq
  # e.g. x-axis boundaries are in first column in lowb/highb 
  #      x-index "i" is found in last position, direction is + 
  #             -> fastest counting, increasing 
  #             -> along data rows, lowes x left highest x right
  #             
  if ( ! ( "sequ" %in% names(header) ) ) {
    stop (paste(file,"does not contain information on index order"))
  }
  sequ=unlist(strsplit(header$sequ,","))
  if (length(sequ) != dims) {
    stop (paste("number of indices does not match number of dimensions in",file))
  }
  if (dims==1||dims==2||dims==3) {
    # index names
    inam=c("i","j","k")
    # direction of each index in sequence
    # take scond character of sequence entry,
    # assume "+" if 2nd character is missing
    seqdir=substr(paste(sequ[1:dims],"+",sep=""),2,2) 
    seqind=substr(sequ[1:dims],1,1) 
    # position of each index in sequence
    ipos=rep(NA,dims)
    idir=rep("",dims)
    for (i in 1:dims){
      if (inam[i] %in% seqind){
        ipos[i]=grep(inam[i],seqind)
        idir[i]=seqdir[ipos[i]]
      }  
    }
  } else {
    stop (paste(dims," dimensions are not supported by this version"))
  }
  #
  # index boundaries
  #
  if ( ! ( "hghb" %in% names(header) &&  "lowb" %in% names(header))) {
    stop (paste(file,"does not contain information on grid size"))
  }
  lowb=as.integer(read.table(header=F,text=header$lowb))
  hghb=as.integer(read.table(header=F,text=header$hghb))
  ilen=hghb-lowb+1
  seqlen=rep(NA,dims)
  seqlen[ipos]=ilen
  numrec=prod(ilen)
  if (debug) {
    print(paste('ipos',ipos))
    print(paste('idir',idir))
    print(paste('ilen',ilen))
    print(paste('seqdir',seqdir))
    print(paste('seqlen',seqlen))
  }
  #
  # ascii or binary ?
  #
  if ( "mode" %in% names(header)) {
    mode=header["mode"]
  } else {
    mode="text"
  }
  if (debug) print(paste('mode:',mode))
  #
  # how many values per data record
  #
  if ( "form" %in% names(header)) {
    forms=unlist(strsplit(header$form," "))
    nval=length(forms)
    valnams=unlist(lapply(strsplit(forms,"%"), function(l) l[[1]]))
  } else {
    nval=1
  }
  if (debug) print (paste("nval:",nval))
  #
  # select value
  #
  if (nval>1) {
    if (is.character(val)){
      if (val %in% valnams) {
        ov=which(valnams==val)
      } else {
        stop(paste('ERROR: value name ',val,' nof in file ',file))
      }
    } else {
      if (val %in% 1:nval) {
        ov=val
      } else {
        stop(paste('ERROR: value number ',val,' not in file ',file))
      }
    }
    if (length(ov)>1) {
      stop(paste('ERROR: value specification ',val,' not unique in file ',file))
    }
  } else {
    if (val!=1){
      warning(paste('value specification ',val,' not needed in file ',file))
    }
    ov=1
  }  
  #
  # read the data 
  #
  if ( mode == "text" ){
    if ( dims > 1 & nval > 1 ){
      stop (paste("ERROR: number of formats >1 for dimensions >1 not implemented "))
    }
    numbers=as.matrix(read.table(file,skip=header$lines,
                                 header=FALSE,comment.char = "*",
                                 as.is=TRUE,blank.lines.skip=TRUE))
    values=c(t(numbers))
  }else if ( mode == "binary" ){
    binfile=gsub("dmna","dmnb",file)
    to.read <- file(binfile, "rb")
    values=readBin(to.read, "numeric", n=nval*numrec, 
                   size = 4, endian = "little")
    close(to.read)
  }else{
    stop (paste("ERROR: unsopported mode ",mode))
  }
  # data in the file are in FORTRAN order i.e. last indext is counting fastest
  # in the array command, the first index is the fastest counting
  # => use reverse order of indixes 
  value=list()
  for ( l in 1:nval ){
    value[[l]]=array(values[l+seq(0,(nval*numrec-1),nval)],rev(seqlen))
  }
  #
  # rearrange axes to the R-convetion (equivalent to sequence "i+,j+,k+")
  # i.e. un-reverse order of of indices
  #
  dat=list()
  for ( l in 1:nval ){
    dat[[l]]=aperm(value[[l]],rev(ipos))
  }
  #
  # reverse order of values if an index was counting backwards
  #
  for ( l in 1:nval ){
    dat[[l]]=.revdim(dat[[l]],which(idir=='-'))
  }
  #
  # make output
  #
  if (dims==1) {
    out=as.data.frame(do.call(cbind, dat))
    names(out)=valnams
    #
    # if timeseries: find time column and convert to POSIXct
    #
    for (i in names(out)){
      if (artp=="ZA" & i=="te") {
        # "." in format string breaks processing !?
        out$te=as.POSIXct(gsub("\\."," ",out$te), "%Y-%m-%d %H:%M:%S", tz="")
      } else {
        if (is.factor(out[[i]])){
          out[[i]] = as.character(out[[i]])
        } 
        if (is.character(out[[i]])){
          if (! any(is.na(suppressWarnings(as.numeric(out[[i]]))))) {
            out[[i]] = as.numeric(out[[i]])
          }
        }
        out[[i]] = type.convert(out[[i]], as.is = FALSE)
      }
    }
    #
    # try harder to output numbers not factors
    #
    
  } else {
    out=drop(dat[[ov]])
  }
  return (out)
}
# ------------------------------------------------------------------------
#
# define a helper function following a discussion at
# http://r.789695.n4.nabble.com/distributing-the-values-of-data-frame-to-a-vector-based-on-td837896.html
#
#' @NoRd
.revdim <- function(arr, di=NA) {
  if ( is.vector(arr) ){
    return(rev(arr))
  } else {
    ndim = dim(arr)
    dims = length(ndim)
    idim = 1:dims
    rdim = idim %in% di 
    if ( sum(rdim) > 0 ) {
      rind <- function(i){
        if (! i %in% idim) {
          return(0)
        } else if (rdim[i]) {
          return(rev(1:ndim[i]))
        } else {
          return(1:ndim[i])
        }
      }
      do.call("[", c(list(arr), sapply(idim, rind, simplify=FALSE),drop=FALSE))
    } else {
      return(arr)
    }
  }
}

