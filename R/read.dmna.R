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
  if (debug) print (paste("artp:",artp))
  # get index order and orientation
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
    # take second character of sequence entry,
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
    form=header$form
    form=gsub("[\t ]+", " ", form)
    forms=unlist(strsplit(form," "))
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
  if (debug) print (paste("header lines: ",header$lines))
  #
  # read the data 
  #
  if ( mode == "text" ){
    if ( dims > 1 & nval > 1 ){
      stop (paste("ERROR: number of values >1 for dimensions >1 not implemented "))
    }
    table=read.table(file,skip=header$lines,
                     header=FALSE,
                     comment.char = "*",quote=NULL,
                     as.is=TRUE,blank.lines.skip=TRUE)  
    # if data section contains in-line comments
    # remove supposed comment column and following columns
    for (i in 1:dim(table)[2]){
      if (any(table[,i]=="'")){
        table = table[,1:(i-1)]
        break
      }
    }
    # warp table
    numbers=as.matrix(table)
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
  # data in the file are in FORTRAN order i.e. last index is counting fastest
  # in the array command, the first index is the fastest counting
  # => use reverse order of indexes 
  value=list()
  for ( l in 1:nval ){
    value[[l]]=array(values[l+seq(0,(nval*numrec-1),nval)],rev(seqlen))
  }
  #
  # rearrange axes to the R-convention (equivalent to sequence "i+,j+,k+")
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
      } 
    }
    #
    # try harder to output numbers not factors 
    #
    for (i in names(out)){
      # exclude time column
      if (i!="te") {
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
  } else if (artp=="CZ") {
    #
    # artp CZ is a timeseries though described as as 2D array (why?)
    #
    if (debug) print ('filetype CZ: adding time column')
    t1 = .parsedifftime(header$t1)  # "00:00:00"
    t2 = .parsedifftime(header$t2)  # "366.00:00:00"
    dt = .parsedifftime(header$dt)  # "01:00:00"
    rd = gsub("([0-9]{4}-[0-9]{2}-[0-9]{2})[ T.]([0-9]{2}:[0-9]{2}:[0-9]{2}).*",
              "\\1T\\2",
              header$rdat) # "2000-01-01T00:00:00+0100" or "2000-01-01.00:00:00+0100" or "2000-01-01 00:00:00"
    rdat = as.POSIXct(rd, format="%Y-%m-%dT%H:%M:%S") 
    if (debug) print (paste(rdat + t1, rdat + t2 - dt, dt))
    te = seq(from = rdat + t1, to = rdat + t2 - dt, by = dt)
    if (debug) print (paste('... ',strftime(te[1],"%F %T"),' -- ',strftime(te[length(te)],"%F %T")))
    out = cbind(as.data.frame(dat), data.frame(te=te))
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
#
#' @NoRd
.parsedifftime <- function(s){
  # parse difftime format [ddd.]hh:mm:ss
  #
  # split into days and hh:mm:ss
  if (grepl(".", s, fixed=TRUE)) {
    # if days are present, use them
    p <-  unlist(strsplit(s, ".", fixed=TRUE))
    d <- p[1]
    s <- p[2]
  } else {
    # if not, use zero days
    d <- "0" 
  }
  re1 <- as.difftime(as.numeric(d), units = "days")
  re2 <- as.difftime(s, format="%H:%M:%S")
  return(re1 + re2)
}