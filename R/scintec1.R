#
# read header "header"
#
scintec1.header <- function(file) {
  #
  # read the file as text lines and check magic
  #
  header=list()
  lines <- readLines(file)
  if (lines[1] == "FORMAT-1") {
    header["version"]="1.0"
    header["starttime"]=as.POSIXct(substr(lines[2],1,19))
    header["filecount"]=as.numeric(substr(lines[2],21,99))
    header["instrument"]=trimws(lines[3])
    header[c("commentlines","variables","heightlevels")]=unlist(lapply(strsplit(lines[4],"\ +"),as.numeric))
    header["fixedlines"]=5
    header["datatype"]="Main Data"
  } else if (lines[1] == "FORMAT-1.1") {
    header["version"]="1.1"
    header["starttime"]=as.POSIXct(substr(lines[2],1,20))
    header["instrument"]=trimws(lines[3])
    header[c("commentlines","variables")]=unlist(lapply(strsplit(lines[4],"\ +"),as.numeric))
    header["fixedlines"]=5
    typeline=header$fixedlines+1+header$commentlines
    header["datatype"]=trimws(lines[typeline])
  } else {
      stop (paste(file,"is not in Scintec format-1.x format"))
  }
  return(header)
}  
#
# read the comments
#
scintec1.comments <- function(file,header=list()) {
  if (length(header)==0){
    header=scintec1.header(file)
  }
  comments=list()
  firstline=header$fixedlines+1
  lastline=header$fixedlines+1+header$commentlines-1
  #
  lines <- readLines(file)
  lines=lines[!grepl("^\ *#",lines)]
  #
  for (pointer in firstline:lastline){
    # ignore commented lines
    if ( grepl("^\ *#",lines[pointer]) == FALSE ) {
      fields = unlist(strsplit(lines[pointer],":"))
      if ( length(fields) != 2 ){
        warning(paste("cannot comprehend comment line: ",lines[pointer]))
      }
      key=trimws(fields[1])
      value=trimws(fields[2])
      comments[key]=value
    }
  }
  return(comments)
}  
#
# read variable definitions
#
scintec1.variables <- function(file,header=list()) {
  #
  # get header and read the file as text lines
  #
  if (length(header)==0){
    header=scintec1.header(file)
  }
  lines <- readLines(file)
  lines=lines[!grepl("^\ *#",lines)]
  #
  # read actual definitions
  #
  firstline=header$fixedlines+header$commentlines+2
  lastline=header$fixedlines+header$commentlines+header$variables+2
  for (pointer in firstline:lastline){
    fields = trimws(unlist(strsplit(lines[pointer],"#")))
    # fix symbol for error code:
    if ( fields[1] == "error code" & fields[2] != "error") {
      fields[3] = fields[2]
      fields[2] = "error"
      fields[6] = NA
    }
    # fix Time line:
    if ( fields[1] == "Time" & length(fields) == 5 ){
      fields[6] = NA
    }
    names(fields)=c("label","symbol","unit","type","error.mask","gap.value")
    if ( ! exists("vars") ) {
      vars=data.frame(t(fields))
    } else {
      vars=rbind(vars,data.frame(t(fields)))
    }
  }
  # remove dummy column
  vars=vars[, !(names(vars) %in% "dummy")]
  # make shure entries have proper type
  for ( c in c("label","symbol","unit","type","error.mask") ) {
    vars[[c]]=as.character(vars[[c]])
  }
  vars$gap.value=suppressWarnings(as.numeric(as.character(vars$gap.value)))
  #  
  return(vars)
}
#
# read non-profile data
#
scintec1.nonprofile <- function(file,header=list(),vars=list()) {
  #
  # get info about variables / diemensions
  #
  if (length(header)==0){
    header=scintec1.header(file)
  }
  if (length(vars)==0){
    vars=scintec1.variables(file)
  }
  #
  # read the file as text lines
  #
  lines <- readLines(file)
  lines=lines[!grepl("^\ *#",lines)]
  #
  # switch Format versions
  #
  if (header["version"]=="1.0"){
    #
    # Scintec Format-1
    #
    # if nonprofile variuable were recorded
    #
    if ( "NS" %in% vars$type ) {
      #
      # loop data
      #
      pointer=header$fixedlines+1
      continue = TRUE
      while ( continue ) {
        if ( pointer <= length(lines) ){
          # look for date/time
          if ( grepl("^....-..-.. ",lines[pointer]) ) {
            # get date/time
            field=unlist(strsplit(trimws(lines[pointer]),"\ +"))
            datetime=as.numeric(as.POSIXct(paste(field[1],field[2])))
            message(paste("   ... reading nonprofile data:",datetime))
            line=gsub("^\ *#","",lines[pointer+1])
            names=vars$symbol[vars$type == 'NS']
            values=suppressWarnings(as.numeric(unlist(strsplit(trimws(lines[pointer+2]),"\ +"))))
            v=as.list(values)
            names(v)=names
            v["time"]=datetime
            df=data.frame(v)
            if ( ! exists("npdata") ) {
              npdata=df
            } else {
              npdata=rbind(npdata,df)
              #rm(df)
            }
          }
          pointer=pointer+1
        } else {
          continue = FALSE
        }
      }
      npdata$time=as.POSIXct(npdata$time,origin="1970-01-01")
      return(npdata)
    } else {
      return(NULL)
    }
  } else if (header["version"]=="1.1"){
    #  
    # Scintec Format-1.1
    #
    pointer=header$fixedlines+1+header$commentlines+1+header$variables+1+1
    continue = TRUE
    while ( continue ) {
      if ( pointer <= length(lines) ){
        # get date/time
        field=unlist(strsplit(trimws(lines[pointer]),"[\ \t]+"))
        timestr=unlist(strsplit(field[1],"/"))
        datetime=as.numeric(as.POSIXct(timestr[2]))
        message(paste("   ... reading nonprofile data:",paste(datetime)))
        values=suppressWarnings(as.numeric(field))
        v=as.list(values)
        if (length(v) != length(vars$symbol)){
          warning( sprintf("incomplete line #%i",pointer))
        } else{
          names(v)=vars$symbol
          v[["time"]]=datetime
          df=data.frame(v)
          if ( ! exists("npdata") ) {
            npdata=df
          } else {
            npdata=rbind(npdata,df)
            #rm(df)
          }
        }
        pointer=pointer+1
      } else {
        continue = FALSE
      }
    }
    npdata$time=as.POSIXct(npdata$time,origin="1970-01-01")
    return(npdata)
    
        
  } else {
    stop(paste(c('unknown Format version',header["version"],' reading nonprofile data ')))
  }
}
#
# read profile data
#
scintec1.profile <- function(file,header=list(),vars=list()) {
  #
  # read the file as text lines and check magic
  #
  lines <- readLines(file)
  lines=lines[!grepl("^\ *#",lines)]
  #
  # get info about variables / dimensions
  #
  if (length(header)==0){
    header=scintec1.header(file)
  }
  if (length(vars)==0){
    vars=scintec1.variables(file)
  }
  #
  # switch Format versions
  #
  if (header["version"]=="1.0"){
    #
    # Scintec Format-1
    #
    #
    if ( "NS" %in% vars$type ) {
      npvars=TRUE
    } else {
      npvars=FALSE
    }
    levels=header$heightlevels
    # initialize
    pointer=6
    continue = TRUE
    # loop data
    while ( pointer <= length(lines) ){
      # look for date/time
      if ( grepl("^....-..-.. ",lines[pointer]) ) {
          # get date/time
          field=unlist(strsplit(trimws(lines[pointer]),"\ +"))
          datetime=as.POSIXct(paste(field[1],field[2]))
          message(paste("   ... reading profile data:",paste(field[1],field[2])))
          # jump over non-profile data
          if ( npvars ) {
            pointer=pointer+2
          } else {
            pointer=pointer+1
          }
          # get the lines that form the block of numbers
          block=paste(lines[pointer:(pointer+levels)],"\n")
          # read the numbers
          df=read.table(textConnection(block))
          # get names of the profile variables
          names(df)=vars$symbol[vars$type != 'NS']
          if ( ! exists("fields") ) {
            # initialize output on first pass
            fields=list()
            for ( c in names(df) ) {
              fields[[c]]=as.matrix(t(df[[c]]))
            }
            times=c(datetime)
          } else {
            # append read data to output
            for ( c in names(df) ) {
              nf=as.vector((df[[c]]))
              fields[[c]]=rbind(fields[[c]],nf)
              rm(nf)
            }
            times=c(times,datetime)
          }
          # jump to end of block
          pointer=pointer+levels-1
      } else {
        # break, if end of file is reached
        continue = FALSE
      }
      # continue search for next block
      pointer=pointer+1
    }
    # replace gap values in each field by NA
    if ( exists("fields") ) {
      for ( c in names(fields) ) {
        fields[[c]][fields[[c]]==as.numeric(vars$gap.value[which(vars$symbol==c)])]=NA
      }
    }
    # replace field of height by vector of heights
    fields[["z"]]=fields[["z"]][1,]
    # add vector of times 
    fields[["time"]]=times
  } else if (header["version"]=="1.1"){
    #  
    # Scintec Format-1.1
    #
    # (no grid variables)
    fields=list()
  } else {
    stop(paste(c('unknown Format version',header["version"],' reading nonprofile data ')))
  }
  return(fields)
}

read.scintec1 <- function(files) {
  for ( file in files ) {
    message(paste("opening file",file))
    header=scintec1.header(file)
    vars=scintec1.variables(file,header)
    if (length(vars) > 0){
      if ( length(files) == 1 | file==files[1]) {
         fields=scintec1.profile(file,header,vars)
         series=scintec1.nonprofile(file,header,vars)
       } else {
         fmore=scintec1.profile(file,header,vars)
         for ( c in names(fields) ) {
           if ( c == "z" ) {
             if ( ! all(fields$z == fmore$z )) {
               stop(paste("ERROR: level heights change in file",file))
             } 
           } else if ( c == "time" ) {
             fields[[c]]=c(fields[[c]],fmore[[c]])
           } else if (is.matrix(fields[[c]])) {
             fields[[c]]=rbind(fields[[c]],fmore[[c]])
           } else {
             stop(paste("ERROR: dont know how to handle variable",c))
           }
         }
         smore=scintec1.nonprofile(file,header,vars)
         series=rbind(series,smore)
       }
    }  
  }
  if (exists("fields")) {
    # sort data by time
    if ( is.unsorted(fields$time) ) {
      perm=order(fields$time)
      for ( c in names(fields) ) {
        if ( c == "z" ) {
          # do nothing
        } else if ( c == "time" ) {
          fields[[c]]=fields[[c]][perm]
        } else if (is.matrix(fields[[c]])) {
          fields[[c]]=fields[[c]][perm,]
        } else {
          stop(paste("ERROR: dont know how to handle variable",c))
        }
      }
      series=series[perm,]
    }
  } else {
    fielsd=list()
  }
  # add non-profile data as vectors to list
  for (c in names(series)){
    if ( ! ( c == "time" & "time" %in% names(fields))){
      fields[[c]]=series[[c]]
    }
  }
  return(fields)
}
