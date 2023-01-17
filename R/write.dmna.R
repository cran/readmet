write.dmna <- function( filename, values, axes=NULL, name=NULL,
                        types=NULL, vldf="V", debug=FALSE) {
  
  known_types = c('c', 'd', 'x', 'f', 'e', 't')
  implemented_types = c('d', 'f', 'e', 't')
  known_keys = c("cset", "prgm", "artp", "axes", "idnt", 
                 "t1", "t2", "dt", "dtbnummax", "index", 
                 "groups", "xmin", "ymin", "delta", "refx", 
                 "refy", "ggcs", "zscl", "sscl", "sk", 
                 "name", "unit", "vldf", "valid", "locl", 
                 "form", "refv", "exceed", "sequ", "file", 
                 "dims", "size", "lowb", "hghb")
  if (data.class(values) == "matrix"){
    #
    # cast matrix to single-element list
    # variable name
    #
    if (is.null(name)){
      stop("name must be given if values is a matrix")
    }
    values = list(values)
    names(values) = c(name)
  }
  #
  # apply variable types names given separately 
  #
  vartype = list()
  if (! is.null(types)){
    if (length(setdiff(types,values)) != 0 ){
      stop('names of list types must match names of list values')
    } else {
      for (var in names(values)){
        if (! types[[var]] %in% known_types){
          stop(sprintf('unkown type "%s" for variable: %s', types[[var]], var))
        } else if (! types[[var]] %in% implemented_types){
          stop(sprintf('type not implemented: %s', types[[var]]))
        } else {
          vartype[[var]] = types[[var]]  
        }
      }
    }
  } 
  #
  # type-specific processing
  #
  if (data.class(values) == "list"){
    filetype='grid'
    #
    #  break down axis values
    #
    if (is.null(axes)){
      stop('no axes values given')
    }
    if (! 'y' %in% names(axes)) {
      delta = unique(diff(axes$x))
      xmin=min(axes$x)
    } else {
      delta = unique(c(diff(axes$x),diff(axes$y)))
      xmin=min(axes$x)
      ymin=min(axes$y)
    }
    if (length(delta) > 1){
      stop('horizontal grid spacing not unique')
    } else {
      delta = delta[1]
    }
    if ('sk' %in% names(axes)) {
      sk = axes$sk
    } else if ('z' %in% names(axes)) {
      sk = axes$z
    } else {
      sk = NULL
    }
    #
    # cast all matrices to three dimensions
    #
    alldim = NULL
    vartype = list()
    for (var in names(values)){
      #
      # check that all matrices have same dims
      #
      if (is.null(alldim)){
        alldim = dim(values[[var]])
        ndim = length(alldim)
      } else {
        if (dim(values) != alldim){
          stop(sprintf('variable does not have identical dimension: %s', var))
        }
      }
      #
      # raise number of dims to three
      #
      if ( ndim == 1) {
        dim(values[[var]]) = c(alldim, 1, 1)
      } else if ( ndim == 2) {
        dim(values[[var]]) = c(alldim, 1)
      } else if ( ndim != 3) {
        stop(sprintf('illegal number of dimensions: %d', ndim))
      }
      #
      #  remember dimensions
      #
      filedims = dim(values[[var]])    
    }
  } else if (data.class(values) == "data.frame") {
    filetype='timeseries'
    values$te=strftime(values$te, "  %Y-%m-%d.%H:%M:%S", tz="")
  } else {
    stop(sprintf('dont know how to handle value class: %s', data.class(values)))
  }    
  #
  # determine variable format and range
  #
  form = list()
  size = -1
  filefmt = c()
  for (var in names(values)){
    #
    # auto-determine variable type and field width
    #
    if (! var %in% names(vartype)){
      if (var == "te") {
        vartype[[var]] = "t"
      } else if (filetype == 'timeseries' & grepl("\\.",var)){
        # prefer exp form for source strenghts in timeseries
        vartype[[var]] = "e"
      } else if (all((values[[var]] - floor(values[[var]])) == 0)) {
        vartype[[var]] = "d"
      } else {
        digits = max(ceiling(log10(abs(values[[var]]))))
        if (digits > 7 | digits < 0){
          vartype[[var]] = "e"
        } else {
          vartype[[var]] = "f"
        }
      }  
    }
    #
    # determine variable format
    #
    if (vartype[[var]] == "d"){
      digits = max(ceiling(log10(abs(values[[var]]))))
      digits = max(digits, 4)
      fmt = sprintf('%%%dhd',digits)
      flen = digits
    } else if (vartype[[var]] == "f"){
      digits = max(ceiling(log10(abs(values[[var]]))))
      if (all(values[[var]] - floor(values[[var]]) == 0)){
        precision = 0 
        digits = max(digits, 5)
      } else {
        precision = 1
        digits = max(digits+2, 7)
      }
      fmt = sprintf('%%%d.%df',digits,precision)
      flen = digits
    } else if (vartype[[var]] == "e"){
      fmt = "%10.3e"
      flen = 10
    } else if (vartype[[var]] == "t"){
      fmt = "%20lt"
      flen = 0
    } else {
      stop(sprintf('wrong type for matrix: %s', vartype[[var]]))
    }
    form[[var]] = paste(tolower(var), fmt, sep="")
    filefmt = append(filefmt, fmt)
    size = size + 1 + flen
  }
  filefmt = paste(filefmt, collapse = " ")
  filefmt = gsub("l([fe])", "\\1", filefmt)
  filefmt = gsub("h([de])", "\\1", filefmt)
  
  #
  # assemble header info
  #
  header = list(file=gsub(".akterm$","",filename))
  header[['mode']] = "text"
  header[['cset']] = "UTF-8"
  header[['form']] = c(as.vector(unlist(form),'character'))
  if (filetype == 'grid'){
    header[['dims']] = 3
    header[['lowb']] = c(rep(1, 3))
    header[['hghb']] = c(filedims)
    header[['sequ']] = "k+,j-,i+"
    header[['size']] = size
    header[['xmin']] = xmin
    header[['ymin']] = ymin
    header[['delta']] = delta
    
  } else if (filetype == 'timeseries') {
    header[['dims']] = 1
    header[['lowb']] = 1
    header[['hghb']] = nrow(values)
    header[['sequ']] = "i+"
    header[['artp']] = "ZA"
  }
  #
  #  write file
  #
  # write header
  con = file(filename,"w")
  for (key in known_keys){
    if (key %in% names(header)){
      if (data.class(header[[key]]) == 'numeric'){
        # numbers without quotes
        line = paste(key,paste(header[[key]], collapse = ' '), sep='  ')
      } else {
        # characters surrounded by quotes
        line = paste(key,paste(shQuote(header[[key]], type='cmd'), collapse = ' '), sep=' ')
      }
      writeLines(line, con, sep = "\r\n")
    }
  }
  writeLines("*", con, sep = "\r\n")
  close(con)
  #
  # write data body (type specific)
  if (filetype == 'grid'){ 
    con = file(filename,"a")
    #
    # write block for each layer (3. dim)
    for (k in 1:filedims[3]){
      #
      # block separator
      if (k > 1){
        writeLines("*", con, sep = "\r\n")
      }
      #
      # write lines for each y grid line (2. dim)
      lines = c()
      for (j in rev(1:filedims[2])) {
        #
        # write group for each x grid line (1. dim)
        groups = c()
        for (i in 1:filedims[1]){
          groupvals = c()
          #
          # write sequence of all variables in each group
          for (var in names(values)){
            groupvals = append(groupvals, values[[var]][i,j,k])
          }
          groupstr = sprintf(filefmt, groupvals)
          groups = append(groups, groupstr)
        }
        line = paste(groups, collapse = " ")
        lines = c(lines, line)
      }
      writeLines(lines, con, sep = "\r\n")
    }
    close(con)
  } else if (filetype == 'timeseries') {
    out = values
    for (var in names(values)) {
      if (var != "te"){
        # get the format part of "form" header
        fmt = gsub(".*%", "%", form[[var]])
        fmt = gsub("l([fe])", "\\1", fmt)
        fmt = gsub("h([de])", "\\1", fmt)
        out[[var]] = unlist(lapply(as.array(values[[var]]), sprintf, fmt=fmt))
      }
    }
    write.table(out, filename, append = TRUE,
                col.names = FALSE, row.names = FALSE,
                sep = " ", eol = "\r\n", quote = FALSE)
  }
  #
  # write footer
  con = file(filename,"a")
  writeLines(c("","***"), con, sep = "\r\n")
  close(con)
}  
