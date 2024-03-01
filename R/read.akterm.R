read.akterm <- function(file) {
  # Define column names based on whether precipitation data is needed
  akt_columns <- c('KENN', 'STA', 'JAHR', 'MON', 'TAG', 'STUN', 'NULL', 'QDD', 'QFF', 'DD', 'FF', 'QQ1', 'KM', 'QQ2', 'HM', 'QQ3')
  prec_columns <- c('PP', 'QPP')
  
  # Read file, skip comment lines
  lines <- readLines(file)
  # Skipping header and anemometer heights lines
  data_lines <- lines[!grepl("^\\*|^\\+", lines)]
  
  # Split each line by spaces and convert to a matrix
  data <- read.table(text = data_lines)
  if (ncol(data) == length(akt_columns)){
    names(data) = akt_columns
    prec = FALSE
  } else if (ncol(data) == (length(akt_columns) + length(prec_columns))){
    names(data) = c(akt_columns, prec_columns)
    prec = TRUE
  } else {
    stop ("this file does not have the correct number of columns")
  }
  
  # Apply quality flags for wind speed and direction
  data$FF[data$QFF == 0] <- data$FF[data$QFF == 0] * 0.514
  data$FF[data$QFF == 1] <- data$FF[data$QFF == 1] * 0.1
  data$FF[data$QFF == 9] <- NA  # NA for missing data
  data$DD[data$QDD == 0] <- data$DD[data$QDD == 0] * 10
  data$DD[data$QDD == 9] <- NA  # NA for missing data
  if (prec){
    data$PP[data$QPP == 9] <- NA  # NA for missing data
  }
  # Combine date and time into a single POSIXct datetime object
  data$Time <- as.POSIXct(with(data, paste(JAHR, MON, TAG, STUN)), format = "%Y %m %d %H")
  
  # Drop columns that are not needed for the analysis
  data <- data[, !(names(data) %in% c('KENN', 'JAHR', 'MON', 'TAG', 'STUN', 'NULL'))]
  
  return(data)
}