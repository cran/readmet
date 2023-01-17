dmna.header <- function(file) {
  #
  # read the file as text lines and find the divider line "*"
  #
  lines <- readLines(file)
  divider=which(lines == "*")
  if (length(divider) != 1) {
    stop (paste(file,"is not in DMNA format"))
  }
  #
  # convert the file header into named list
  #
  head_lines=lines[1:(divider-1)]
  # remove empty lines
  head_lines=head_lines[head_lines!=""]
  # convert space behind line tag into tab (if not already present)
  head_lines=sub("\ +","\t",head_lines)
  # 1st field is name 2nd and on is content
  head_names=sub("\t.*","",head_lines)
  head_values=sub("[a-z0-9]+\t","",head_lines)
  #remove tabs and quotes
  head_values=gsub("\t"," ",head_values)
  head_values=gsub("\\\"","",head_values)
  #make named list
  names(head_values)=head_names
  head=as.list( head_values )
  head$lines=divider
  return(head)
}
