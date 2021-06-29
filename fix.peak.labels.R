fix.peak.labels <- function(arg1) {
  #takes a vector of strings of metabolite/peak names (that were altered by R's makenames())
  #and fixes names to match the format in data sheets
  ##George Crowley
  out <- substr(arg1, 2, length(arg1))
  out <- sub(".z", "/z", out)
  # out <- substr(out, 1, length(out)-2)
  # out <- paste0(out, "/z")
  return(out)
}