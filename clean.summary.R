clean.summary <- function(arg1, warn=TRUE){
  ##takes output of summary() and returns an html table with median(iqr) for continuous (numeric) variables
  arg1 <- data.frame(arg1)
  quantiles <- c("1st Qu.:", "Median :", "3rd Qu.:")
  vars <- unique(arg1[,"Var2"])
  vals <- matrix(nrow=length(vars), ncol=(length(quantiles)+1))
  for (i in 1:length(vars)){
    vals[i,1] = trimws(vars[i], which="both")
    for (j in 1:length(quantiles)){
      replacer <- arg1[(vars[i]==arg1[,"Var2"]) & grepl(pattern=quantiles[j], x=arg1[,"Freq"]),"Freq"]
      if(!length(replacer)){
        replacer <- "NA"
      }
      if(length(replacer > 1) & warn){
        warning("This function assumes that duplicate column names correspond to duplicate data. Reported values will be wrong if this is not the case. To fix, force unique column names in data frame used that is passed to summary().")
      }
      vals[i,(j+1)] = replacer[1]
    }
  }
  vals <- data.frame(vals)
  colnames(vals)[2:ncol(vals)] <- quantiles
  colnames(vals)[1] <- "vars"
  for (i in 1:length(quantiles)){
    vals[,(i+1)] <- trimws(gsub(quantiles[i], '',vals[,(i+1)]), which="both")
  }
  
  vals$med.iqr <- paste("<strong>", vals[,"Median :"], "</strong>", "<br>", "(", vals[,"1st Qu.:"], "-", vals[,"3rd Qu.:"], ")", sep="")
  out.html.table <- kable(vals[,c("vars", "med.iqr")], "html", escape=FALSE)
  return(out.html.table)
}