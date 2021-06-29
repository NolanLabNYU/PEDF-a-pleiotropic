uniqueEl <- function(arg1) {
  ##computes the number of unique elements in a pair of sets, for all pairs of columns in a matrix
  ##George Crowley
  if(!is.matrix(arg1)){
    stop("input argument is not a matrix")
  }
  n = ncol(arg1)
  unique.mat <- matrix(,ncol= n, nrow = n)
  for (i in 1:n){
    for (j in 1:n) {
      unique.mat[i,j] <- length(setdiff(union(arg1[,i], arg1[,j]), intersect(arg1[,i], arg1[,j])))
    }
  }
  return(unique.mat)
}