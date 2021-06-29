ham.dist <- function(arg1) {
  ##function ham.dist computes all pairwise hamming distances among the columns of arg1, a matrix
  ##George Crowley
  if(!is.matrix(arg1)) {
    stop("input argument is not a matrix")
  }
  n = ncol(arg1)
  dist.mat <- matrix(, ncol = n, nrow = n)
  for (i in 1:n) {
    for (j in 1:n) {
      dist.mat[i,j] <-  sum(arg1[,i] != arg1[,j])
    }
  }
  return(dist.mat)
}