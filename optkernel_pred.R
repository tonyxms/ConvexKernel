optkernel_pred<-function(x,optK,newx){
  ### x is the training data and newx is the test data
  ### optK is the optimal kernel we get
  
  
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.matrix(newx)) newx <- as.matrix(newx)
  n <- nrow(x)
  p <- ncol(x)
  m <- nrow(newx)
  
  if (optK$normx$x.norm) {
    x <- scale(x,optK$normx$x.mean,optK$normx$x.scale)
    newx <- scale(newx,optK$normx$x.mean,optK$normx$x.scale)
  }
  
  L <- length(optK$weight)
  ### Compute the pairwised distance between newx and x
  mzmax <- m*n
  D <- matrix(0,nrow=mzmax,ncol=p)
  for (k in 1:m) {
    lines_ix <- ((k-1)*n+1):(k*n)
    D[lines_ix,] <- (newx[rep(k,n),]-x)^2
  }
  
  ### Get the prediction for the test data
  r <- matrix(0,nrow=m,ncol=n)
  for (i in 1:m) {
    line_ix <- ((i-1)*n+1):(i*n)
    for (j in 1:L) {
      valid_id <- optK$kernel_input[j,]
      theta <- optK$theta[j]
      if (sum(valid_id)==1) 
        add_r <- optK$weight[j]*exp(-theta*D[line_ix,valid_id])
      else 
        add_r <- optK$weight[j]*exp(-theta*rowSums(D[line_ix,valid_id]))
      
      r[i,] <- r[i,] + drop(add_r)
    }
  }
  predy <- r %*% drop(optK$coeff)
  return(drop(predy))
}