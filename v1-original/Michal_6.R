michal_6 <- function(xx)
{
  ##########################################################################
  #
  # MICHALEWICZ FUNCTION
  #
  # The version used in Haaland and Qian (2011)
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx is an N-by-2 matrix of input, in the range of [0,1]^2
  #
  ##########################################################################
  
  if (is.vector(xx)) {
    xx <- matrix(xx,nrow=1,ncol=6)
  }
  y <- sin(pi*xx[,1])*(sin(pi*xx[,1]^2))^20+sin(pi*xx[,2])*(sin(2*pi*xx[,2]^2))^20+sin(pi*xx[,3])*(sin(pi*xx[,3]^2))^20+sin(pi*xx[,4])*(sin(pi*xx[,4]^2))^20+sin(pi*xx[,5])*(sin(2*pi*xx[,5]^2))^20+sin(pi*xx[,6])*(sin(pi*xx[,6]^2))^20
  return(y)
}