optkernel_fun_loocv=function(x,y,mu,theta=seq(from=0,to=10,length=50),max_active=2,heredity=c("strong","weak"),drop.tol=0.01,drop.tol2=0.05,stop.tol=0.1,stop2.tol=0.1,max_iter_1=100,max_iter_2=100,x.norm=TRUE){
  ### mu is the set of tuning parameters for the nugget effect \eta
  ### other parameters are the same in optkernel_method_3.
  n=length(mu)
  m=length(y)
  loocv=rep(0,n) # store the loocv error for each nugget effect
  optK=rep(list(NULL),n)
  y=as.matrix(y)
  for (i in 1:n) {
    optK[[i]]=optkernel3_fun(x,y,mu=mu[i],theta=theta,max_active=max_active,heredity=heredity,drop.tol=drop.tol,drop.tol2=drop.tol2,stop.tol=stop.tol,stop2.tol=stop2.tol,max_iter_1=max_iter_1,max_iter_2=max_iter_2,x.norm=x.norm)
    y_pred=(optK[[i]]$S)%*%y
    for (j in 1:m){
      loocv[i]=loocv[i]+((y[j]-y_pred[j])/(1-optK[[i]]$S[j,j]))^2
    }
    loocv[i]=loocv[i]/m
  }
  optimal=which.min(loocv)
  return(optK[[optimal]])
}