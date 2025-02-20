
library(abind)
library(foreach)
library(doParallel)

library(lattice)
library(ggplot2)
library(GPfit)
library(lhs)
library(laGP)
library(mlegp)
library(MRFA)

rm(list=ls())

source('header.R')

### Borehole example for input dimension p=40, n_train=200, n_test=1000
B=50
RMSPE_sd <- data.frame(laGP=numeric(B),Opt_kernel=numeric(B),mlegp=numeric(B),MRFA=numeric(B))
theta=seq(from=1,to=9,length=5)
theta=c(theta*10^-2,theta*10^-1,theta,theta*10,theta*10^2)
mu=c(0.005,0.01,0.02,0.05,0.1,0.5)
p=40
m=1000
n=200
fp_mlegp=0
fn_mlegp=0
fp_opt=0
fn_opt=0
fp_MRFA=0
fn_MRFA=0
active_true=c(1,2,3,4,5,6,7,8)
t_lagp=0
t_mlegp=0
t_opt=0
t_MRFA=0
lower<-c(0.05,100,63070,990,63.1,700,1120,9855,200,500)
lower=rep(lower,4)
upper<-c(0.15,50000,115600,1110,116,820,1680,12045,1000,5000)
upper=rep(upper,4)

for (k in 1:B){
  
  
  ## Training data
  xx<-maximinLHS(n,p)
  minxx<-apply(xx,2,min)
  maxxx<-apply(xx,2,max)
  xx<-scale(xx,minxx,maxxx-minxx)  ## xx is in the range of [0,1]^p, we use the standardized data to train the models
  x<-scale(xx,F,1/(upper-lower))
  x<-scale(x,-lower,F)    ## x is in the range of the phyiscal model
  
  y<-2*pi*x[,3]*(x[,4]-x[,6])/(log(x[,2]/x[,1])*(1+2*x[,7]*x[,3]/(log(x[,2]/x[,1])*x[,1]^2*x[,8])+x[,3]/x[,5]))   ## No noise involved. 
  
  y <- y+ rnorm(n, mean = 0, sd=0.02) ## Add noise
  
  ## Testing data
  x.new<-randomLHS(m,p)
  new<-matrix(0,ncol=p,nrow=m)
  for (i in 1:m)
    new[i,]<-lower+x.new[i,]*(upper-lower)
  
  truey<-2*pi*new[,3]*(new[,4]-new[,6])/(log(new[,2]/new[,1])*(1+2*new[,7]*new[,3]/(log(new[,2]/new[,1])*new[,1]^2*new[,8])+new[,3]/new[,5]))
  
  
  eps=sqrt(.Machine$double.eps)
  t1=Sys.time()
  gpi=newGPsep(xx,y,d=0.1,g=0.1*var(y),dK=TRUE)
  mle=mleGPsep(gpi,param="both",tmin=c(eps,eps),tmax=c(10,var(y)))
  w=predGPsep(gpi,x.new)
  t2=Sys.time()
  RMSPE_sd$laGP[k] <- sqrt(mean((w$mean-truey)^2))/sd(truey)
  t_lagp=t_lagp+(t2-t1)
  
  t1=Sys.time()
  mlegp_model=mlegp(xx,y,nugget=1,parallel = TRUE)
  y_pred=predict(mlegp_model,x.new)
  t2=Sys.time()
  RMSPE_sd$mlegp[k]=sqrt(mean((y_pred-truey)^2))/sd(truey)
  t_mlegp=t_mlegp+(t2-t1)
  
  active_mlegp=which(mlegp_model$beta>0.01)
  fp_mlegp=fp_mlegp+length(setdiff(active_mlegp,active_true))
  fn_mlegp=fn_mlegp+length(setdiff(active_true,active_mlegp))
  
  t1=Sys.time()
  optK=optkernel_fun_loocv(xx,y,mu=mu,theta=theta,max_active=4,heredity="strong",drop.tol=0.05,drop.tol2=0.05,stop.tol=0.005,stop2.tol=0.005,max_iter_1=1000,max_iter_2=1000,x.norm=FALSE)
  optK_pred <- optkernel_pred(xx,optK,newx=x.new)
  t2=Sys.time()
  RMSPE_sd$Opt_kernel[k] <- sqrt(mean((optK_pred-truey)^2))/sd(truey)
  t_opt=t_opt+(t2-t1)
  
  active_opt=unique(which(optK$kernel_input,arr.ind = TRUE)[ ,2])
  fp_opt=fp_opt+length(setdiff(active_opt,active_true))
  fn_opt=fn_opt+length(setdiff(active_true,active_opt))
  
  t1=Sys.time()
  MRFA_model=MRFA_fit(xx,y,verbose=TRUE)
  y_pred=predict(MRFA_model,x.new,lambda=min(MRFA_model$lambda))$y_hat
  t2=Sys.time()
  RMSPE_sd$MRFA[k] <- sqrt(mean((y_pred-truey)^2))/sd(truey)
  t_MRFA=t_MRFA+(t2-t1)
  
  active_MRFA=NULL
  for (i in 1:length(MRFA_model$active.group)){
    for (j in 1:length(MRFA_model$active.group[[i]])){
      active_MRFA=c(active_MRFA,MRFA_model$active.group[[i]][[j]]$effect)
    }
  }
  active_MRFA=unique(active_MRFA)
  fp_MRFA=fp_MRFA+length(setdiff(active_MRFA,active_true))
  fn_MRFA=fn_MRFA+length(setdiff(active_true,active_MRFA))
}

RMSPE_sd_box<-data.frame(RMSPE=c(RMSPE_sd$laGP, RMSPE_sd$Opt_kernel,RMSPE_sd$mlegp,RMSPE_sd$MRFA),type=rep(colnames(RMSPE_sd),each=B))
ggplot(RMSPE_sd_box, aes(x=as.factor(type), y=RMSPE)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  xlab("method")+ylab("Standardized RMSPE")

#ggsave('Borehole_40dim_200_50.pdf')

RMSPE_sd_lagp=mean(RMSPE_sd$laGP)
RMSPE_sd_mlegp=mean(RMSPE_sd$mlegp)
RMSPE_sd_opt=mean(RMSPE_sd$Opt_kernel)
RMSPE_sd_MRFA=mean(RMSPE_sd$MRFA)
fp_mlegp=fp_mlegp/B
fn_mlegp=fn_mlegp/B
fp_opt=fp_opt/B
fn_opt=fn_opt/B
fp_MRFA=fp_MRFA/B
fn_MRFA=fn_MRFA/B
t_lagp=t_lagp/B
t_mlegp=t_mlegp/B
t_opt=t_opt/B
t_MRFA=t_MRFA/B
dim_10=data.frame(method=c('lagp','mlegp','MRFA','optK'),rmspe_sd=c(RMSPE_sd_lagp,RMSPE_sd_mlegp,RMSPE_sd_MRFA,RMSPE_sd_opt),time=c(t_lagp,t_mlegp,t_MRFA,t_opt),fp=c('/',fp_mlegp,fp_MRFA,fp_opt),fn=c('/',fn_mlegp,fn_MRFA,fn_opt))

#write.csv(dim_10,file='Borehole_40dim_200_50.csv')
save.image('result/Borehole_40dim_200_50.RData')

rm(list=ls())

source('header.R')

### Borehole example for input dimension p=40, n_train=500, n_test=1000
B=20
RMSPE_sd <- data.frame(laGP=numeric(B),Opt_kernel=numeric(B),mlegp=numeric(B),MRFA=numeric(B))
theta=seq(from=1,to=9,length=5)
theta=c(theta*10^-2,theta*10^-1,theta,theta*10,theta*10^2)
mu=c(0.005,0.01,0.02,0.05,0.1,0.5)
p=40
m=1000
n=500
fp_mlegp=0
fn_mlegp=0
fp_opt=0
fn_opt=0
fp_MRFA=0
fn_MRFA=0
active_true=c(1,2,3,4,5,6,7,8)
t_lagp=0
t_mlegp=0
t_opt=0
t_MRFA=0
lower<-c(0.05,100,63070,990,63.1,700,1120,9855,200,500)
lower=rep(lower,4)
upper<-c(0.15,50000,115600,1110,116,820,1680,12045,1000,5000)
upper=rep(upper,4)

for (k in 1:B){
  
  
  ## Training data
  xx<-maximinLHS(n,p)
  minxx<-apply(xx,2,min)
  maxxx<-apply(xx,2,max)
  xx<-scale(xx,minxx,maxxx-minxx)  ## xx is in the range of [0,1]^p, we use the standardized data to train the models
  x<-scale(xx,F,1/(upper-lower))
  x<-scale(x,-lower,F)    ## x is in the range of the phyiscal model
  
  y<-2*pi*x[,3]*(x[,4]-x[,6])/(log(x[,2]/x[,1])*(1+2*x[,7]*x[,3]/(log(x[,2]/x[,1])*x[,1]^2*x[,8])+x[,3]/x[,5]))   ## No noise involved. 
  
  y <- y+ rnorm(n, mean = 0, sd=0.02) ## Add noise
  
  ## Testing data
  x.new<-randomLHS(m,p)
  new<-matrix(0,ncol=p,nrow=m)
  for (i in 1:m)
    new[i,]<-lower+x.new[i,]*(upper-lower)
  
  truey<-2*pi*new[,3]*(new[,4]-new[,6])/(log(new[,2]/new[,1])*(1+2*new[,7]*new[,3]/(log(new[,2]/new[,1])*new[,1]^2*new[,8])+new[,3]/new[,5]))
  
  
  eps=sqrt(.Machine$double.eps)
  t1=Sys.time()
  gpi=newGPsep(xx,y,d=0.1,g=0.1*var(y),dK=TRUE)
  mle=mleGPsep(gpi,param="both",tmin=c(eps,eps),tmax=c(10,var(y)))
  w=predGPsep(gpi,x.new)
  t2=Sys.time()
  RMSPE_sd$laGP[k] <- sqrt(mean((w$mean-truey)^2))/sd(truey)
  t_lagp=t_lagp+(t2-t1)
  
  t1=Sys.time()
  mlegp_model=mlegp(xx,y,nugget=1,parallel = TRUE)
  y_pred=predict(mlegp_model,x.new)
  t2=Sys.time()
  RMSPE_sd$mlegp[k]=sqrt(mean((y_pred-truey)^2))/sd(truey)
  t_mlegp=t_mlegp+(t2-t1)
  
  active_mlegp=which(mlegp_model$beta>0.01)
  fp_mlegp=fp_mlegp+length(setdiff(active_mlegp,active_true))
  fn_mlegp=fn_mlegp+length(setdiff(active_true,active_mlegp))
  
  t1=Sys.time()
  optK=optkernel_fun_loocv(xx,y,mu=mu,theta=theta,max_active=4,heredity="strong",drop.tol=0.05,drop.tol2=0.05,stop.tol=0.005,stop2.tol=0.005,max_iter_1=1000,max_iter_2=1000,x.norm=FALSE)
  optK_pred <- optkernel_pred(xx,optK,newx=x.new)
  t2=Sys.time()
  RMSPE_sd$Opt_kernel[k] <- sqrt(mean((optK_pred-truey)^2))/sd(truey)
  t_opt=t_opt+(t2-t1)
  
  active_opt=unique(which(optK$kernel_input,arr.ind = TRUE)[ ,2])
  fp_opt=fp_opt+length(setdiff(active_opt,active_true))
  fn_opt=fn_opt+length(setdiff(active_true,active_opt))
  
  t1=Sys.time()
  MRFA_model=MRFA_fit(xx,y,verbose=TRUE)
  y_pred=predict(MRFA_model,x.new,lambda=min(MRFA_model$lambda))$y_hat
  t2=Sys.time()
  RMSPE_sd$MRFA[k] <- sqrt(mean((y_pred-truey)^2))/sd(truey)
  t_MRFA=t_MRFA+(t2-t1)
  
  active_MRFA=NULL
  for (i in 1:length(MRFA_model$active.group)){
    for (j in 1:length(MRFA_model$active.group[[i]])){
      active_MRFA=c(active_MRFA,MRFA_model$active.group[[i]][[j]]$effect)
    }
  }
  active_MRFA=unique(active_MRFA)
  fp_MRFA=fp_MRFA+length(setdiff(active_MRFA,active_true))
  fn_MRFA=fn_MRFA+length(setdiff(active_true,active_MRFA))
}

RMSPE_sd_box<-data.frame(RMSPE=c(RMSPE_sd$laGP, RMSPE_sd$Opt_kernel,RMSPE_sd$mlegp,RMSPE_sd$MRFA),type=rep(colnames(RMSPE_sd),each=B))
ggplot(RMSPE_sd_box, aes(x=as.factor(type), y=RMSPE)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  xlab("method")+ylab("Standardized RMSPE")

#ggsave('Borehole_40dim_500_20.pdf')

RMSPE_sd_lagp=mean(RMSPE_sd$laGP)
RMSPE_sd_mlegp=mean(RMSPE_sd$mlegp)
RMSPE_sd_opt=mean(RMSPE_sd$Opt_kernel)
RMSPE_sd_MRFA=mean(RMSPE_sd$MRFA)
fp_mlegp=fp_mlegp/B
fn_mlegp=fn_mlegp/B
fp_opt=fp_opt/B
fn_opt=fn_opt/B
fp_MRFA=fp_MRFA/B
fn_MRFA=fn_MRFA/B
t_lagp=t_lagp/B
t_mlegp=t_mlegp/B
t_opt=t_opt/B
t_MRFA=t_MRFA/B
dim_10=data.frame(method=c('lagp','mlegp','MRFA','optK'),rmspe_sd=c(RMSPE_sd_lagp,RMSPE_sd_mlegp,RMSPE_sd_MRFA,RMSPE_sd_opt),time=c(t_lagp,t_mlegp,t_MRFA,t_opt),fp=c('/',fp_mlegp,fp_MRFA,fp_opt),fn=c('/',fn_mlegp,fn_MRFA,fn_opt))

#write.csv(dim_10,file='Borehole_40dim_500_20.csv')
save.image('result/Borehole_40dim_500_20.RData')





