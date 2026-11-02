
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

### Satellite example for N_2, input dimension is 14, n_train=200
B=50
RMSPE_sd <- data.frame(laGP=numeric(B),Opt_kernel=numeric(B),mlegp=numeric(B),MRFA=numeric(B))
theta=seq(from=1,to=9,length=5)
theta=c(theta*10^-2,theta*10^-1,theta,theta*10,theta*10^2)
mu=c(0.005,0.01,0.02,0.05,0.1,0.5)
fp_mlegp=0
fn_mlegp=0
fp_opt=0
fn_opt=0
fp_MRFA=0
fn_MRFA=0
active_true=c(1,2,3,4,5,6,7)
t_lagp=0
t_mlegp=0
t_opt=0
t_MRFA=0

train <- read.table("CD_GRACE_1000_N2.dat")
test <- read.table("CD_GRACE_100_N2.dat")

r <- apply(rbind(train, test)[,1:7], 2, range)

xx_train <- train[,1:7]
xx_test <- test[,1:7]
### standardize the input data
for(j in 1:ncol(xx_train)) {
  xx_train[,j] <- xx_train[,j] - r[1,j]
  xx_test[,j] <- xx_test[,j] - r[1,j]
  xx_train[,j] <- xx_train[,j]/(r[2,j] - r[1,j])
  xx_test[,j] <- xx_test[,j]/(r[2,j] - r[1,j])
}
xx_train=as.matrix(xx_train)
xx_test=as.matrix(xx_test)
yy_train=train[ ,8]
yy_test=test[ ,8]

for (k in 1:B){
  
  
  x_train=xx_train[1:200, ]
  x_test=xx_test
  y_train=yy_train[1:200]
  y_test=yy_test
  
  ### add the fake dimension
  for (i in 1:7){
    m_train=maximinLHS(dim(x_train)[1],7)
    m_test=maximinLHS(dim(x_test)[1],7)
    x_train=cbind(x_train,rowSums(m_train*x_train[ ,1:7]))
    x_test=cbind(x_test,rowSums(m_test*x_test[ ,1:7]))
  }
  
  r <- apply(rbind(x_train,x_test), 2, range)
  for(j in 1:ncol(x_train)) {
    x_train[,j] <- x_train[,j] - r[1,j]
    x_test[,j] <- x_test[,j] - r[1,j]
    x_train[,j] <- x_train[,j]/(r[2,j] - r[1,j])
    x_test[,j] <- x_test[,j]/(r[2,j] - r[1,j])
  }
  
  
  eps=sqrt(.Machine$double.eps)
  t1=Sys.time()
  gpi=newGPsep(x_train,y_train,d=0.1,g=0.1*var(y_train),dK=TRUE)
  mle=mleGPsep(gpi,param="both",tmin=c(eps,eps),tmax=c(10,var(y_train)))
  w=predGPsep(gpi,x_test)
  t2=Sys.time()
  RMSPE_sd$laGP[k] <- sqrt(mean((w$mean-y_test)^2))/sd(y_test)
  t_lagp=t_lagp+(t2-t1)
  
  t1=Sys.time()
  mlegp_model=mlegp(x_train,y_train,nugget=1,parallel = TRUE)
  y_pred=predict(mlegp_model,x_test)
  t2=Sys.time()
  RMSPE_sd$mlegp[k]=sqrt(mean((y_pred-y_test)^2))/sd(y_test)
  t_mlegp=t_mlegp+(t2-t1)
  
  active_mlegp=which(mlegp_model$beta>0.01)
  fp_mlegp=fp_mlegp+length(setdiff(active_mlegp,active_true))
  fn_mlegp=fn_mlegp+length(setdiff(active_true,active_mlegp))
  
  t1=Sys.time()
  optK=optkernel_fun_loocv(x_train,y_train,mu=mu,theta=theta,max_active=4,heredity="strong",drop.tol=0.05,drop.tol2=0.05,stop.tol=0.005,stop2.tol=0.005,max_iter_1=1000,max_iter_2=1000,x.norm=FALSE)
  optK_pred <- optkernel_pred(x_train,optK,newx=x_test)
  t2=Sys.time()
  RMSPE_sd$Opt_kernel[k] <- sqrt(mean((optK_pred-y_test)^2))/sd(y_test)
  t_opt=t_opt+(t2-t1)
  
  active_opt=unique(which(optK$kernel_input,arr.ind = TRUE)[ ,2])
  fp_opt=fp_opt+length(setdiff(active_opt,active_true))
  fn_opt=fn_opt+length(setdiff(active_true,active_opt))
  
  t1=Sys.time()
  MRFA_model=MRFA_fit(x_train,y_train,verbose=FALSE)
  y_pred=predict(MRFA_model,x_test,lambda=min(MRFA_model$lambda))$y_hat
  t2=Sys.time()
  RMSPE_sd$MRFA[k] <- sqrt(mean((y_pred-y_test)^2))/sd(y_test)
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

#ggsave('Satellite_n2_200.pdf')

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

#write.csv(dim_10,file='Satellite_n2_200.csv')
save.image('result/Satellite_n2_200.RData')

rm(list=ls())

source('header.R')

### Satellite example for N_2, input dimension is 14, n_train=500
B=20
RMSPE_sd <- data.frame(laGP=numeric(B),Opt_kernel=numeric(B),mlegp=numeric(B),MRFA=numeric(B))
theta=seq(from=1,to=9,length=5)
theta=c(theta*10^-2,theta*10^-1,theta,theta*10,theta*10^2)
mu=c(0.005,0.01,0.02,0.05,0.1,0.5)
fp_mlegp=0
fn_mlegp=0
fp_opt=0
fn_opt=0
fp_MRFA=0
fn_MRFA=0
active_true=c(1,2,3,4,5,6,7)
t_lagp=0
t_mlegp=0
t_opt=0
t_MRFA=0

train <- read.table("CD_GRACE_1000_N2.dat")
test <- read.table("CD_GRACE_100_N2.dat")

r <- apply(rbind(train, test)[,1:7], 2, range)

xx_train <- train[,1:7]
xx_test <- test[,1:7]
### standardize the input data
for(j in 1:ncol(xx_train)) {
  xx_train[,j] <- xx_train[,j] - r[1,j]
  xx_test[,j] <- xx_test[,j] - r[1,j]
  xx_train[,j] <- xx_train[,j]/(r[2,j] - r[1,j])
  xx_test[,j] <- xx_test[,j]/(r[2,j] - r[1,j])
}
xx_train=as.matrix(xx_train)
xx_test=as.matrix(xx_test)
yy_train=train[ ,8]
yy_test=test[ ,8]

for (k in 1:B){
  
  
  x_train=xx_train[1:500, ]
  x_test=xx_test
  y_train=yy_train[1:500]
  y_test=yy_test
  
  ### add the fake dimension
  for (i in 1:7){
    m_train=maximinLHS(dim(x_train)[1],7)
    m_test=maximinLHS(dim(x_test)[1],7)
    x_train=cbind(x_train,rowSums(m_train*x_train[ ,1:7]))
    x_test=cbind(x_test,rowSums(m_test*x_test[ ,1:7]))
  }
  
  r <- apply(rbind(x_train,x_test), 2, range)
  for(j in 1:ncol(x_train)) {
    x_train[,j] <- x_train[,j] - r[1,j]
    x_test[,j] <- x_test[,j] - r[1,j]
    x_train[,j] <- x_train[,j]/(r[2,j] - r[1,j])
    x_test[,j] <- x_test[,j]/(r[2,j] - r[1,j])
  }
  
  
  eps=sqrt(.Machine$double.eps)
  t1=Sys.time()
  gpi=newGPsep(x_train,y_train,d=0.1,g=0.1*var(y_train),dK=TRUE)
  mle=mleGPsep(gpi,param="both",tmin=c(eps,eps),tmax=c(10,var(y_train)))
  w=predGPsep(gpi,x_test)
  t2=Sys.time()
  RMSPE_sd$laGP[k] <- sqrt(mean((w$mean-y_test)^2))/sd(y_test)
  t_lagp=t_lagp+(t2-t1)
  
  t1=Sys.time()
  mlegp_model=mlegp(x_train,y_train,nugget=1,parallel = TRUE)
  y_pred=predict(mlegp_model,x_test)
  t2=Sys.time()
  RMSPE_sd$mlegp[k]=sqrt(mean((y_pred-y_test)^2))/sd(y_test)
  t_mlegp=t_mlegp+(t2-t1)
  
  active_mlegp=which(mlegp_model$beta>0.01)
  fp_mlegp=fp_mlegp+length(setdiff(active_mlegp,active_true))
  fn_mlegp=fn_mlegp+length(setdiff(active_true,active_mlegp))
  
  t1=Sys.time()
  optK=optkernel_fun_loocv(x_train,y_train,mu=mu,theta=theta,max_active=4,heredity="strong",drop.tol=0.05,drop.tol2=0.05,stop.tol=0.005,stop2.tol=0.005,max_iter_1=1000,max_iter_2=1000,x.norm=FALSE)
  optK_pred <- optkernel_pred(x_train,optK,newx=x_test)
  t2=Sys.time()
  RMSPE_sd$Opt_kernel[k] <- sqrt(mean((optK_pred-y_test)^2))/sd(y_test)
  t_opt=t_opt+(t2-t1)
  
  active_opt=unique(which(optK$kernel_input,arr.ind = TRUE)[ ,2])
  fp_opt=fp_opt+length(setdiff(active_opt,active_true))
  fn_opt=fn_opt+length(setdiff(active_true,active_opt))
  
  t1=Sys.time()
  MRFA_model=MRFA_fit(x_train,y_train,verbose=FALSE)
  y_pred=predict(MRFA_model,x_test,lambda=min(MRFA_model$lambda))$y_hat
  t2=Sys.time()
  RMSPE_sd$MRFA[k] <- sqrt(mean((y_pred-y_test)^2))/sd(y_test)
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

#ggsave('Satellite_n2_500.pdf')

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

#write.csv(dim_10,file='Satellite_n2_500.csv')
save.image('result/Satellite_n2_500.RData')