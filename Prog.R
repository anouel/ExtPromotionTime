# rm(list=ls())


#----------------------------------------------------------------
# needed R packages.
#----------------------------------------------------------------
library('survival')
library("optimx")
# optimx options
local_opts <- list( "algorithm" = "NLOPT_LD_MMA","xtol_rel"  = 1.0e-14 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG","xtol_rel"  = 1.0e-14,"maxeval"   = 1000,"local_opts" = local_opts ) 
#----------------------------------------------------------------


#----------------------------------------------------------------
# Some utility functions
#----------------------------------------------------------------
Q=function(Vec,Delta,W)
{
  n=length(Delta)
  vw=(cumsum((W*Vec*Delta)[n:1]))[n:1]+sum(W*Vec*(1-Delta))
  return(vw/n)
}
# Log-likelihood
PrfLogLik=function(betA,delta,Delta,X,W)
{
  gX=apply(X,1,g,betA=betA)
  QX=Q(Vec=gX,Delta=Delta,W=W)
  gX[(delta==0)|(W==0)]=1;QX[(delta==0)|(W==0)]=1
  LL=sum(W*delta*log(gX/QX),na.rm=TRUE)
  return(LL)
}
# Gradient
dPrfLogLik=function(betA,delta,Delta,X,W)
{
  gX=apply(X,1,g,betA=betA)
  dgX=apply(X,1,dg,betA=betA)
  dgX=t(dgX)
  QX=Q(gX,Delta,W)
  dQX=apply(dgX,2,Q,Delta,W)
  gX[(delta==0)|(W==0)]=1;QX[(delta==0)|(W==0)]=1
  dLL=apply(W*delta*(dgX/gX-dQX/QX),2,sum,na.rm=TRUE)
  return(dLL)
}
# Hessian
ddPrfLogLik=function(betA,delta,Delta,X,W)
{
  d=ncol(X)
  gX=apply(X,1,g,betA=betA)
  dgX=t(apply(X,1,dg,betA=betA))
  QX=Q(gX,Delta,W)
  dQX=apply(dgX,2,Q,Delta,W)
  ddetaBX=apply(X,1,g,betA=betA,ETA=ddeta)
  gX[(delta==0)|(W==0)]=1;QX[(delta==0)|(W==0)]=1
  ddLL=matrix(NA,ncol=d,nrow=d)
  for(k in 1:d){
    for(l in 1:k)
    {
      ddgX.kl=ddetaBX*X[,k]*X[,l]
      ddQX.kl=Q(ddgX.kl,Delta,W) 
      ddLL.kl=W*delta*((ddgX.kl*gX-dgX[,k]*dgX[,l])/(gX^2)-(ddQX.kl*QX-dQX[,k]*dQX[,l])/(QX^2))
      ddLL[k,l]=sum(ddLL.kl,na.rm=TRUE)
    }}
  ddLL[upper.tri(ddLL)]=ddLL[lower.tri(ddLL)]
  return(ddLL)
}
PPrfLogLik=function(...){-PrfLogLik(...)}
dPPrfLogLik=function(...){-dPrfLogLik(...)}
ddPPrfLogLik=function(...){-ddPrfLogLik(...)}
#----------------------------------------------------------------


#----------------------------------------------------------------
# The main function
#----------------------------------------------------------------
estF<-function(data)
{
  OT<-data[,1];delta<-data[,2];X<-data[,-c(1,2)]
  tau <- max(OT[delta==1])
  # "observed" uncure indicator
  ODelta <- (OT<=tau) 
  OTO=order(OT)
  delta=delta[OTO]
  X=as.matrix(X[OTO,])
  OT=OT[OTO]
  ODelta=ODelta[OTO]
  W=rep(1,length(OTO))
  d=ncol(X);n=nrow(X)
  # starting value
  betAI=coxph(Surv(OT,delta)~.,data=data)$coef
  # The estimation of gamma's
  An=optimx(fn=PPrfLogLik,gr=dPPrfLogLik,hess=ddPPrfLogLik,par=betAI,delta=delta,Delta=ODelta,X=X,W=W,method=c('nlminb','newuoa','bobyqa')) 
  vAn=-An$value; an=coef(An)
  bmt=which.max(vAn)
  betAE=as.vector(an[bmt,]) # gamma's
  # The estimation of theta
  gX=apply(X,1,g,betA=betAE)
  QX=Q(gX,ODelta,W)
  thetaE=mean((W*delta)/QX,na.rm=TRUE) # theta
  # VarCov gamma's
  dgX=t(apply(X,1,dg,betA=betAE))
  dQX=apply(dgX,2,Q,ODelta,W)
  D=NULL;H=NULL;HQ=NULL
  I1=matrix(0,ncol=d,nrow=d)
  for(i in 1:n){
    di=(dgX[i,]/gX[i])*delta[i];D=rbind(D,di=di)
    hi=(dQX[i,]/QX[i])*delta[i];H=rbind(H,hi=hi)
    hqi=hi/QX[i];HQ=rbind(HQ,hqi=hqi)
    I1=I1+((di-hi)%*%t(di-hi))
  }
  I1=try(solve(I1/n))   #VarCov(gamma's)
  hq=apply(HQ,2,mean)
  vth1=mean(delta/QX^2)+t(hq)%*%I1%*%hq #Var(theta)
  return(list(Gamma=betAE,theta=thetaE,VarCov_Gamma=I1,Var_theta=vth1))
}
#----------------------------------------------------------------
