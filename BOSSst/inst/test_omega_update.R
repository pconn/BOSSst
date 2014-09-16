#test omega update
library(mvtnorm)

#simulate some data
N=100000
S=1000
n.iter=100000
n.sampled=100
Which.sampled=sample(c(1:S),n.sampled)
Cov=rnorm(S,0,1)
tau.eps=20
Epsilon=rnorm(S,0,1/sqrt(tau.eps))
Omega=Cov-Cov^2+Epsilon #implied beta = 1, -1
Omega.obs=Omega[Which.sampled]
Pi=exp(Omega)/sum(exp(Omega))
Pi.obs.true=exp(Omega.obs)/sum(exp(Omega.obs))
N.s=rmultinom(1,N,Pi)
r=0.2
C=rbinom(n.sampled,N.s[Which.sampled],r)


#estimation - here, for this one, no unobserved component
X=matrix(0,n.sampled,2)
X[,1]=Cov[Which.sampled]
X[,2]=X[,1]^2
XpXinv=solve(crossprod(X))
XpXinvXp=solve(crossprod(X))%*%t(X)
MCMC.beta=matrix(0,2,n.iter)
MCMC.omega=matrix(0,n.sampled,n.iter)
MCMC.omega[,1]=Omega.obs+2+rnorm(n.sampled,0,0.4)
MCMC.beta[,1]=t(rmvnorm(1,XpXinvXp%*%MCMC.omega[,1],XpXinv/tau.eps+0.001))
MCMC.tau=rep(0,n.iter)
MCMC.tau[1]=1
Pi.obs=exp(MCMC.omega[,1])/sum(exp(MCMC.omega[,1]))
Accept=rep(0,n.sampled)
MH.omega=0.5
MH.Sigma=diag(n.sampled)
Omega.current=matrix(0,n.sampled,100)
Omega.current[,1]=Omega.obs
Old.accept=0
for(iiter in 2:n.iter){
  if((iiter%%1000)==1)cat(paste("iteration ",iiter," of ",n.iter,"\n"))
  #update fixed effect parameters
  MCMC.beta[,iiter]=t(rmvnorm(1,XpXinvXp%*%MCMC.omega[,iiter-1],XpXinv/tau.eps+0.001))
  
  #update omega
  
  #   #unobserved not needed for this one
  #     #observed
  #     MCMC.omega[,iiter]=MCMC.omega[,iiter-1]
  #     Pred=X%*%MCMC.beta[,iiter]
  #     Dt.old=0.5*MH.omega^2*d_logP_omega(Omega=MCMC.omega[,iiter],Count=C,Mu=Pred,tau=tau.eps,G.sum=sum(C),Cur.thin=rep(1,100))
  #     Prop=MCMC.omega[,iiter]+Dt.old+rnorm(100,0,MH.omega)
  #     Prop.exp=exp(Prop)
  #     Pi.prop=Prop.exp/sum(Prop.exp)    
  #     post.new=sum(dnorm(Prop,Pred,1/sqrt(tau.eps),log=TRUE))+sum(C*log(Pi.prop))
  #     post.old=sum(dnorm(MCMC.omega[,iiter],Pred,1/sqrt(tau.eps),log=TRUE))+sum(C*log(Pi.obs))
  #     Temp=Prop-MCMC.omega[,iiter]-Dt.old
  #     jump.old.to.new=-0.5/MH.omega^2*crossprod(Temp)
  #     Dstar=0.5*MH.omega^2*d_logP_omega(Omega=Prop,Count=C,Mu=Pred,tau=tau.eps,G.sum=sum(C),Cur.thin=rep(1,100))
  #     Temp=MCMC.omega[,iiter]-Prop-Dstar
  #     jump.new.to.old=-0.5/MH.omega^2*crossprod(Temp)
  #     if(runif(1)<exp(post.new-post.old+jump.new.to.old-jump.old.to.new)){
  #         Pi.obs=Pi.prop
  #         MCMC.omega[,iiter]=Prop
  #         Accept[1]=Accept[1]+1
  #     }  
  
  #multi - par MH
  #   sd=1/sqrt(MCMC.tau[iiter-1])
  #   MCMC.omega[,iiter]=MCMC.omega[,iiter-1]
  #   Pred=X%*%MCMC.beta[,iiter]
  #   Prop=MCMC.omega[,iiter-1]+0.03*rnorm(n.sampled)
  #   Pi.obs.prop=exp(Prop)/sum(exp(Prop))
  #   Cur.sum=sum(Pi.obs.prop*r)
  #   newL=sum(dnorm(Prop,Pred,sd,log=TRUE))+sum(C*log(Pi.obs.prop))-sum(C)*log(Cur.sum)
  #   Cur.sum=sum(Pi.obs*r)
  #   oldL=sum(dnorm(MCMC.omega[,iiter],Pred,sd,log=TRUE))+sum(C*log(Pi.obs))-sum(C)*log(Cur.sum)
  #   if(runif(1)<exp(newL-oldL)){
  #     Pi.obs=Pi.obs.prop
  #     MCMC.omega[,iiter]=Prop
  #     Accept[1]=Accept[1]+1
  #   }
  
  #adapting multi-par MH
  sd=1/sqrt(MCMC.tau[iiter-1])
  MCMC.omega[,iiter]=MCMC.omega[,iiter-1]
  Pred=X%*%MCMC.beta[,iiter]
  Prop=MCMC.omega[,iiter-1]+MH.omega*rmvnorm(1,rep(0,n.sampled),MH.Sigma,method='chol')
  Pi.obs.prop=exp(Prop)/sum(exp(Prop))
  Cur.sum=sum(Pi.obs.prop*r)
  newL=sum(dnorm(Prop,Pred,sd,log=TRUE))+sum(C*log(Pi.obs.prop))-sum(C)*log(Cur.sum)
  Cur.sum=sum(Pi.obs*r)
  oldL=sum(dnorm(MCMC.omega[,iiter],Pred,sd,log=TRUE))+sum(C*log(Pi.obs))-sum(C)*log(Cur.sum)
  if(runif(1)<exp(newL-oldL)){
    Pi.obs=Pi.obs.prop
    MCMC.omega[,iiter]=Prop
    Accept[1]=Accept[1]+1
  }  
  cur.j=iiter%%100
  if(cur.j==0)cur.j=100
  Omega.current[,cur.j]=MCMC.omega[,iiter]
  
  #update tau
  Diff=MCMC.omega[,iiter]-Pred
  tau.eps <- rgamma(1,n.sampled/2 + 1, as.numeric(crossprod(Diff,Diff))*0.5 + 0.01)
  MCMC.tau[iiter]=tau.eps
  
  if(iiter%%100==0){
    #adapt Omega proposal covariance matrix using Alg 2 of Shaby and Wells (2010)
    R=(Accept[1]-Old.accept)/100
    Sample.Cov=cov(t(Omega.current))
    gamma.1=(iiter/100)^-0.8
    MH.omega=sqrt(MH.omega^2*exp(gamma.1*(R-0.234)))
    MH.Sigma=MH.Sigma+gamma.1*(Sample.Cov-MH.Sigma)
    if(iiter==100)MH.Sigma=MH.Sigma+0.1*diag(n.sampled)
    Old.accept=Accept[1] 
  }  
}

# Omega.est=
# #in this version, model binomial contribution
# for(iiter in 2:n.iter){
#   #update fixed effect parameters
#   MCMC.beta[,iiter]=t(rmvnorm(1,XpXinvXp%*%MCMC.omega[,iiter-1],XpXinv/tau.eps+0.001))
#   
#   #update tau
#   tau.eps=20
#   
#   #update omega
#   
#   #unobserved not needed for this one
#   
#   #observed
#   MCMC.omega[,iiter]=MCMC.omega[,iiter-1]
#   Pred=X%*%MCMC.beta[,iiter]
#   Prop=MCMC.omega[,iiter-1]
#   for(ipar in 1:100){
#     prop=MCMC.omega[ipar,iiter-1]+rnorm(1,0,0.8)
#     Prop[ipar]=prop
#     Pi.obs.prop=exp(Prop)/sum(exp(Prop))
#     oldL=dnorm(MCMC.omega[ipar,iiter-1],Pred[ipar],1/sqrt(tau.eps),log=TRUE)+sum(C*log(Pi.obs))
#     newL=dnorm(prop,Pred[ipar],1/sqrt(tau.eps),log=TRUE)+sum(C*log(Pi.obs.prop))
#     if(runif(1)<exp(newL-oldL)){
#       Pi.obs=Pi.obs.prop
#       MCMC.omega[ipar,iiter]=prop
#       Accept[ipar]=Accept[ipar]+1
#     }
#     else Prop[ipar]=MCMC.omega[ipar,iiter-1]
#   }
# }
