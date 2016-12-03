#=============================================================
#  R CODE FOR
#     "PREDICTION IN HETEROSCEDASTIC NESTED ERROR
#      REGRESSION MODELS WITH RANDOM DISPERSIONS"
#
#     BY TATSUYA KUBOKAWA, SHONOSUKE SUGASAWA,
#          MALAY GHOSH AND SANJAY CHAUDHURI
#
#        PUBLISHED IN STATISTICA SINICA 2016
#=============================================================

#-------------------------------------------------------------
#   FUNCTIONS FOR PARAMETER ESTIMATION AND COMPUTING EBLUP
#-------------------------------------------------------------
# INPUT
# y: responce vector
# X: (N,p) design matrix
# ni: vector of number sample sizes in areas
# C: (m,p) matrix for EBLUP
# maxr: maximum number of iterations for computing MLE

# OUTPUT
# est: vector of parameter estimates
# pred: vector of EBLUP

RHNERM=function(y,X,ni,C,maxr=100){

  m=length(ni); N=sum(ni)
  p=dim(X)[2]
  cum=c(0,cumsum(ni))

  areamean=function(y){
    marea=c()
    for(i in 1:m){ marea[i]=mean(y[(cum[i]+1):cum[i+1]]) }
    return(marea)
  }

  logam=function(z){
    l=length(z); res=c()
    for(i in 1:l){
      x=z[i]
      if(x<170){res[i]=log(gamma(x))}
      else{res[i]=log(2*pi*(x-1))/2+(x-1)*(log(x-1)-1)}
    }
    return(res)
  }

  Q=function(y,mu,la){
    Qi=c()
    for(i in 1:m){
      term=(cum[i]+1):cum[i+1]
      a1=sum((y[term]-mu[term])^2)
      a2=-ni[i]^2*la/(ni[i]*la+1)*(mean(y[term])-mean(mu[term]))^2
      Qi[i]=a1+a2
    }
    return(Qi)
  }

  like.RNER=function(y,para){
    beta=para[1:p]; la=para[p+1]; t1=para[p+2]; t2=para[p+3]
    mu=as.vector(X%*%beta)
    loglike=m*t1*log(t2)+2*sum(logam((ni+t1)/2))-2*m*logam(t1/2)-sum(log(ni*la+1))-sum((ni+t1)*log(Q(y,mu,la)+t2))
    return(-loglike)
  }

  EBLUP=function(y,para){
    beta=para[1:p]; la=para[p+1]
    mu=as.vector(X%*%beta); muC=as.vector(C%*%beta)
    vh=ni*la/(ni*la+1)*(areamean(y)-areamean(mu))
    return(muC+vh)
  }

  dd=1; r=1
  hbeta=solve(t(X)%*%X)%*%t(X)%*%y
  hvar=c(1,3,1)

  while(dd>0.001 & r<=maxr){
    like1=function(para){like.RNER(y,c(hbeta,para))}
    out=optim(par=hvar, fn=like1, gr=NULL, method="L-BFGS-B",
              lower=c(0.01,2.01,0.01), upper=c(1000,500,1000), control=list(),hessian=T)
    new.hvar=out$par
    dd=sum((new.hvar-hvar)^2); hvar=new.hvar
    like2=function(para){like.RNER(y,c(para,hvar))}
    out=optim(par=hbeta, fn=like2, gr=NULL, method="L-BFGS-B",
              lower=rep(-100,p), upper=rep(100,p), control=list(),hessian=T)
    hbeta=out$par
    r=r+1
  }

  name=c(paste0("beta",1:p),"lam","tau1","tau2")
  est=c(hbeta,hvar); names(est)=name

  pred=EBLUP(y,est)
  Res=list(est,pred)
  names(Res)=c("MLE","EB")
  return(Res)
}



#-------------------------------------------------------------
#   FUNCTIONS FOR MSE ESTIMATION
#-------------------------------------------------------------
# INPUT
# y: responce vector
# X: (N,p) design matrix
# ni: vector of number sample sizes in areas
# C: (m,p) matrix for EBLUP
# maxr: maximum number of iterations for computing MLE
# B: number of bootstrap iterations

# OUTPUT
# estMSE: vector of MSE estimates

mseRHNERM=function(y,X,ni,C,maxr=100,B=100){
  m=length(ni); N=sum(ni)
  p=dim(X)[2]
  cum=c(0,cumsum(ni))

  est=RHNERM(y,X,ni,C,maxr=maxr)[[1]]
  est.beta=est[1:p]; est.la=est[p+1]; est.t1=est[p+2]; est.t2=est[p+3]

  areamean=function(y){
    marea=c()
    for(i in 1:m){ marea[i]=mean(y[(cum[i]+1):cum[i+1]]) }
    return(marea)
  }

  g1=function(para){
    la=para[p+1]; t1=para[p+2]; t2=para[p+3]; gam=1/(1+ni*la)
    return((1-gam)*t2/(ni*(t1-2)))
  }

  gen.boot=function(para){
    beta=para[1:p]; la=para[p+1]; t1=para[p+2]; t2=para[p+3]
    mu=X%*%beta; eff=c(); obs=c()
    sig=1/rgamma(m,t1/2,t2/2); eff=rnorm(m,0,sqrt(la*sig))
    for(i in 1:m){
      term=(cum[i]+1):cum[i+1]
      obs[term]=mu[term]+eff[i]+rnorm(ni[i],0,sqrt(sig[i]))
    }
    return(obs)
  }

  boot.g1=matrix(NA,B,m)
  V.beta=array(NA,c(p,p,B)); V.la=c()

  for(b in 1:B){
    boot=gen.boot(est)
    boot.est=RHNERM(boot,X,ni,C,maxr=maxr)[[1]]
    boot.beta=boot.est[1:p]; boot.la=boot.est[p+1]; boot.t1=boot.est[p+2]; boot.t2=boot.est[p+3]
    V.beta[,,b]=(boot.beta-est.beta)%*%t(boot.beta-est.beta)
    V.la[b]=(boot.la-est.la)^2
    boot.g1[b,]=g1(boot.est)
  }

  V.beta=apply(V.beta,c(1,2),mean); V.la=mean(V.la)
  boot.g1=apply(boot.g1,2,mean)
  est.gam=1/(1+ni*est.la)
  est.g1=2*g1(est)-boot.g1
  for(i in 1:m){ est.g1[i]=max(est.g1[i],0) }

  estMSE=est.g1+est.gam^2*diag(C%*%V.beta%*%t(C))+ni*est.gam^3*est.t2/(est.t1-2)*V.la
  return(estMSE)
}



#-------------------------------------------------------------
#   FUNCTIONS FOR CONDITIONAL MSE (CMSE) ESTIMATION
#-------------------------------------------------------------
# INPUT
# y: responce vector
# X: (N,p) design matrix
# ni: vector of number sample sizes in areas
# C: (m,p) matrix for EBLUP
# k: area number
# maxr: maximum number of iterations for computing MLE
# B: number of bootstrap iterations

# OUTPUT
# estCMSE: conditional MSE estimate in the kth area

cmseRHNERM=function(y,X,ni,C,k=1,maxr=100,B=100){
  m=length(ni); N=sum(ni)
  p=dim(X)[2]
  cum=c(0,cumsum(ni))

  est=RHNERM(y,X,ni,C)[[1]]
  est.beta=est[1:p]; est.la=est[p+1]; est.t1=est[p+2]; est.t2=est[p+3]

  areamean=function(y){
    marea=c()
    for(i in 1:m){ marea[i]=mean(y[(cum[i]+1):cum[i+1]]) }
    return(marea)
  }

  EBLUP=function(y,para){
    beta=para[1:p]; la=para[p+1]
    mu=as.vector(X%*%beta); muC=as.vector(C%*%beta)
    vh=ni*la/(ni*la+1)*(areamean(y)-areamean(mu))
    return(muC+vh)
  }

  g1c=function(para,num){
    beta=para[1:p]; la=para[p+1]; t1=para[p+2]; t2=para[p+3]
    term=(cum[num]+1):cum[num+1]
    mu=as.vector(X%*%beta)[term]; cond=y[term]
    Qi=sum((cond-mu)^2)-ni[num]^2*la/(ni[num]*la+1)*(mean(cond)-mean(mu))^2
    return(la*(ni[num]*la+1)^(-1)*(Qi+t2)/(ni[num]+t1-2))
  }

  gen.boot=function(para){
    beta=para[1:p]; la=para[p+1]; t1=para[p+2]; t2=para[p+3]
    mu=X%*%beta; eff=c(); obs=c()
    sig=1/rgamma(m,t1/2,t2/2); eff=rnorm(m,0,sqrt(la*sig))
    for(i in 1:m){
      term=(cum[i]+1):cum[i+1]
      obs[term]=mu[term]+eff[i]+rnorm(ni[i],0,sqrt(sig[i]))
    }
    return(obs)
  }

  boot.g2c=c(); boot.g1c=c()
  term=(cum[k]+1):cum[k+1]; cond=y[term]
  for(b in 1:B){
    boot=gen.boot(est); boot[term]=cond
    boot.res=RHNERM(boot,X,ni,C,maxr=maxr)
    boot.est=boot.res[[1]]; boot.pred=boot.res[[2]]
    pred=EBLUP(boot,est)
    boot.g2c[b]=(pred[k]-boot.pred[k])^2
    boot.g1c[b]=g1c(boot.est,k)
  }

  estCMSE=2*g1c(est,k)-mean(boot.g1c)+mean(boot.g2c)
  names(estCMSE)=paste0("cmse-",k)
  return(estCMSE)
}






