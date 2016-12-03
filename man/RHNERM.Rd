\name{RHNERM}
\alias{RHNERM}


\title{
Estimation of random heteroscedastic nested error regression models
}
\description{
Calculates the maximum likelihood estimates of the model parameters in random heteroscedastic nested error regression models. The empirical Bayes estimates of area-level parameters with random effects are also given.
}
\usage{
RHNERM(y, X, ni, C, maxr=100)
}


\arguments{
  \item{y}{
N*1 vector of response values.
}
  \item{X}{
N*p matrix containing N*1 vector of 1 in the first column and vectors of covariates in the rest of columns.
}
  \item{ni}{
m*1 vector of sample sizes in each area.
}
  \item{C}{
m*p matrix of area-level covariates included in the area-level parameters.
}
  \item{maxr}{
maximum number of iteration for computing the maximum likelihood estimates.
}
}


\value{
The function returns a list with the following objects:

\item{MLE}{
(p+3)*1 vector of maximum likelihood estimates of the model parameters.
}
\item{EB}{
m*1 vector of empirical Bayes estimates of the area-level parameters.
}
}

\references{
Kubokawa, K., Sugasawa, S., Ghosh, M. and Chaudhuri, S. (2016). Prediction in Heteroscedastic nested error regression models with random dispersions. Statistica Sinica, 26, 465-492.
}

\author{
Shonosuke Sugasawa
}


\examples{
#generate data
set.seed(1234)
beta=c(1,1); la=1; tau=c(8,4)
m=20; ni=rep(3,m); N=sum(ni)
X=cbind(rep(1,N),rnorm(N))

mu=beta[1]+beta[2]*X[,2]
sig=1/rgamma(m,tau[1]/2,tau[2]/2); v=rnorm(m,0,sqrt(la*sig))
y=c()
cum=c(0,cumsum(ni))
for(i in 1:m){
  term=(cum[i]+1):cum[i+1]
  y[term]=mu[term]+v[i]+rnorm(ni[i],0,sqrt(sig[i]))
}

#fit the random heteroscedastic nested error regression
C=cbind(rep(1,m),rnorm(m))
fit=RHNERM(y,X,ni,C)
fit
}

