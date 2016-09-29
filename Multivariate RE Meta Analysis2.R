#Multivariate Meta Analytic Function
library(maxLik)
dat<-read.table("C:/Users/Ruben/Desktop/Dropbox/Research/7. Q3/Kalaian & Raudenbush/DATA CSV.csv",sep=",",header=TRUE)
missing1<-is.na(dat$SATV)
missing2<-is.na(dat$SATM)
dat$SATV[missing1]=0
dat$SATM[missing2]=0
#dat<-na.omit(dat)

x<-as.matrix(cbind(dat$SATV,dat$SATM))
v1<-((dat$ne+dat$nc)/(dat$ne*dat$nc))+((dat$SATV^2)/(2*(dat$ne+dat$nc)))
v2<-((dat$ne+dat$nc)/(dat$ne*dat$nc))+((dat$SATM^2)/(2*(dat$ne+dat$nc)))
v1[missing1]=1000
v2[missing2]=1000

p<-1
i<-1
cor<-.75
cov1<-matrix(nrow=length(x[,1]), ncol=1)
cov2<-matrix(nrow=length(x[,1]), ncol=1)
for(i in 1:length(x[,i])){
  cov1[i]<-((1/dat$ne[i])+(1/dat$nc[i]))*cor + ((.5*( ( x[i,1]%*%t(x[i,2]) ) * cor^2 )) /(dat$ne[i]+dat$nc[i]))
  #cov2[i]<-((1/dat$ne[i])+(1/dat$nc[i]))*cor + (.5*( ( x[i,2]%*%t(x[i,2]) ) * cor^2 ) /(dat$ne[i]+dat$nc[i]))
  cov<-cbind(cov1,cov2)
}
cov


#matrix of vi's 
v<-as.matrix(cbind(v1,v2))

LLMVN<-function (pars) {
  #if (pars[3]>pars[4]) #constrain var to be bigger than cov
  mu  = c(
    pars[1],
    pars[2]
  )
  
  sigma = matrix(c(
    pars[3],
    pars[4],
    pars[4],
    pars[5]
  ),nrow=2,ncol=2,byrow=TRUE)
  
  ll=0
  i=1
  for(i in 1:length(x[,1])){
    mat = matrix(c(
      v[i,1],
      cov[i,1],
      cov[i,1],
      v[i,2]                            
    ),nrow=2,ncol=2,byrow=TRUE)
    
    result=1/2*(log(det(sigma+mat)))+1/2*(t(x[i,]-mu)%*%solve(sigma+mat)%*%(x[i,]-mu))
    ll=result+ll
  } 
  print(c(round(pars,6),round(ll)))
  return(ll)
}

par<-c(.10,.09, .007, -.008, .02)

new<-optim(par,LLMVN,lower=c(-Inf,-Inf,0,-Inf,0),method=('L-BFGS-B'))
new
new<-maxLik(par,LLMVN)
## use constraints: x + y >= 1
A <- matrix(c(1,1,1,1,1), 5, 5)
B <- c(-1,-1,-1,-1,-1)

A <- matrix(c(1,1,1,1,1), , 5)
B <- c(-1)
maxLik(LLMVN,start=par), constraint=list(ineqA=A, ineqB=B))



#Generate Starting Values
par<-c(.10,.09, .007, -.008, .02)
par<-c(0.100, 0.090, 0.070, 0.008, 0.001) #log Lik of 74
par<-c(.103, .099, .00768, -.00835, .02848)
res<-nlminb(par,LLMVN,lower=c(-Inf,-Inf,0,-Inf,0),control=list(iter.max=500))
round(res$par,4)


names(res$par)<-c('SAT-Verbal','SAT-Math','Variance 1','Covariance','Variance 2')
round(res$par,4)      
res

old<-res



new<-optim(par,LLMVN,lower=c(-Inf,-Inf,0,-Inf,0),method=('L-BFGS-B'))
new<-optim(LLMVN,par,lower=c(-Inf,-Inf,0,-Inf,0),method=('L-BFGS-B'))
help(optim)


























