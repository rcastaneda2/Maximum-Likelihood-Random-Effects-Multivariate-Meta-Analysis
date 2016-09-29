#Multivariate Meta Analytic Function

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

#matrix of vi's 
v<-as.matrix(cbind(v1,v2))

#Create matrix of correlations
p<-1
i<-1
cor<-.7
for(i in 1:length(x[,i])){
  cov1[i]<-((1/dat$ne[i])+(1/dat$nc[i]))*cor + (( ( x[i,1]%*%t(x[i,1]) ) * cor^2 ) /(dat$ne[i]+dat$nc[i]))
  cov2[i]<-((1/dat$ne[i])+(1/dat$nc[i]))*cor + (( ( x[i,2]%*%t(x[i,2]) ) * cor^2 ) /(dat$ne[i]+dat$nc[i]))
  cov<-cbind(cov1,cov2)
}
cov

LLMVN<-function (pars,x,n) {
  mu  = c(
    pars[1],
    pars[2]
          )

  ll=0
  i=1
  for(i in 1:length(x[,1])){
    mat = matrix(c(
      v[i,1],
      cov[i,1],
      cov[i,1],
      v[i,2]                            
    ),nrow=2,ncol=2,byrow=TRUE)
    
      result=1/2*(log(det(mat)))+1/2*(t(x[i,]-mu)%*%solve(mat)%*%(x[i,]-mu))
      print(log(det(mat)))
      ll=result+ll
   
  } 
  print(ll)
  return(ll)
}





#Generate Starting Values
par<-c(mean(x[1,]),mean(x[2,]))

#0 indicates successful convergence
optim(par,LLMVN)
res<-nlminb(par,LLMVN)
names(res$par)<-c('SAT-Verbal','SAT-Math')
res$par




