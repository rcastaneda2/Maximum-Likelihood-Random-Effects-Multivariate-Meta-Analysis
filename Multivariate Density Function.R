library(mvtnorm)

#Generate data with known mean and covariance
#First generate the covariance matrix that I want (I want .5, 0, 0, .5)

sigma<-matrix(c(1.5,0.6,0.6,1.5),nrow=2,ncol=2,byrow=TRUE)




#next specify the multivariate normal I want with mean of 1 for both vars and var of .5 for both vars. 

x <- rmvnorm(n=1000, mean=c(1,1),sigma=sigma)  #1000 observations



#covariance structure is sigma which is also

round(cov(rmvnorm(n=100000,mean=c(1,1),sigma=sigma)),2) 


#now that I have generated the data that I wanted to. It is time to retrieve the data using my function. 

MVN<-function (pars) {
  mu  = c(
    pars[1],
    pars[2]
          )
  mat = matrix(c(
    pars[3],
    pars[4],
    pars[4],
    pars[3]
              ),nrow=2,ncol=2,byrow=TRUE)
  ll=0
  i=1
  for(i in 1:length(x[,1])){
      result=1/2*(log(det(mat)))+1/2*(t(x[i,]-mu)%*%solve(mat)%*%(x[i,]-mu))
      ll=result+ll
  } 
  return(ll)
}




#Generate Starting Values
par<-c(0,0,.05,.04)


nlminb(par,MVN)












