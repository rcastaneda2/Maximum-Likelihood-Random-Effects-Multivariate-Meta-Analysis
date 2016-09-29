
#Multivariate Meta Analytic Function

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












