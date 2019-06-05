uniHawkesCompensator<-function(lambda0,alpha,beta,t){
  n<-length(t)
  delta.t<-t-c(0,t[-n])
  Lambda<-rep(0,n)
  A<-0
  Lambda[1]<-lambda0*(t[1])*2
  for(i in 2:n){
    A<-1+exp(-beta*(delta.t[i-1]))*A
    Lambda[i]<-lambda0*(delta.t[i])*2+alpha/beta*(1-exp(-beta*delta.t[i]))*A
  }
  return(Lambda)
}