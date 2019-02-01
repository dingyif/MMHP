iAndSI<-function(M){
  n=ncol(M)
  result=rep(0,2)
  k=(M-t(M))/2
  k[upper.tri(k)]=0
  result[1] = 0
  a=length(which(k>0))
  y=which(k>0)%%n
  x=(which(k>0)-1)%/%n+1
  y[y==0]=y[y==0]+n
  if (a>0){
    result[1]=a
    result[2]=sum(y-x)
  }
  return(result)
}