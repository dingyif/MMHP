#---------This function is used to compute compensator for mmpp
mmppCompensator<-function(params=list(lambda0,c,q1,q2), t, pzt){
  lambda0 <- params$lambda0
  c <- params$c
  q1 <- params$q1
  q2 <- params$q2
  n <- length(t) - 1
  interevent <- t[-1]-t[-(n+1)]
  
  Lambda_mixed <- lambda0*(1+c)*interevent*pzt + lambda0*interevent*(1-pzt)
  return(Lambda_mixed)
}
mmppModifiedLatentTrajectory <- function(params=list(lambda0,c,q1,q2), interevent, zt, start=0){
  # input pars
  #       interevent: length N
  #       zt: inferred latent state, 1 or 2 , length N
  # output  z.hat: states of Markov process
  #         x.hat: time of each transition of Markov process
  lambda0 <- params$lambda0
  c <- params$c
  q1 <- params$q1
  q2 <- params$q2
  N <- length(interevent)
  temp.t <- cumsum(interevent)
  if(length(unique(zt))==1){
    z.hat <- unique(zt)
    x.hat <- 0
  }else{
    z.hat <- rep(NA,sum(diff(zt)!=0)+1)
    x.hat <- rep(NA,sum(diff(zt)!=0)+1)
    z.hat[1] <- zt[1]
    x.hat[1] <- 0
    temp.count <- 1
    for(l in 2:N){
      if(zt[l]==1 & zt[l-1]==2){
        temp.count <- temp.count + 1
        x.hat[temp.count] <- temp.t[l]
        z.hat[temp.count] <- 1
      }
      if(zt[l]==2 & zt[l-1]==1){
        temp.count <- temp.count + 1
        x.hat[temp.count] <- temp.t[l-1]
        z.hat[temp.count] <- 2
      }
    }
  }
  return(list(x.hat=c(x.hat,tail(temp.t,1))+start,z.hat=c(z.hat,tail(z.hat,1))))
}