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

myViterbiMMPP <- function(events, param, termination = NULL, if.print = FALSE){
  ## initialize
  if(is.null(termination)){
    termination <- tail(events, 1)+0.01
  }
  N <- length(events)
  interevent <- diff(c(0,events,termination))
  
  zt_v <- rep(NA,N)
  
  m_matrix <- matrix(0,nrow=N,ncol=2)
  b_matrix <- matrix(0,nrow=N,ncol=2)
  m_termination <- rep(0,2)
  b_termination <- rep(0,2)
  
  probs_1_matrix <- matrix(0,nrow=N+1,ncol=2) ## Probability vector for transition to state 1
  probs_2_matrix <- matrix(0,nrow=N+1,ncol=2) ## Probability vector for transition to state 0
  logp1_1 <- rep(0,N+1)
  logp1_2 <- rep(0,N+1)
  logp2_1 <- rep(0,N+1)
  logp2_2 <- rep(0,N+1)
  
  ## --- Calculation for log probability of Markov transition logP_ij(t)
  for(n in 1:(N+1)){
    probs_1_matrix[n,1] <- log(param$q2/(param$q1+param$q2)+param$q1/(param$q1+param$q2)*exp(-(param$q1+param$q2)*interevent[n])) # 1->1
    probs_2_matrix[n,2] <- log(param$q1/(param$q1+param$q2)+param$q2/(param$q1+param$q2)*exp(-(param$q1+param$q2)*interevent[n])) # 2->2
    probs_1_matrix[n,2] <- log(1-exp(probs_2_matrix[n,2])) # 2->1
    probs_2_matrix[n,1] <- log(1-exp(probs_1_matrix[n,1])) # 1->2
  }
  
  
  ## --- consider n = 1
  #------------ transit to z(t_n) = 1
  logp1_1[1] <- probs_1_matrix[1,1] # 1->1
  logp2_1[1] <- probs_1_matrix[1,2] # 2->1
  if(logp1_1[1]>logp2_1[1]){
    m_matrix[1,1] <- logp1_1[1] + log(param$lambda0*(1+param$c)) - param$lambda0*(1+param$c)*interevent[1]
    b_matrix[1,1] <- 1
  }else{
    m_matrix[1,1] <- logp2_1[1] + log(param$lambda0*(1+param$c)) - param$lambda0*(1+param$c)*interevent[1]
    b_matrix[1,1] <- 2
  }
  
  #------------ transit to z(t_1) = 2
  logp1_2[1] <- probs_2_matrix[1,1]  # 1->2
  logp2_2[1] <- probs_2_matrix[1,2]  # 2->2
  if(logp1_2[1]>logp2_2[1]){
    m_matrix[1,2] <- logp1_2[1] + log(param$lambda0) - param$lambda0*interevent[1]
    b_matrix[1,2] <- 1
  }else{
    m_matrix[1,2] <- logp2_2[1] + log(param$lambda0) - param$lambda0*interevent[1]
    b_matrix[1,2] <- 2
  }
  
  if(N>1){
    ## interate n from 2 to N
    for(n in 2:N){
      #------------ transit to z(t_n) = 1
      logp1_1[n] <- m_matrix[n-1,1] + probs_1_matrix[n,1] # 1->1
      logp2_1[n] <- m_matrix[n-1,2] + probs_1_matrix[n,2] # 2->1
      if(logp1_1[n]>logp2_1[n]){
        m_matrix[n,1] <- logp1_1[n] + log(param$lambda0*(1+param$c)) - param$lambda0*(1+param$c)*interevent[n]
        b_matrix[n,1] <- 1
      }else{
        m_matrix[n,1] <- logp2_1[n] + log(param$lambda0*(1+param$c)) - param$lambda0*(1+param$c)*interevent[n]
        b_matrix[n,1] <- 2
      }
      
      #------------ transit to z(t_n) = 2
      logp1_2[n] <- m_matrix[n-1,1] + probs_2_matrix[n,1] # 1->2
      logp2_2[n] <- m_matrix[n-1,2] + probs_2_matrix[n,2]  # 2->2
      if(logp1_2[n]>logp2_2[n]){
        m_matrix[n,2] <- logp1_2[n] + log(param$lambda0) - param$lambda0*interevent[n]
        b_matrix[n,2] <- 1
      }else{
        m_matrix[n,2] <- logp2_2[n] + log(param$lambda0) - param$lambda0*interevent[n]
        b_matrix[n,2] <- 2
      }
    }
  }
  
  ## termination
  #------------ transit to z(T) = 1
  logp1_1[N+1] <- m_matrix[N,1] + probs_1_matrix[N+1,1] 
  logp2_1[N+1] <- m_matrix[N,2] + probs_1_matrix[N+1,2] 
  if(logp1_1[N+1]>logp2_1[N+1]){
    m_termination[1] <- logp1_1[N+1] - param$lambda0*(1+param$c)*interevent[N+1]
    b_termination[1] <- 1
  }else{
    m_termination[1] <- logp2_1[N+1] - param$lambda0*(1+param$c)*interevent[N+1]
    b_termination[1] <- 2
  }
  
  #------------ transit to z(T) = 2
  logp1_2[N+1] <- m_matrix[N,1] + probs_2_matrix[N+1,1] 
  logp2_2[N+1] <- m_matrix[N,2] + probs_2_matrix[N+1,2] 
  if(logp1_2[N+1]>logp2_2[N+1]){
    m_termination[2] <- logp1_2[N+1]- param$lambda0*interevent[N+1]
    b_termination[2] <- 1
  }else{
    m_termination[2] <- logp2_2[N+1]- param$lambda0*interevent[N+1]
    b_termination[2] <- 2
  }
  
  # ----calulate zt_v for global decoding
  if(m_termination[1]>m_termination[2]){
    termination_state <- 1
  }else{
    termination_state <- 2
  }
  
  zt_v[N] <-  b_termination[termination_state]
  
  if(N > 1){
    for(k in 1:(N-1)){
      zt_v[N-k] = b_matrix[N-k+1,zt_v[N-k+1]];
    }
  }
  
  if(if.print){
    print_matrix <- rbind(events,t(m_matrix),t(b_matrix),zt_v)
    write.csv(print_matrix, "test.csv")
  }
  
  return(list(zt_v=zt_v, initial_state=b_matrix[1,zt_v[1]], termination_state=termination_state))
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
  N <- length(interevent)+1
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
  if(exists("temp.count")){
    return(list(x.hat=c(x.hat[1:temp.count],tail(temp.t,1))+start,z.hat=c(z.hat[1:temp.count],3-z.hat[temp.count])))
  }else{
    return(list(x.hat=x.hat,z.hat=z.hat))
  }
}

mmppIntensityNumeric <-function(params=list(lambda0,c,q1,q2), latent.vec){
  ## latent.vec is vector with same length as time.vec, each entry is the probability at state 1
  lambda.t <- params$lambda0*(1+params$c)*latent.vec + params$lambda0*(1-latent.vec)
  return(lambda.t)
}