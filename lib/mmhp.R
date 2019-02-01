
################################################################################################
#This R file is source code to including functions related to Markov Modulated Hawkes Process (MMHP)
#1. simulated.mmhp: to simulate MMHP using thinning theorem
#2. simulate.hp: to simulate Hawkes process during state 1 (used in function 'simulated.mmhp')
#3. simulate.mmhp.given.state: to simulate MMHP given hidden state trajectory
#4. drawUniMMHPIntensity: plot intensity function for MMHP
#5. drawHPIntensity: a helper function for 'drawUniMMHPIntensity', to draw Hawkes process intensity
#6. mmhpNegLogLik: compute negative log likelihood for mmhp
#7. mmhpCompensator: compute vector of compensators for mmhp
#8. mmhpViterbi: infer the probability of zt=1 for each event using Viterbi algorithm
#9. mmhpLocalState: infer the probability of zt=1 for each event using local decoding
#10. modifiedLatentTrajectory: interpolate the latent state trajectory given estimated states of events
################################################################################################

##Simulate Markov Modulated Hawkes Process (including all the history)
simulate.mmhp <- function(object, nsim=1, seed=NULL, ...){
  # VALID FOR 2 STATES: 1: active, 2: not active
  # input q1,q2: transition probability
  #       delta: initial state probability; 
  #       lambda1, alpha, beta : parameters for Hawkes process
  #       lambda0: parameters for homogeneous Poisson process
  #       nsim: preset number of events
  # output  z: states of Markov process
  #         x: time of each transition of Markov process
  #         zt: state at each event
  #         tau: times of Poisson events
  if (!is.null(seed)) set.seed(seed)
  m <- 2
  #------------------------
  if (sum(object$delta)!=1) stop("Invalid delta")
  if (any(object$delta==1))
    initial <- (1:m)[as.logical(object$delta)]
  else
    initial <- sample(m, 1, prob=object$delta)
  #------------------------
  Q <- matrix(c(-object$q1,object$q1,object$q2,-object$q2),ncol=m,byrow=TRUE)
  lambda0 <- object$lambda0
  lambda1 <- object$lambda1
  alpha <- object$alpha
  beta <- object$beta
  
  Pi <- diag(m) - diag(1/diag(Q)) %*% Q
  zt <- rep(NA, nsim+1)
  tau <- rep(NA, nsim+1)
  #------------------------ initialization for Markov process
  #    the length of x and z may be too short
  #    gets extended later if required
  x <- rep(NA, nsim*10)
  z <- rep(NA, nsim*10)
  z[1] <- zt[1] <- initial
  x[1] <- tau[1] <- 0
  lambda.max <- 0
  i <- 1 #index for state
  j <- 2 #index for event
  #------------------------ initialization for Hawkes process
  
  while (j < nsim+2){
    i <- i+1
    #   extend x and z if too short
    if (i > length(x)){
      x <- c(x, rep(NA, nsim*10))
      z <- c(z, rep(NA, nsim*10))
    }
    #   sim time spent in Markov state y[i-1]
    z[i] <- sample(x=1:m, size=1, prob=Pi[(z[i-1]),])
    x[i] <- x[i-1] + rexp(1, rate=-Q[z[i-1], z[i-1]])
    t0 <- x[i-1]
    
    if(z[i-1]==1){
      #   sim times of Hawkes Poisson events
      simulate.result <- simulate.hp(lambda1,alpha,beta,x[i-1],x[i],tau[1:(j-1)])
      hp <- simulate.result$t
      lambda.max <- ifelse(lambda.max>simulate.result$lambda.max,lambda.max,simulate.result$lambda.max)
      if(!hp[1]==0){
        tau[j:(j+length(hp)-1)] <- hp
        zt[j:(j+length(hp)-1)] <- z[i-1]
        j <- j + length(hp)
      }
    }
    
    if(z[i-1]==2){
      while(j < nsim+2){
        #   sim times of Poisson events
        ti <- t0 + rexp(1, rate=lambda0)
        if (ti < x[i]){
          tau[j] <- t0 <- ti
          zt[j] <- z[i-1]
          j <- j + 1
        }
        else break
      }
    }
  }
  return(list(x=x[1:i],z=z[1:i],tau=tau[1:(nsim+1)],zt=zt[1:(nsim+1)],lambda.max=lambda.max))
}

##Simulate Hawkes process during active state (including all the histoty)
simulate.hp <- function(lambda0,alpha,beta,start,horizon,history){
  j0<-length(history)+1
  lambda.star<-ifelse(j0==2,lambda0,lambda0+alpha*sum(exp(-beta*(rep(start,j0-2)-history[2:(j0-1)]))))
  lambda.max <- lambda.star
  t<-numeric(10)
  n<-1
  U<-runif(1)
  s<--log(U)/lambda.star
  ti<-start+s
  repeat {
    if (ti > horizon){
      break
    } 
    
    lambda.star<-lambda.star+alpha
    t[n]<-ti
    if (length(t) < n+1) t <- c(t, numeric(10)) 
    
    repeat{
      U<-runif(1)
      s<-s-log(U)/lambda.star
      ti<-start+s
      lambda.s<-lambda0+alpha*sum(exp(-beta*c(rep(ti,n)-t[1:n],rep(ti,j0-1)-history[1:j0-1])))
      D<-runif(1)
      if(D<=lambda.s/lambda.star){
        lambda.star<-lambda.s
        lambda.max <- ifelse(lambda.max>lambda.star,lambda.max,lambda.star)
        break
      }
      lambda.star<-lambda.s
      lambda.max <- ifelse(lambda.max>lambda.star,lambda.max,lambda.star)
    }
    
    n <- n + 1
  }
  
  return(list(t=t[1:(n-1)],lambda.max=lambda.max))
}

##Simulate Markov Modulated Hawkes Process given hidden state trajectory
simulate.mmhp.given.state <- function(object, states, ending, max.nsim=1, seed=NULL, ...){
  # VALID FOR 2 STATES: 1: active, 2: not active
  # input q1,q2: transition probability
  #       delta: initial state probability; 
  #       lambda1, alpha, beta : parameters for Hawkes process
  #       lambda0 : parameters for homogeneous poisson process
  #       states: z: states of Markov process; x: time of each transition of Markov process
  #       max.nsim: preset max number of events
  #       ending: preset ending time for the process
  # output  zt: state at each event
  #         tau: times of Poisson events
  if (!is.null(seed)) set.seed(seed)
  m <- 2
  #------------------------
  Q <- matrix(c(-object$q1,object$q1,object$q2,-object$q2),ncol=m,byrow=TRUE)
  lambda0 <- object$lambda0
  lambda1 <- object$lambda1
  alpha <- object$alpha
  beta <- object$beta
  
  Pi <- diag(m) - diag(1/diag(Q)) %*% Q
  zt <- rep(NA, max.nsim+1)
  tau <- rep(NA, max.nsim+1)
  #------------------------ initialization for Markov process
  #    the length of x and z may be too short
  #    gets extended later if required
  x <- states$x
  z <- states$z
  zt[1] <- z[1] 
  tau[1] <- 0
  lambda.max <- 0
  i <- 1 #index for state
  j <- 2 #index for event
  #------------------------ initialization for Hawkes process
  
  while (tau[j-1] <= ending&i<length(x)){
    i <- i+1
    t0 <- x[i-1]
    
    if(z[i-1]==1){
      #   sim times of Hawkes Poisson events
      simulate.result <- simulate.hp(lambda1,alpha,beta,x[i-1],x[i],tau[1:(j-1)])
      hp <- simulate.result$t
      lambda.max <- ifelse(lambda.max>simulate.result$lambda.max,lambda.max,simulate.result$lambda.max)
      if(!hp[1]==0){
        tau[j:(j+length(hp)-1)] <- hp
        zt[j:(j+length(hp)-1)] <- z[i-1]
        j <- j + length(hp)
      }
    }
    
    if(z[i-1]==2){
      while(tau[j-1] <= ending){
        #   sim times of Poisson events
        ti <- t0 + rexp(1, rate=lambda0)
        if (ti < x[i]){
          tau[j] <- t0 <- ti
          zt[j] <- z[i-1]
          j <- j + 1
        }
        else break
      }
    }
  }
  return(list(tau=tau[1:(j-1)][tau[1:(j-1)]<=ending],zt=zt[1:(j-1)][tau[1:(j-1)]<=ending],lambda.max=lambda.max))
}

##This function plots the intensity function 
drawUniMMHPIntensity<-function(object,simulation,yupper=10,add=FALSE,color=1,given_main="Intensity Plot of MMHP"){
  # input object: the parameter list used for generating mmhp
  #       simulation: simulation result from simulate.mmhp 
  
  t <- simulation$tau
  state <- simulation$z
  state_time <- simulation$x
  lambda0 <- object$lambda0
  lambda1 <- object$lambda1
  alpha <- object$alpha
  beta <- object$beta
  
  n <- length(t)
  m <- length(state)
  
  if(add==FALSE){
    plot(0,0,xlim=c(0,state_time[m]),ylim=c(0,yupper),type="n",xlab="Time",ylab="Intensity",
         main=given_main)
    points(t[-1],rep(lambda0/2,n-1),cex=0.6,pch=ifelse(simulation$zt[-1]==1,16,1),col="blue")
    points(state_time,rep(lambda0,m),cex=0.6,pch=4,col="red")
  }
  for(i in 1:(m-1)){
    if(state[i]==1){
      hawkes_time <- t[t>=state_time[i]&t<state_time[i+1]]
      if(i==1) hawkes_time <- hawkes_time[-1]
      history <- t[t<state_time[i]]
      drawHPIntensity(lambda1,i,alpha,beta,state_time[i],state_time[i+1],history[-1],hawkes_time,color=color)
    }else{
      segments(x0=state_time[i],x1=state_time[i+1], y0=lambda0, lty=2,col=color)
    }
  }
  
  if(add==FALSE){
    legend("topleft",c("Hawkes event","Poisson process event","state change point"),col = c("blue","blue","red"), 
           pch = c(16,1,4))
  }else{
    legend("topright", c("True","Estimation"),col=c("black",color),lty=c(1,1))
  }
}

##This function is a helper function for 'drawUniMMHPIntensity'
drawHPIntensity<-function(lambda0,i,alpha,beta,start,end,history,hawkes_time,color=1){
  n <- length(hawkes_time)
  m <- length(history)
  
  if(n==0){
    if(i==1){
      segments(x0=start,x1=end,y0=lambda0)
    } else{
      lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,m)-history)))
      new.lambda.n<-Vectorize(lambda.n)
      segments(x0=start,y0=lambda0,y1=lambda.n(end),lty=2,col=color)
      curve(new.lambda.n, from=start, to=end, add=TRUE,col=color)
    }
  }else{
    if(i==1){
      segments(x0=start,x1=hawkes_time[1],y0=lambda0,col=color)
      segments(x0=hawkes_time[1],y0=lambda0,y1=lambda0+alpha,col=color)
    } else{
      lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,m)-history)))
      new.lambda.n<-Vectorize(lambda.n)
      segments(x0=start,y0=lambda0,y1=lambda.n(start),lty=2,col=color)
      curve(new.lambda.n, from=start, to=hawkes_time[1], add=TRUE,col=color)
      segments(x0=hawkes_time[1],y0=lambda.n(hawkes_time[1]),y1=lambda.n(hawkes_time[1])+alpha,col=color)
    }
    if(n>1){
      for(j in 1:(n-1)){
        lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,m+j)-c(history,hawkes_time[1:j]))))
        new.lambda.n<-Vectorize(lambda.n)
        curve(new.lambda.n, from=hawkes_time[j], to=hawkes_time[j+1], add=TRUE,col=color)
        segments(x0=hawkes_time[j+1],y0=lambda.n(hawkes_time[j+1]),y1=lambda.n(hawkes_time[j+1])+alpha,col=color)
      }
    }
    lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,m+n)-c(history,hawkes_time[1:n]))))
    new.lambda.n<-Vectorize(lambda.n)
    curve(new.lambda.n, from=hawkes_time[n], to=end, add=TRUE,col=color)
    segments(x0=end,y0=lambda.n(end),y1=lambda0,lty=2,col=color)
  }
}

##Compute negative log likelihood for mmhp
mmhpNegLogLik <- function(params=list(lambda0,lambda1,alpha,beta,q1,q2), t) {
  # t is event time, t[1]=0
  lambda0 <- params[[1]]
  lambda1 <- params[[2]]
  alpha <- params[[3]]
  beta <- params[[4]]
  q1 <- params[[5]]
  q2 <- params[[6]]
  n <- length(t) - 1
  interevent <- t[-1]-t[-(n+1)]
  
  forward <- matrix(0,ncol=2,nrow=n)
  probs_1 <- matrix(0,ncol=2,nrow=n) #Probability vector for transition to state 1 (active state)
  probs_2 <- matrix(0,ncol=2,nrow=n) #Probability vector for transition to state 2 (inactive state)
  r <- rep(0,n)
  
  integ1 <- interevent[1]*lambda1
  integ2 <- interevent[1]*lambda0
  
  probs_1[,1] = -q1*interevent
  probs_2[,2] = -q2*interevent
  probs_1[,2] = log(1-exp(probs_2[,2]))
  probs_2[,1] = log(1-exp(probs_1[,1]))
  
  forward[1,1] = log(lambda1) - integ1
  forward[1,2] = log(lambda0) - integ2
  
  for(i in 2:n){
    integ1 <- interevent[i]*lambda1
    integ2 <- interevent[i]*lambda0
    r[i] <- exp(-beta*interevent[i])*(r[i-1]+1)
    a <- min(forward[i-1,] + probs_1[i-1,])
    forward[i,1] <- a + log(sum(exp(forward[i-1,] + probs_1[i-1,] - a))) + log(lambda1+alpha*exp(-beta*interevent[i])*(r[i-1]+1)) - integ1 + 
      alpha/beta*(r[i]-r[i-1]-1)
    a <- min(forward[i-1,] + probs_2[i-1,])
    forward[i,2] <- a + log(sum(exp(forward[i-1,] + probs_2[i-1,] - a))) + log(lambda0) - integ2
  }
  return(-sum(forward[n,]))
}

##Compute vector of compensators for mmhp
mmhpCompensator<-function(params=list(lambda0,lambda1,alpha,beta,q1,q2), t, pzt, if.pzt=TRUE){
  lambda0 <- params$lambda0
  lambda1 <- params$lambda1
  alpha <- params$alpha
  beta <- params$beta
  q1 <- params$q1
  q2 <- params$q2
  n <- length(t) - 1
  interevent <- t[-1]-t[-(n+1)]
  if(if.pzt){
    pzt <- mmhpLocalState(params=params, interevent)$pzt
  }
  
  ##compute compensator for Hawkes process
  Lambda <- rep(0,n)
  A<-0
  Lambda[1]<-lambda0*(interevent[1])
  for(i in 2:n){
    A<-1+exp(-beta*(interevent[i-1]))*A
    Lambda[i]<-lambda1*(interevent[i])+alpha/beta*(1-exp(-beta*interevent[i]))*A
  }
  
  Lambda_mixed <- Lambda*pzt + lambda0*interevent*(1-pzt)
  return(Lambda_mixed)
}

##Infer the probability of zt=1 for each event using Viterbi algorithm
mmhpViterbi <- function(params=list(lambda0,lambda1,alpha,beta,q1,q2), interevent){
  lambda0 <- params$lambda0
  lambda1 <- params$lambda1
  alpha <- params$alpha
  beta <- params$beta
  q1 <- params$q1
  q2 <- params$q2
  N <- length(interevent)
  
  zt <- rep(0,N)
  r <- rep(0,N) 
  mm <- matrix(0,nrow=N,ncol=2)
  bb <- matrix(0,nrow=N,ncol=2)
  probs_1 <- matrix(0,ncol=2,nrow=N) #Probability vector for transition to state 1 (active state)
  probs_2 <- matrix(0,ncol=2,nrow=N) #Probability vector for transition to state 2 (inactive state)
  
  integ1 <- interevent[1]*lambda1
  integ2 <- interevent[1]*lambda0
  
  probs_1[,1] <- -q1*interevent
  probs_2[,2] <- -q2*interevent
  probs_1[,2] <- log(1-exp(probs_2[,2]))
  probs_2[,1] <- log(1-exp(probs_1[,1]))
  
  #calculate m and log_p_z_star => zt*[N]
  integral <- interevent[1]*lambda0
  mm[1,1] <- log(lambda1) - interevent[1]*lambda1
  mm[1,2] <- log(lambda0) - integral #calculate forward variables, uniform initial distribution for latent state
  r[1] <- 0
  for(n in 2:N){
    r[n] <- exp(-beta*interevent[n])*(r[n-1]+1)
    
    logp1 <- mm[n-1,1] + probs_1[n,1]
    logp2 <- mm[n-1,2] + probs_1[n,2]
    if(logp1>logp2){
      mm[n,1] <- logp1 + log(max(lambda1+alpha*exp(-beta*interevent[n])*(r[n-1]+1),0)) - interevent[n]*lambda1 + alpha/beta*(r[n]-r[n-1]-1)
      bb[n,1] <- 1
    }else{
      mm[n,1] <- logp2 + log(max(lambda1+alpha*exp(-beta*interevent[n])*(r[n-1]+1),0)) - interevent[n]*lambda1 + alpha/beta*(r[n]-r[n-1]-1)
      bb[n,1] <- 2
    }
    
    logp1 <- mm[n-1,1] + probs_2[n,1]
    logp2 <- mm[n-1,2] + probs_2[n,2]
    if(logp1>logp2){
      mm[n,2] <- logp1 + log(lambda0) - interevent[n]*lambda0
      bb[n,2] <- 1
    }else{
      mm[n,2] <- logp2 + log(lambda0) - interevent[n]*lambda0
      bb[n,2] <- 2
    }
  }
  
  if(mm[N,1]>mm[N,2]){
    zt[N] <- 1
  }else{
    zt[N] <- 2
  }
  
  for(n in 1:(N-1)){
    zt[N-n] <- bb[N-n+1,zt[N-n+1]]
  }
  return(zt)
}

##Infer the probability of zt=1 for each event using local decoding
mmhpLocalState <- function(params=list(lambda0,lambda1,alpha,beta,q1,q2), interevent){
  lambda0 <- params$lambda0
  lambda1 <- params$lambda1
  alpha <- params$alpha
  beta <- params$beta
  q1 <- params$q1
  q2 <- params$q2
  N <- length(interevent)
  
  pzt <- rep(0,N)
  zt <- rep(0,N)
  r <- rep(0,N) 
  intensities <- rep(0,2)
  integrals <- rep(0,2)
  forward <- matrix(0,nrow=N,ncol=2)
  backward <- matrix(0,nrow=N,ncol=2)
  probs_1 <- matrix(0,nrow=N,ncol=2) #Probability vector for transition to state 1 (active state)
  probs_2 <- matrix(0,nrow=N,ncol=2) #Probability vector for transition to state 2 (inactive state)
  
  probs_1[,1] <- -q1*interevent
  probs_2[,2] <- -q2*interevent
  probs_1[,2] <- log(1-exp(probs_2[,2]))
  probs_2[,1] <- log(1-exp(probs_1[,1]))
  
  #calculate forward and log_p_z_star => zt*[N]
  forward[1,1] <- log(lambda0) - interevent[1]*lambda1
  forward[1,2] <- log(lambda0) - interevent[1]*lambda0 #calculate forward variables, uniform initial distribution for latent state
  r[1] <- 0
  
  for(n in 2:N){
    r[n] <- exp(-beta*interevent[n])*(r[n-1]+1)
    a <- max(forward[n-1,] + probs_1[n,])
    forward[n,1] <- a + log(sum(exp(forward[n-1,] + probs_1[n,] - a))) + log(lambda1+alpha*exp(-beta*interevent[n])*(r[n-1]+1)) - interevent[n]*lambda1 + 
      alpha/beta*(r[n]-r[n-1]-1)
    a <- max(forward[n-1,] + probs_2[n-1,])
    forward[n,2] <- a + log(sum(exp(forward[n-1,] + probs_2[n,] - a))) + log(lambda0) - interevent[n]*lambda0
  }
  
  #calculate backward variables
  backward[N,1] <- 0
  backward[N,2] <- 0
  intensities[2] <- lambda0
  
  for(n in 2:N){
    m <- N-n+1
    intensities[1] <- lambda1 + alpha*exp(-beta*interevent[m+1])*(r[m]+1)
    integrals[2] <- lambda0*interevent[m]
    integrals[1] <- lambda1*interevent[m] - alpha/beta*(r[m+1]-r[m]-1)
    a <-  max(backward[m+1,]+probs_1[m+1,]+log(intensities)-integrals)
    backward[m,1] <- a+log(sum(exp(backward[m+1,]+probs_1[m+1,]+log(intensities)-integrals-a)))
    a <-  max(backward[m+1,]+probs_2[m+1,]+log(intensities)-integrals)
    backward[m,2] <- a+log(sum(exp(backward[m+1,]+probs_2[m+1,]+log(intensities)-integrals-a)))
  }
  
  #infer the probability of zt=1
  for(n in 1:N){
    pzt[n] <- 1/(1+exp(forward[n,2]-forward[n,1]+backward[n,2]-backward[n,1]))
    if(forward[n,2] + backward[n,2] > forward[n,1] + backward[n,1]){
      zt[n] = 2;
    }else{
      zt[n] = 1;
    }
  }
  return(list(pzt=pzt,zt=zt))
}

##Interpolate the latent state trajectory given estimated states of events
modifiedLatentTrajectory <- function(params=list(lambda0,lambda1,alpha,beta,q1,q2), interevent, zt, start=0){
  # input pars
  #       interevent: length N
  #       zt: inferred latent state, 1 or 2 , length N
  # output  z.hat: states of Markov process
  #         x.hat: time of each transition of Markov process
  lambda0 <- params$lambda0
  lambda1 <- params$lambda1
  alpha <- params$alpha
  beta <- params$beta
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
        #x.hat[temp.count] <- temp.t[l]
        #z.hat[temp.count] <- 1
        if((alpha<0)|(q2-q1<0)){
          x.hat[temp.count] <- temp.t[l-1]
          z.hat[temp.count] <- 1
        }else{
          temp.A <- alpha/beta*sum(exp(beta*temp.t[1:(l-1)]))
          temp.delta.t <- -temp.t[l-1]-log((q2-q1)/(temp.A*beta))/beta
          x.hat[temp.count] <- ifelse(temp.delta.t>0&(temp.delta.t<temp.t[l]-temp.t[l-1]),temp.t[l-1] + temp.delta.t,
                                      ifelse(temp.delta.t<=0,temp.t[l-1],temp.t[l]))
          z.hat[temp.count] <- 1
        }
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
