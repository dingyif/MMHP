drawUniMMHPIntensityPaper<-function(object,simulation,yupper=10,add=FALSE,color=1,line.width=1,title_name=""){
  # input object: the parameter list used for generating mmhp
  #       simulation: simulation result from simulate.mmhp 
  
  t <- simulation$tau
  state <- simulation$z
  state_time <- simulation$x
  lambda0 <- object$lambda0
  alpha <- object$alpha
  beta <- object$beta
  
  n <- length(t)
  m <- length(state)
  
  #plot(1,4,xlim=c(0.27,23),ylim=c(-0.5,9), type="n",xlab="Time",ylab="",xaxt="n",yaxt="n",cex.lab=1.8,bty="n")
  
  if(add==FALSE){
    plot(0,0,xlim=c(0.4,state_time[m]),ylim=c(-1,yupper),type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="l",cex.lab=2)
    title(title_name,line=1,cex.main=9)
    axis(2,at=seq(0,30,5),labels=FALSE,cex.lab=2,lwd=0,lwd.ticks=1)
    axis(1,at=seq(0,12,4),labels=FALSE,cex.lab=2,lwd=0,lwd.ticks=1)
    text(-0.8, y=seq(0,30,5), labels=seq(0,30,5), cex=1.6, srt=0, xpd=TRUE)
    text(seq(0,12,4), y=-3.4, labels=seq(0,12,4), cex=1.7, srt=0, xpd=TRUE)
    points(t[-1],rep(-0.5,n-1),cex=1.3,pch=ifelse(simulation$zt[-1]==1,16,1),col="blue")
    points(state_time,rep(lambda0,m),cex=1.5,pch=4,col="red",lwd=2.5)
  }
  for(i in 1:(m-1)){
    if(state[i]==1){
      hawkes_time <- t[t>=state_time[i]&t<state_time[i+1]]
      if(i==1) hawkes_time <- hawkes_time[-1]
      history <- t[t<state_time[i]]
      drawHPIntensityPaper(lambda0,i,alpha,beta,state_time[i],state_time[i+1],history[-1],hawkes_time,color=color,line.width=line.width)
    }else{
      segments(x0=state_time[i],x1=state_time[i+1], y0=lambda0, lty=1,col=color,lwd=line.width)
    }
  }
  
  # if(add==FALSE){
  #   legend(9,32,c("State 1/0 events","State change point","True intensity","Estimated MMHP intensity",
  #                  "Estimated MMHPSD intensity","Estimated MMPP intensity"),
  #          col = c(NA,"red","black",rgb(0,190,0,maxColorValue = 255),
  #                  rgb(221,160,221,maxColorValue = 255),rgb(169,169,169,maxColorValue = 255)),
  #          y.intersp=0.88,x.intersp=0.15,
  #          pch = c(NA,4,NA,NA,NA,NA), pt.cex = c(NA,1.5,NA,NA,NA,NA),
  #          lty = c(NA,NA,1,1,1,1), lwd = c(NA,2.5,5,5,5,5),cex=1.6)
  #   points(9.43,30.09,pch=1,cex=1.8,col="blue")
  #   points(9.2,30.09,pch=16,cex=1.9,col="blue")
  # }
}

drawHPIntensityPaper<-function(lambda0,i,alpha,beta,start,end,history,hawkes_time,color=1,line.width=1){
  n <- length(hawkes_time)
  m <- length(history)
  
  if(n==0){
    if(i==1){
      segments(x0=start,x1=end,y0=lambda0,lwd=line.width)
    } else{
      lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,m)-history)))
      new.lambda.n<-Vectorize(lambda.n)
      segments(x0=start,y0=lambda0,y1=lambda.n(end),lty=2,col=color,lwd=line.width)
      curve(new.lambda.n, from=start, to=end, add=TRUE,col=color,lwd=line.width)
    }
  }else{
    if(i==1){
      segments(x0=start,x1=hawkes_time[1],y0=lambda0,col=color,lwd=line.width)
      segments(x0=hawkes_time[1],y0=lambda0,y1=lambda0+alpha,col=color,lwd=line.width)
    } else{
      lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,m)-history)))
      new.lambda.n<-Vectorize(lambda.n)
      segments(x0=start,y0=lambda0,y1=lambda.n(start),lty=2,col=color,lwd=line.width)
      curve(new.lambda.n, from=start, to=hawkes_time[1], add=TRUE,col=color,lwd=line.width)
      segments(x0=hawkes_time[1],y0=lambda.n(hawkes_time[1]),y1=lambda.n(hawkes_time[1])+alpha,col=color,lwd=line.width)
    }
    if(n>1){
      for(j in 1:(n-1)){
        lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,m+j)-c(history,hawkes_time[1:j]))))
        new.lambda.n<-Vectorize(lambda.n)
        curve(new.lambda.n, from=hawkes_time[j], to=hawkes_time[j+1], add=TRUE,col=color,lwd=line.width)
        segments(x0=hawkes_time[j+1],y0=lambda.n(hawkes_time[j+1]),y1=lambda.n(hawkes_time[j+1])+alpha,col=color,lwd=line.width)
      }
    }
    lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,m+n)-c(history,hawkes_time[1:n]))))
    new.lambda.n<-Vectorize(lambda.n)
    curve(new.lambda.n, from=hawkes_time[n], to=end, add=TRUE,col=color,lwd=line.width)
    segments(x0=end,y0=lambda.n(end),y1=lambda0,lty=2,col=color,lwd=line.width)
  }
}

plotIntensityThreeModel <- function(data.path, no.events=100, plot.idx=16, my.mar=c(0,0,1,1)+0.1){
  par(mfrow=c(3,1),tcl=0.2,mgp=c(0.5,0,0),oma=c(3,3.2,2.5,0)+0.5,
      mar=my.mar,xpd=TRUE)
  
  ### 1. mmhp
  load(paste(data.path,"fixed_state_stan_result_",no.events,".Rdata",sep=""))
  mmhp_par_name <- c("lambda0","alpha","beta","q1","q2")
  mmhp_par_est <- numeric(length(mmhp_par_name))
  test.mmhp <- fixed.state.processes[plot.idx][[1]]
  test.mmhp$lambda.max <- 35
  for(j in 1:5){
    mmhp_par_est[j] <- mean(unlist(fixed.state.stan.result[plot.idx][[1]][mmhp_par_name[j]])[1001:2000])
  }
  object_hat <- list(lambda0=mmhp_par_est[1], alpha=mmhp_par_est[2], beta=mmhp_par_est[3], q1=mmhp_par_est[4], q2=mmhp_par_est[5])
  test.latent <- stepfun(test.mmhp.fix$x[-1],test.mmhp.fix$z)
  time.segment <- seq(test.mmhp$tau[2],tail(test.mmhp$tau,1),length.out=2000)
  time.segment.delta <- time.segment[2]-time.segment[1]
  
  fixed.state.est.latent.mmhp <- ifelse(apply(fixed.state.stan.result[plot.idx][[1]]$zt[1001:2000,],2,function(x) sum(x==1)) >500,rep(1,1000),rep(2,1000))
  fixed.state.est.latent.mmhp <- modifiedLatentTrajectory(params=object_hat, interevent = test.mmhp$tau[-1]-test.mmhp$tau[-length(test.mmhp$tau)], zt = fixed.state.est.latent.mmhp)
  
  drawUniMMHPIntensityPaper(object,simulation=list(x=test.mmhp.fix$x,z=test.mmhp.fix$z,tau=test.mmhp$tau,zt=test.mmhp$zt,lambda.max=test.mmhp$lambda.max),yupper=test.mmhp$lambda.max+2, color="black",line.width=3.5,
                            title_name="")
  drawUniMMHPIntensityPaper(object_hat,simulation=list(x=fixed.state.est.latent.mmhp$x.hat,z=fixed.state.est.latent.mmhp$z.hat,tau=test.mmhp$tau,zt=test.mmhp$zt,lambda.max=test.mmhp$lambda.max),
                            yupper=test.mmhp$lambda.max+1,add=TRUE,color=rgb(0,190,0,alpha=180,maxColorValue=255),line.width=4.2)
  
  legend(-0.31,44.6,c("State 1/0 events","State change point"),
         col = c(NA,"red"),
         y.intersp=0.85,x.intersp=0.01,bty = "n",
         pch = c(NA,4), pt.cex = c(NA,3), cex=4,
         lty = c(NA,NA), lwd=c(NA,4))
  legend(5.65,44,c("True intensity","Estimated MMHP intensity"),
         col = c("black",rgb(0,190,0,maxColorValue = 255)),
         y.intersp=0.85,x.intersp=0.15,bty = "n",
         lty = c(1,1), lwd = c(5,5), cex=4)
  points(0.3,33.6,pch=1,cex=3,col="blue")
  points(0.1,33.6,pch=16,cex=3,col="blue")
  
  ### 2. mmhpsd
  load(paste(data.path,"fixed_state_mmhpsd_result_",no.events,".Rdata",sep=""))
  event_time <- fixed.state.processes[plot.idx][[1]]$tau
  time_horizon <- tail(test.mmhp$tau,1)
  time_segment <- 0.02
  intensity_estimation <- numeric(1+round(time_horizon/time_segment))
  state_estimation <- numeric(1+round(time_horizon/time_segment))
  
  drawUniMMHPIntensityPaper(object,simulation=list(x=test.mmhp.fix$x,z=test.mmhp.fix$z,tau=test.mmhp$tau,zt=test.mmhp$zt,lambda.max=test.mmhp$lambda.max),yupper=test.mmhp$lambda.max+2, color="black",line.width=3.5)
  
  for(j in 0:(round(time_horizon/time_segment))){
    est_result <- tryCatch(estSInt(tims = j*time_segment,
                                   ti = event_time,
                                   lamb0 = mmhpsd_par_result[plot.idx][[1]]$lamb,
                                   nu0 = mmhpsd_par_result[plot.idx][[1]]$nu,
                                   eta0 = mmhpsd_par_result[plot.idx][[1]]$eta,
                                   Q0 = mmhpsd_par_result[plot.idx][[1]]$Q,
                                   pai0 = mmhpsd_par_result[plot.idx][[1]]$pai, fortran = TRUE),
                           warning = function(w) {0},
                           error = function(e) {0})
    intensity_estimation[j+1] <- ifelse(class(est_result)=="list",est_result$lambdat,est_result)
    state_estimation[j+1] <- ifelse(class(est_result)=="list",est_result$estStat[2],est_result)
  }
  
  state_estimation <- state_estimation>0.5
  change_point <- numeric(0)
  for(j in c(1:length(state_estimation[-1]))){
    if(state_estimation[j]==FALSE&state_estimation[j+1]==TRUE){
      change_point <- c(change_point,j+1)
    }
    if(state_estimation[j]==TRUE&state_estimation[j+1]==FALSE){
      change_point <- c(change_point,j)
    }
  }
  time_vec <- c(0:(round(time_horizon/time_segment)))*time_segment
  j <- 1
  line_type <- 1
  for(k in c(1:length(change_point))){
    lines(x=time_vec[c(j:change_point[k])], y=intensity_estimation[c(j:change_point[k])],
          lwd=4.2, col=rgb(218,112,214,alpha=180,maxColorValue=255),lty=line_type)
    j <- change_point[k]
  }
  lines(x=time_vec[c(j:length(state_estimation[-1]))], y=intensity_estimation[c(j:length(state_estimation[-1]))],
        lwd=4.2, col=rgb(218,112,214,alpha=180,maxColorValue=255),lty=line_type)
  legend(5.65,42,c("Estimated MMHPSD intensity"),
         col = c(rgb(221,160,221,maxColorValue = 255)),
         y.intersp=0.88,x.intersp=0.15,bty = "n",
         lty = c(1,1), lwd = c(5,5), cex=4,seg.len=1.5)
  
  ### 3.mmpp
  load(paste(data.path,"mmpp_fixed_state_stan_result_",no.events,".Rdata",sep=""))
  mmpp_par_name <- c("lambda0","c","q1","q2")
  mmpp_par_est <- numeric(4)
  fixed.state.latent.fun <- matrix(rep(list(),50*2),nrow=50,ncol=2)
  
  drawUniMMHPIntensityPaper(object,simulation=list(x=test.mmhp.fix$x,z=test.mmhp.fix$z,tau=test.mmhp$tau,zt=test.mmhp$zt,lambda.max=test.mmhp$lambda.max),yupper=test.mmhp$lambda.max+2, color="black",line.width=3.5)
  
  for(j in 1:4){
    mmpp_par_est[j] <- mean(unlist(fixed.state.stan.result[plot.idx][[1]][mmpp_par_name[j]])[1001:2000])
  }
  
  object_hat <- list(lambda0=mmpp_par_est[1], c=mmpp_par_est[2], q1=mmpp_par_est[3],
                     q2=mmpp_par_est[4])
  fixed.state.est.latent <- ifelse(apply(fixed.state.stan.result[plot.idx][[1]]$zt_v[1001:2000,],2,function(x) sum(x==1)) >500,rep(1,1000),rep(2,1000))
  fixed.state.est.latent <- mmppModifiedLatentTrajectory(params=object_hat, interevent = test.mmhp$tau[-1]-test.mmhp$tau[-length(test.mmhp$tau)], zt = fixed.state.est.latent)
  fixed.state.latent.fun <- stepfun(fixed.state.est.latent$x.hat[-1],fixed.state.est.latent$z.hat)
  
  time.segment <- seq(0,tail(test.mmhp$tau,1),length.out = 5000)
  fun.result <- fixed.state.latent.fun(time.segment)
  lines(ifelse(fun.result==2,object_hat$lambda0,object_hat$lambda0*(1+object_hat$c))~time.segment, type="l",col=rgb(169,169,169,alpha=180,maxColorValue=255), lwd=4.2)
  
  legend(5.65,42,c("Estimated MMPP intensity"),
         col = c(rgb(169,169,169,alpha=180,maxColorValue=255)),
         y.intersp=0.88,x.intersp=0.15,bty = "n",
         lty = c(1,1), lwd = c(5,5),cex=4)
  mtext(text="(b)",side=3,line=-1,outer=TRUE,cex=3,font=2)
  mtext(text="Time",side=1,line=2.1,outer=TRUE,cex=2.5)
  mtext(text="Intensity",side=2,line=0.4,outer=TRUE,cex=2.5)
}
