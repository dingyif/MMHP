plotSimMMHP <- function(data.path, my.mar=c(3.5,6.7,3,0.8)){
  par(mfrow=c(1,1),tcl=0.2,mgp=c(1.7,0,0),
      mar=my.mar)
  plot(1,4,xlim=c(0.27,23),ylim=c(-0.5,9), type="n",xlab="Time",ylab="",xaxt="n",yaxt="n",cex.lab=2.5,bty="n")
  title("(a)",line=1,cex.main=3)
  axis(2,at=c(0,1,2.5,3.5,5,6,7.5,8.5),labels=FALSE,cex.lab=2,las=2,lwd=0,lwd.ticks=1)
  axis(1,at=seq(0,20,length.out=5),labels=FALSE,cex.lab=2,las=2,lwd=0,lwd.ticks=1,line=-0.7)
  text(-1.8, y=c(0,1,2.5,3.5,5,6,7.5,8.5), labels=paste(rep(c("        ","state"),4),rep(c(0,1),4)," "), cex=2.5, srt=0, xpd=TRUE)
  text(seq(0,20,length.out=5), y=-1, labels=seq(0,20,length.out=5), cex=2, srt=0, xpd=TRUE)
  
  for(l in c(1:4)){
    no_events <- no_vec_events[l]
    load(paste(data.path,"fixed_state_stan_result_",no_events,".Rdata",sep=""))
    fixed.state.array <- matrix(NA,nrow=n_sim,ncol=5)
    fixed.state.est.latent <- matrix(rep(list(),n_sim*3),nrow=n_sim,ncol=3)
    fixed.state.latent.fun <- matrix(rep(list(),n_sim*2),nrow=n_sim,ncol=2)
    temp_par_name <- c("lambda0","alpha","beta","q1","q2") 
    
    for(i in 1:n_sim){
      for(j in c(1:5)){
        fixed.state.array[i,j] <- mean(unlist(fixed.state.stan.result[i][[1]][temp_par_name[j]])[1001:2000])
      }
      #estimate latent state
      test.mmhp <- fixed.state.processes[i][[1]]
      object_hat <- list(lambda0=fixed.state.array[i,1], lambda1=fixed.state.array[i,1], alpha=fixed.state.array[i,2], beta=fixed.state.array[i,3], q1=fixed.state.array[i,4], q2=fixed.state.array[i,5])
      test.latent <- stepfun(test.mmhp.fix$x[-1],test.mmhp.fix$z)
      time.segment <- seq(test.mmhp$tau[2],tail(test.mmhp$tau,1),length.out=2000)
      time.segment.delta <- time.segment[2]-time.segment[1]
      
      #2 estimation methods
      #fixed.state.est.latent[i,1][[1]] <- apply(fixed.state.stan.result[i][[1]]$zt_v[1001:2000,],2,mean)
      fixed.state.est.latent[i,2][[1]] <- ifelse(apply(fixed.state.stan.result[i][[1]]$zt_v[1001:2000,],2,function(x) sum(x==1)) >500,rep(1,1000),rep(2,1000))
      
      #preperation for step function
      fixed.state.est.latent[i,3][[1]] <- modifiedLatentTrajectory(params=object_hat, interevent = test.mmhp$tau[-1]-test.mmhp$tau[-length(test.mmhp$tau)], zt = fixed.state.est.latent[i,2][[1]])
      
      #interpolated function for each method
      fixed.state.latent.fun[i,1][[1]] <- approxfun(x=test.mmhp$tau[-1], y=fixed.state.est.latent[i,1][[1]],rule=2)
      fixed.state.latent.fun[i,2][[1]] <- stepfun(fixed.state.est.latent[i,3][[1]]$x.hat[-1],fixed.state.est.latent[i,3][[1]]$z.hat)
    }
    
    #plot
    time.segment <- seq(0,tail(test.mmhp.fix$tau,1),length.out = 5000)
    draw.start<-1
    draw.end<-length(test.mmhp.fix$tau)
    draw.true.x<-which(test.mmhp.fix$x>test.mmhp.fix$tau[draw.start]&test.mmhp.fix$x<test.mmhp.fix$tau[draw.end])
    draw.true.time<-c(test.mmhp.fix$tau[draw.start],test.mmhp.fix$x[draw.true.x],test.mmhp.fix$tau[draw.end])
    draw.true.state<-3-c(test.mmhp.fix$z[head(draw.true.x,1)-1],test.mmhp.fix$z[draw.true.x],test.mmhp.fix$z[tail(draw.true.x,1)])
    
    fun.result <- array(NA,dim=c(n_sim,5000,2))
    for(i in 1:n_sim){
      for(j in 1:2){
        fun.result[i,,j] <- fixed.state.latent.fun[i,j][[1]](time.segment)
      }
    }
    
    change <- 9-l*2.5
    cur_method <- 2
    for(i in 1:n_sim){
      state_mmhp <- fun.result[i,,2] 
      lines(3-state_mmhp+change~time.segment, type="l", 
            col=rgb(169,169,169,maxColorValue=255), lwd=0.35)
    }
    
    lower.ci <- apply(fun.result[,,cur_method],2,mean)-apply(fun.result[,,2],2,sd)
    lower.ci <- ifelse(lower.ci<1,rep(1,length(lower.ci)),lower.ci)
    upper.ci <- apply(fun.result[,,2],2,mean)+apply(fun.result[,,2],2,sd)
    upper.ci <- ifelse(upper.ci>2,rep(2,length(upper.ci)),upper.ci)
    lines(3-apply(fun.result[,,2],2,mean)+change~time.segment,type="l",col="blue",lwd=4)
    polygon(c(time.segment,rev(time.segment)),
            c(3-lower.ci+change,rev(3-upper.ci+change)),
            col=rgb(0,0,225,alpha=50,maxColorValue=255),border=FALSE)
    
    lines(draw.true.state+change~draw.true.time,type="s",col="grey27",lwd=4.5)
    segments(x0=-0.5,y0=7-(l-1)*2.5,y1=9-(l-1)*2.5)
    segments(x0=c(22.5,16,14.3,11)[5-l]+1.25,y0=7-(l-1)*2.5,y1=9-(l-1)*2.5)
    segments(x0=-0.5,x1=c(22.5,16,14.3,11)[5-l]+1.25,y0=7-(l-1)*2.5)
    segments(x0=-0.5,x1=c(22.5,16,14.3,11)[5-l]+1.25,y0=9-(l-1)*2.5)
  }
  
  legend(15,10,c("truth","MMHP mean","confidence band","individual process"),
         lty=1, lwd=c(5,5,20,5), x.intersp=0.3,
         col=c(rgb(0,0,0),rgb(0,0,225,maxColorValue=255),
               rgb(0,0,225,alpha=100,maxColorValue=255),rgb(169,169,169,maxColorValue=255)),
         bg=rgb(0,0,0,alpha=0), bty="n",cex=2.75,seg.len=0.85)
  text(c(22.5,15,12.5,9.5),c(0.5,3,5.5,8),paste("M=",rev(no_vec_events),sep=""),cex=2.5)
  
}

plotSimMMPP <- function(data.path, my.mar=c(3.5,6.7,3,0.8)){
  par(mfrow=c(1,1),tcl=0.2,mgp=c(1.7,0,0),
      mar=my.mar)
  plot(1,4,xlim=c(0.27,23),ylim=c(-0.5,9), type="n",xlab="Time",ylab="",xaxt="n",yaxt="n",cex.lab=2.5,bty="n")
  title("(c)",line=1,cex.main=3)
  axis(2,at=c(0,1,2.5,3.5,5,6,7.5,8.5),labels=FALSE,cex.lab=2,las=2,lwd=0,lwd.ticks=1)
  axis(1,at=seq(0,20,length.out=5),labels=FALSE,cex.lab=2,las=2,lwd=0,lwd.ticks=1,line=-0.7)
  text(-1.8, y=c(0,1,2.5,3.5,5,6,7.5,8.5), labels=paste(rep(c("        ","state"),4),rep(c(0,1),4)," "), cex=2.5, srt=0, xpd=TRUE)
  text(seq(0,20,length.out=5), y=-1, labels=seq(0,20,length.out=5), cex=2, srt=0, xpd=TRUE)
  mmpp_par_name <- c("lambda0","c","q1","q2")
  
  for(l in c(1:4)){
    no_events <- no_vec_events[l]
    load(paste(data.path,"mmpp_fixed_state_stan_result_",no_events,".Rdata",sep=""))
    fixed.state.array <- matrix(NA,nrow=50,ncol=4)
    fixed.state.est.latent <- matrix(rep(list(),50*3),nrow=50,ncol=3)
    fixed.state.latent.fun <- matrix(rep(list(),50*2),nrow=50,ncol=2)
    
    for(i in 1:50){
      for(j in 1:4){
        fixed.state.array[i,j] <- mean(unlist(fixed.state.stan.result[i][[1]][mmpp_par_name[j]])[1001:2000])
      }
      
      #estimate latent state
      test.mmhp <- fixed.state.processes[i][[1]]
      object_hat <- list(lambda0=fixed.state.array[i,1], c=fixed.state.array[i,2], q1=fixed.state.array[i,3], q2=fixed.state.array[i,4])
      test.latent <- stepfun(test.mmhp.fix$x[-1],test.mmhp.fix$z)
      time.segment <- seq(test.mmhp$tau[2],tail(test.mmhp$tau,1),length.out=2000)
      time.segment.delta <- time.segment[2]-time.segment[1]
      
      #2 estimation methods
      fixed.state.est.latent[i,1][[1]] <- apply(fixed.state.stan.result[i][[1]]$zt_v[1001:2000,],2,mean)
      fixed.state.est.latent[i,2][[1]] <- ifelse(apply(fixed.state.stan.result[i][[1]]$zt_v[1001:2000,],2,function(x) sum(x==1)) >500,rep(1,1000),rep(2,1000))
      
      #preperation for step function
      fixed.state.est.latent[i,3][[1]] <- mmppModifiedLatentTrajectory(params=object_hat, interevent = test.mmhp$tau[-1]-test.mmhp$tau[-length(test.mmhp$tau)], zt = fixed.state.est.latent[i,2][[1]])
      
      #interpolated function for each method
      fixed.state.latent.fun[i,1][[1]] <- approxfun(x=test.mmhp$tau[-1], y=fixed.state.est.latent[i,1][[1]],rule=2)
      fixed.state.latent.fun[i,2][[1]] <- stepfun(fixed.state.est.latent[i,3][[1]]$x.hat[-1],fixed.state.est.latent[i,3][[1]]$z.hat)
    }
    
    #plot
    time.segment <- seq(0,tail(test.mmhp.fix$tau,1),length.out = 5000)
    draw.start<-1
    draw.end<-length(test.mmhp.fix$tau)
    draw.true.x<-which(test.mmhp.fix$x>test.mmhp.fix$tau[draw.start]&test.mmhp.fix$x<test.mmhp.fix$tau[draw.end])
    draw.true.time<-c(test.mmhp.fix$tau[draw.start],test.mmhp.fix$x[draw.true.x],test.mmhp.fix$tau[draw.end])
    draw.true.state<-3-c(test.mmhp.fix$z[head(draw.true.x,1)-1],test.mmhp.fix$z[draw.true.x],test.mmhp.fix$z[tail(draw.true.x,1)])
    
    fun.result <- array(NA,dim=c(50,5000,2))
    for(i in 1:50){
      for(j in 1:2){
        fun.result[i,,j] <- fixed.state.latent.fun[i,j][[1]](time.segment)
      }
    }
    
    change <- 9-l*2.5
    
    for(i in 1:50){
      state_mmhp <- fun.result[i,,2] #2-(state_estimation<0.5)
      lines(3-state_mmhp+change~time.segment, type="l", 
            col=rgb(169,169,169,maxColorValue=255), lwd=0.35)
    }
    
    lower.ci <- apply(fun.result[,,2],2,mean)-apply(fun.result[,,2],2,sd)
    lower.ci <- ifelse(lower.ci<1,rep(1,length(lower.ci)),lower.ci)
    upper.ci <- apply(fun.result[,,2],2,mean)+apply(fun.result[,,2],2,sd)
    upper.ci <- ifelse(upper.ci>2,rep(2,length(upper.ci)),upper.ci)
    lines(3-apply(fun.result[,,2],2,mean)+change~time.segment,type="l",col="blue",lwd=4)
    polygon(c(time.segment,rev(time.segment)),
            c(3-lower.ci+change,rev(3-upper.ci+change)),
            col=rgb(0,0,225,alpha=50,maxColorValue=255),border=FALSE)
    
    lines(draw.true.state+change~draw.true.time,type="s",col="grey27",lwd=4.5)
    segments(x0=-0.5,y0=7-(l-1)*2.5,y1=9-(l-1)*2.5)
    segments(x0=c(22.5,16,14.3,11)[5-l]+1.25,y0=7-(l-1)*2.5,y1=9-(l-1)*2.5)
    segments(x0=-0.5,x1=c(22.5,16,14.3,11)[5-l]+1.25,y0=7-(l-1)*2.5)
    segments(x0=-0.5,x1=c(22.5,16,14.3,11)[5-l]+1.25,y0=9-(l-1)*2.5)
  }
  
  legend(15,10,c("truth","MMPP mean","confidence band","individual process"),
         lty=1, lwd=c(5,5,20,5), x.intersp=0.3,
         col=c(rgb(0,0,0),rgb(0,0,225,maxColorValue=255),
               rgb(0,0,225,alpha=100,maxColorValue=255),rgb(169,169,169,maxColorValue=255)),
         bg=rgb(0,0,0,alpha=0), bty="n",cex=2.75,seg.len=0.85)
  text(c(22.5,15,12.5,9.5),c(0.5,3,5.5,8),paste("M=",rev(no_vec_events),sep=""),cex=2.5)
 
}

plotSimMMHPSD <- function(data.path, my.mar=c(3.5,6.7,3,0.8)){
  par(mfrow=c(1,1),tcl=0.2,mgp=c(1.7,0,0),
      mar=my.mar)
  
  plot(1,4,xlim=c(0.27,23),ylim=c(-0.5,9), type="n",xlab="Time",ylab="",xaxt="n",yaxt="n",cex.lab=2.5,bty="n")
  title("(d)",line=1,cex.main=3)
  axis(2,at=c(0,1,2.5,3.5,5,6,7.5,8.5),labels=FALSE,cex.lab=2,las=2,lwd=0,lwd.ticks=1)
  axis(1,at=seq(0,20,length.out=5),labels=FALSE,cex.lab=2,las=2,lwd=0,lwd.ticks=1,line=-0.7)
  text(-1.8, y=c(0,1,2.5,3.5,5,6,7.5,8.5), labels=paste(rep(c("        ","state"),4),rep(c(0,1),4)," "), cex=2.5, srt=0, xpd=TRUE)
  text(seq(0,20,length.out=5), y=-1, labels=seq(0,20,length.out=5), cex=2, srt=0, xpd=TRUE)
  
  
  for(l in c(1:4)){
    no_events <- no_vec_events[l]
    load(paste(data.path,"fixed_state_mmhpsd_result_",no_events,".Rdata",sep=""))
    
    time_segment <- 0.05
    time_horizon <- tail(test.mmhp.fix$tau,1)
    time.segment <- c(0:round(time_horizon/time_segment))*time_segment
    draw.start<-1
    draw.end<-length(test.mmhp.fix$tau)
    draw.true.x<-which(test.mmhp.fix$x>test.mmhp.fix$tau[draw.start]&test.mmhp.fix$x<test.mmhp.fix$tau[draw.end])
    draw.true.time<-c(test.mmhp.fix$tau[draw.start],test.mmhp.fix$x[draw.true.x],test.mmhp.fix$tau[draw.end])
    draw.true.state<-3-c(test.mmhp.fix$z[head(draw.true.x,1)-1],test.mmhp.fix$z[draw.true.x],test.mmhp.fix$z[tail(draw.true.x,1)])
    
    change <- 9-l*2.5
    
    for(i in 1:50){
      state_mmhpsd <- 2-(state_estimation_matrix[i,]<0.5)
      time_mmhpsd <- c(0:(round(time_horizon/time_segment)))*time_segment
      lines(state_mmhpsd+change~time_mmhpsd, type="l", 
            col=rgb(169,169,169,maxColorValue=255), lwd=0.35)
    }
    
    lower.ci <- apply(state_estimation_matrix<0.5,2,mean)+apply(state_estimation_matrix<0.5,2,sd)
    lower.ci <- ifelse(lower.ci>1,rep(1,length(lower.ci)),lower.ci)
    upper.ci <- apply(state_estimation_matrix<0.5,2,mean)-apply(state_estimation_matrix<0.5,2,sd)
    upper.ci <- ifelse(upper.ci<0,rep(0,length(upper.ci)),upper.ci)
    
    lines(2-apply(state_estimation_matrix<0.5,2,mean)+change~time.segment,type="l",col="blue",lwd=4)
    polygon(c(time.segment,rev(time.segment)),
            c(2-lower.ci+change,rev(2-upper.ci+change)),
            col=rgb(0,0,225,alpha=50,maxColorValue=255),border=FALSE)
    lines(draw.true.state+change~draw.true.time,type="s",col="grey27",lwd=4.5)
    segments(x0=-0.5,y0=7-(l-1)*2.5,y1=9-(l-1)*2.5)
    segments(x0=c(22.5,16,14.3,11)[5-l]+1.25,y0=7-(l-1)*2.5,y1=9-(l-1)*2.5)
    segments(x0=-0.5,x1=c(22.5,16,14.3,11)[5-l]+1.25,y0=7-(l-1)*2.5)
    segments(x0=-0.5,x1=c(22.5,16,14.3,11)[5-l]+1.25,y0=9-(l-1)*2.5)
  }
  
  legend(15,10,c("truth","MMHPSD mean","confidence band","individual process"),
         lty=1, lwd=c(5,5,20,5), x.intersp=0.3,
         col=c(rgb(0,0,0),rgb(0,0,225,maxColorValue=255),
               rgb(0,0,225,alpha=100,maxColorValue=255),rgb(169,169,169,maxColorValue=255)),
         bg=rgb(0,0,0,alpha=0), bty="n",cex=2.75,seg.len=0.85)
  text(c(22.5,15,12.5,9.5),c(0.5,3,5.5,8),paste("M=",rev(no_vec_events),sep=""),cex=2.5)
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

plotIntAbsErrorExample <- function(data.path, no.events, plot.idx=37,
                                   my.mar=c(4,4.3,2.8,0.7)){
  ## mmhpsd
  load(paste(data.path,"fixed_state_mmhpsd_result_",no.events,".Rdata",sep=""))
  time_segment <- 0.05
  time_horizon <- tail(test.mmhp.fix$tau,1)
  time.segment <- c(0:round(time_horizon/time_segment))*time_segment
  draw.start<-1
  draw.end<-length(test.mmhp.fix$tau)
  draw.true.x<-which(test.mmhp.fix$x>test.mmhp.fix$tau[draw.start]&test.mmhp.fix$x<test.mmhp.fix$tau[draw.end])
  draw.true.time<-c(test.mmhp.fix$tau[draw.start],test.mmhp.fix$x[draw.true.x],test.mmhp.fix$tau[draw.end])
  draw.true.state<-3-c(test.mmhp.fix$z[head(draw.true.x,1)-1],test.mmhp.fix$z[draw.true.x],test.mmhp.fix$z[tail(draw.true.x,1)])
  
  ## mmhp
  load(paste(data.path,"fixed_state_stan_result_",no.events,".Rdata",sep=""))
  
  fixed.state.array <- array(NA,dim=c(50,5,7))
  fixed.state.est.latent <- matrix(rep(list(),50*3),nrow=50,ncol=3)
  fixed.state.latent.fun <- matrix(rep(list(),50*2),nrow=50,ncol=2)
  par_name <- c("lambda0","alpha","beta","q1","q2")
  for(i in 1:50){
    for(j in 1:5){
      fixed.state.array[i,j,1] <- mean(unlist(fixed.state.stan.result[i][[1]][par_name[j]])[1001:2000])
      fixed.state.array[i,j,2] <- median(unlist(fixed.state.stan.result[i][[1]][par_name[j]])[1001:2000])
    }
    #estimate latent state
    test.mmhp <- fixed.state.processes[i][[1]]
    object_hat <- list(lambda0=fixed.state.array[i,1,1], alpha=fixed.state.array[i,2,1], beta=fixed.state.array[i,3,1], q1=fixed.state.array[i,4,1], q2=fixed.state.array[i,5,1])
    test.latent <- stepfun(test.mmhp.fix$x[-1],test.mmhp.fix$z)
    time.segment <- seq(test.mmhp$tau[2],tail(test.mmhp$tau,1),length.out=2000)
    time.segment.delta <- time.segment[2]-time.segment[1]
    
    #2 estimation methods
    fixed.state.est.latent[i,1][[1]] <- apply(fixed.state.stan.result[i][[1]]$zt_v[1001:2000,],2,mean)
    fixed.state.est.latent[i,2][[1]] <- ifelse(apply(fixed.state.stan.result[i][[1]]$zt_v[1001:2000,],2,function(x) sum(x==1)) >500,rep(1,1000),rep(2,1000))
    
    #preperation for step function
    fixed.state.est.latent[i,3][[1]] <- modifiedLatentTrajectory(params=object_hat, interevent = test.mmhp$tau[-1]-test.mmhp$tau[-length(test.mmhp$tau)], zt = fixed.state.est.latent[i,2][[1]])
    
    #interpolated function for each method
    fixed.state.latent.fun[i,1][[1]] <- approxfun(x=test.mmhp$tau[-1], y=fixed.state.est.latent[i,1][[1]],rule=2)
    fixed.state.latent.fun[i,2][[1]] <- stepfun(fixed.state.est.latent[i,3][[1]]$x.hat[-1],fixed.state.est.latent[i,3][[1]]$z.hat)
  }
  #prepare for plot
  time.segment <- seq(0,tail(test.mmhp.fix$tau,1),length.out = 5000)
  draw.start<-1
  draw.end<-length(test.mmhp.fix$tau)
  draw.true.x<-which(test.mmhp.fix$x>test.mmhp.fix$tau[draw.start]&test.mmhp.fix$x<test.mmhp.fix$tau[draw.end])
  draw.true.time<-c(test.mmhp.fix$tau[draw.start],test.mmhp.fix$x[draw.true.x],test.mmhp.fix$tau[draw.end])
  draw.true.state<-3-c(test.mmhp.fix$z[head(draw.true.x,1)-1],test.mmhp.fix$z[draw.true.x],test.mmhp.fix$z[tail(draw.true.x,1)])
  
  fun.result <- array(NA,dim=c(50,5000,2))
  for(i in c(plot.idx)){
    for(j in 1:2){
      fun.result[i,,j] <- fixed.state.latent.fun[i,j][[1]](time.segment)
    }
  }
  
  par(mfrow=c(1,1),tcl=0.2,mgp=c(2.4,0,0),
      mar=my.mar, oma=c(0,0,0,0))
  
  plot(0.5,0,xlim=c(-0.2,tail(test.mmhp.fix$tau,1)),ylim=c(0.9,2.1), type="n",
       xlab="Time",ylab="",xaxt="n",yaxt="n",cex.lab=1.8,bty="n")
  title("(a)",line=0.5,cex.main=2.2)
  axis(2,at=c(1,2),labels=FALSE,cex.lab=1.8,las=2,lwd=0,lwd.ticks=0)
  axis(1,at=seq(0,40,5),labels=FALSE,cex.lab=1.8,las=2,lwd=1.5,lwd.ticks=1)
  text(-1.15, y=c(1,2), labels=paste(rep(c("state","state"),4),rep(c(0,1),4)," "), cex=1.8, srt=0, xpd=TRUE)
  text(seq(0,40,5), y=0.78, labels=seq(0,40,10), cex=1.8, srt=0, xpd=TRUE)
  segments(x0=0,x1=tail(test.mmhp.fix$tau,1),y0=1,col="lightgrey",lwd=2)
  segments(x0=0,x1=tail(test.mmhp.fix$tau,1),y0=2,col="lightgrey",lwd=2)
  
  state_mmhpsd <- state_estimation_matrix[plot.idx,] #2-(state_estimation<0.5)
  time_mmhpsd <- c(0:(round(time_horizon/time_segment)))*time_segment
  lines(state_mmhpsd+1~time_mmhpsd, type="l", col="blue", lwd=2.5)
  
  state_mmhp <- fun.result[plot.idx,,2] #2-(state_estimation<0.5)
  time_mmhp <- seq(0,tail(test.mmhp.fix$tau,1),length.out = 5000)
  lines(3-state_mmhp~time_mmhp, type="l",col="red", lwd=2.5)
  
  polygon(c(time_mmhpsd,rev(time_mmhpsd)),
          c(state_mmhpsd+1,3-rev(test.latent(time_mmhpsd))),
          col=rgb(0,0,225,alpha=50,maxColorValue=255),border=FALSE)
  polygon(c(time_mmhp,rev(time_mmhp)),
          c(3-state_mmhp,3-rev(test.latent(time_mmhp))),
          col=rgb(225,0,0,alpha=50,maxColorValue=255),border=FALSE) 
  lines(draw.true.state~draw.true.time,type="s",col="black",lwd=4)
}

plotIAEBeeswarm <- function(data.path, model.vec, my.mar=c(4,4.3,2.8,0.7)){
  fun.error <- array(0,dim=c(4,2,50),
                     dimnames = list(c("50events","100events","200events","500events"),
                                     c("mmhp","mmhpsd"),
                                     c(1:50)))
  for(l in c(1:4)){
    no_events <- no_vec_events[l]
    for(j in c(1:2)){
      load(paste(data.path,model.vec[j],no_events,".Rdata",sep=""))
      
      for(i in c(1:50)){
        test.mmhp <- fixed.state.processes[i][[1]]
        test.latent <- stepfun(test.mmhp.fix$x[-1],test.mmhp.fix$z)
        par_name <- c("lambda0","alpha","beta","q1","q2")
        if(j==2){
          time_horizon <- tail(test.mmhp.fix$tau,1)
          time_segment <- 0.05
          time_vec <- time_segment*c(0:round(time_horizon/time_segment))
          fixed.state.est.latent <- state_estimation_matrix[i,]
          fun.error[l,j,i] <- sum(abs(state_estimation_matrix[i,]-test.latent(time_vec)))*time_segment
        }else{
          fixed.state.array <- numeric(5)
          time.segment <- seq(test.mmhp$tau[1],tail(test.mmhp$tau,1),0.05)#length.out=2000)
          time.segment.delta <- time.segment[2]-time.segment[1]
          for(k in 1:5){
            fixed.state.array[k] <- mean(unlist(fixed.state.stan.result[i][[1]][par_name[k]])[1001:2000])
          }
          object_hat <- list(lambda0=fixed.state.array[1], 
                             alpha=fixed.state.array[2], 
                             beta=fixed.state.array[3], 
                             q1=fixed.state.array[4], 
                             q2=fixed.state.array[5])
          fixed.state.est.latent <- ifelse(apply(fixed.state.stan.result[i][[1]]$zt_v[1001:2000,],2,function(x) sum(x==1)) >500,rep(1,1000),rep(2,1000))
          fixed.state.est.latent.new <- modifiedLatentTrajectory(params=object_hat, 
                                                                 interevent = test.mmhp$tau[-1]-test.mmhp$tau[-length(test.mmhp$tau)], 
                                                                 zt = fixed.state.est.latent)
          fixed.state.latent.fun <- stepfun(fixed.state.est.latent.new$x.hat[-1],
                                            fixed.state.est.latent.new$z.hat)
          fun.error[l,j,i] <- sum(abs(fixed.state.latent.fun(time.segment)-test.latent(time.segment)))*time.segment.delta
        }
      }
    }
  }
  
  fun.error.df <- data.frame(method=rep(c(rep("MMHP",50),rep("MMHPSD",50)),4),
                             error=c(fun.error[4,1,],fun.error[4,2,],
                                     fun.error[3,1,],fun.error[3,2,],
                                     fun.error[2,1,],fun.error[2,2,],
                                     fun.error[1,1,],fun.error[1,2,]),
                             number_events=factor(rep(c(500,200,100,50),each=100)),levels=c(500,200,100,50))
  
  f <- ordered(fun.error.df$number_events, levels = c(500,200,100,50))
  par(mfrow=c(1,1),tcl=0.2,mgp=c(2.7,0,0),
      mar=my.mar, mai = c(0.85,0.9,0.4,3.2))
  
  beeswarm(error ~ f, data = fun.error.df, vertical = FALSE,
           method = "center", pch = 16, cex=1, spacing = 1.1,
           pwcol = c(rgb(255,0,0,alpha=190,maxColorValue = 255),
                     rgb(0,0,255,alpha=190,maxColorValue = 255))[ifelse(method=="MMHP",rep(1,100),rep(2,100))],
           xlab = "Integrated absolute error", ylab = "Number of events", corral = "gutter",cex.lab=1.8,
           xaxt="n",yaxt="n",bty = "n")
  text(-2.52, y=c(1:4), labels=c("500", "200", "100", "50"), cex=1.8, srt=0, xpd=TRUE)
  text(seq(0,20,5), y=0.12, labels=seq(0,20,5), cex=1.8, srt=0, xpd=TRUE)
  title("(b)",line=-0.7,cex.main=2.2)
  axis(1,at=seq(0,20,5),labels=FALSE,cex.lab=1.8,las=2,lwd=1.5,lwd.ticks=1.5)
  
  legend(x = 24.5, y = 4, xpd = TRUE, bty = "n",
         legend = c("true state",
                    "MMHP state", "MMHPSD state",
                    "MMHP error area", "MMHPSD error area", 
                    "MMHP error value", "MMHPSD error value"),
         col = c("black","red", "blue",
                 rgb(225,0,0,alpha=50,maxColorValue=255),
                 rgb(0,0,225,alpha=50,maxColorValue=255), 
                 "red", "blue"), 
         lty = c(1,1,1,1,1,NA,NA), lwd = c(4,4,4,20,20,NA,NA),
         pch = c(NA,NA,NA,NA,NA,16,16), pt.cex = c(NA,NA,NA,NA,NA,1.5,1.5), cex=1.5) 
}
