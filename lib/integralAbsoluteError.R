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