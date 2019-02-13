################################################################################################
#This R file is source code to generate intensity figure for the latex file
#There are two functions in this source code: drawUniMMHPIntensityPaper and drawHPIntensityPaper
#Function 'drawHPIntensityPaper' is used to assist main function 'drawUniMMHPIntensityPaper'
################################################################################################

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
