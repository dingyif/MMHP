# ----- Define a function for plotting a matrix ----- #
matrixPlotParameter <- function(x, ...){
  min <- min(x)
  max <- max(x)
  #abs<-6.9
  #min<--abs
  #min<--0.78
  #max<-20
  #min<-0
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- colorRampPalette(rev(brewer.pal(8,"RdBu")))(length(xLabels)*length(yLabels))
     #cm.colors(256) #rgb( seq(0,1,length=256),  # Red
               #      seq(0,1,length=256),  # Green
               #     seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp)) #min max
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  #yLabels <- yLabels[reverse]
  #x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}
# ----- END plot function ----- #

# ----- Define a function for plotting a matrix ----- #
matrix.poly <- function(x, y, z=mat, ROW, COL){
  dist.left <- (x[ROW]-x[ROW-1])/2
  dist.right <- (x[ROW+1]-x[ROW])/2
  if(ROW==1) dist.left <- dist.right
  if(ROW==dim(z)[1]) dist.right <- dist.left
  
  dist.down <- (y[COL]-y[COL-1])/2
  dist.up <- (y[COL+1]-y[COL])/2
  if(COL==1) dist.down <- dist.up
  if(COL==dim(z)[2]) dist.up <- dist.down
  
  xs <- c(x[ROW]-dist.left, x[ROW]-dist.left, x[ROW]+dist.right, x[ROW]+dist.right)
  ys <- c(y[COL]-dist.down, y[COL]+dist.up, y[COL]+dist.up, y[COL]-dist.down)
  poly <- data.frame(x=xs, y=ys)
  
  return(poly)
}
myMatrixPlot <- function(x,emp.COL,emp.ROW, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- colnames(x)
  xLabels <- rownames(x)
  colorPalette <- "RdBu"
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
    if( !is.null(Lst$colorPalette) ){
      colorPalette <- Lst$colorPalette
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- colorRampPalette(brewer.pal(8,colorPalette))(length(xLabels)*length(yLabels))
  #cm.colors(256) #rgb( seq(0,1,length=256),  # Red
  #      seq(0,1,length=256),  # Green
  #     seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp)) #min max
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  #Emphasize on some entries
  #random_N <- sample(1:(length(xLabels)*length(yLabels)),1)
  poly <- matrix.poly(1:length(xLabels),1:length(yLabels),t(x),ROW=emp.COL,COL=length(yLabels)-emp.ROW+1)
  polygon(poly,col="green",border=1)
  box()
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp, xlab="",ylab="", xaxt="n")
  
  layout(1)
}

# ----- Define a function for plotting multiple matrices ----- #
myMultiMatrixPlot <- function(X,xLabels,yLabels,min,max,colorPalette,no_matrix,if.emp=TRUE,emp.ROW,emp.COL,emp.color,title){
  par( oma = c( 0, 0, 3, 0 ) )
  layout(matrix(data=c(1:(2*no_matrix)), nrow=1, ncol=2*no_matrix), widths=rep(c(4,0.6),no_matrix), heights=rep(1,2*no_matrix))
  
  for(i in 1:no_matrix){
    x <- X[[i]]
    min_temp <- min[[i]]
    max_temp <- max[[i]]
    yLabels_temp <- yLabels[[i]]
    xLabels_temp <- xLabels[[i]]
    ColorRamp <- colorRampPalette(brewer.pal(8,colorPalette[[i]]))(length(xLabels_temp)*length(yLabels_temp))
    ColorLevels <- seq(min_temp, max_temp, length=length(ColorRamp)) 
    
    # Reverse Y axis
    reverse <- nrow(x) : 1
    yLabels_temp <- yLabels_temp[reverse]
    x <- x[reverse,]
    
    # Data Map
    par(mar = c(3,5,2.5,0.6))
    image(1:length(xLabels_temp), 1:length(yLabels_temp), t(x), col=ColorRamp, xlab="",
          ylab="", axes=FALSE, zlim=c(min_temp,max_temp))
    axis(BELOW<-1, at=1:length(xLabels_temp), labels=xLabels_temp, cex.axis=0.7)
    axis(LEFT <-2, at=1:length(yLabels_temp), labels=yLabels_temp, las= HORIZONTAL<-1, cex.axis=0.7)
    
    if(if.emp){
      #Emphasize on some entries
      poly <- matrix.poly(1:length(xLabels_temp),1:length(yLabels_temp),t(x),ROW=emp.COL[[i]],COL=length(yLabels_temp)-emp.ROW[[i]]+1)
      polygon(poly,col=ifelse(emp.color[[i]]=="transparent",rgb(1,1,1,0),emp.color[[i]]),border=1)
    }
  
    
    # Color Scale
    par(mar = c(3,0.6,2.5,2))
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
  }
  mtext(title,outer=TRUE)
  layout(1)
}
# ----- END plot function ----- #

decay_adjm <- function(m, new_S_temp){
  t <- new_day_hour[m]
  result_matrix <- matrix(0,ncol=N,nrow=N)
  for(i in 1:N){
    for(j in c(1:N)[-i]){
      index_temp <- new_day_hour[unlist(indicator.each.pair[i,j])]<=t
      if(sum(index_temp)>=1){
        result_matrix[i,j] <- sum(exp(-(t-new_day_hour[unlist(indicator.each.pair[i,j])][index_temp])))
      }
    }
  }
  return(result_matrix)
}