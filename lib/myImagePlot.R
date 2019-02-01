myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  smallScale <- FALSE
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
    if( !is.null(Lst$smallScale) ){
      smallScale <- Lst$smallScale
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(7,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  if( is.null(Lst$colorPalette) ){
    ColorRamp <- colorRampPalette(colorschemes$LightBluetoDarkBlue.7)(length(xLabels)*length(yLabels))
  }else{
    ColorRamp <- c("#FFFFFF",colorRampPalette(brewer.pal(8,colorPalette))(length(xLabels)*length(yLabels)))
  }
  
  #ColorRamp <- colorRampPalette(brewer.pal(8,colorPalette))(length(xLabels)*length(yLabels))
  #cm.colors(256) #rgb( seq(0,1,length=256),  # Red
  #      seq(0,1,length=256),  # Green
  #     seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, ifelse(smallScale,5,max), length=length(ColorRamp)) #min max
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(tcl=-0.2,mgp=c(0,1,0),mar = c(3,3.2,2,1.3))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=1.8)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=1.8)
  
  # Color Scale
  par(tcl=-0.2,mgp=c(1.8,0,0),mar = c(3,1.2,2,1.5))
  if(!smallScale){
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n", yaxt="n")
    axis(LEFT <- 2, at=seq(min,max,length.out = 5), labels=c(0,10,20,30,40),
         cex.axis=1.5)
  }else{
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp[1:(length(ColorLevels))],
          xlab="",ylab="",
          xaxt="n",cex.axis=1.5)
  }
  layout(1)
}