sameTimeHorizonPlot <- function(xl,yl,clean_data,c="black",ltype=1,lwidth=5){
  total_days <- tail(clean_data$day,1)
  previous_no_observe_time <- 0
  previous_x <- 0
  previous_y <- 0
  
  for(day in c(1:total_days)){
    day_idx <- which(clean_data$day == day)
    day_observation_time <- unique(clean_data$new_observation[day_idx])
    if(exists("day_observation_idx")){
      previous_no_observe_time <- previous_no_observe_time + clean_data$time_stamp[min(day_observation_idx)]
    }
    for(observation_time in c(1:length(day_observation_time))){
      day_observation_idx <- which(clean_data$new_observation == day_observation_time[observation_time])
      observation_start_time <- clean_data$day_hour[min(day_observation_idx)]
      observation_end_time <- clean_data$day_hour[max(day_observation_idx)]
      observation_duration <- observation_end_time - observation_start_time
      
      if(max(xl)<observation_start_time|min(xl)>observation_end_time){
        next
      }
      
      start_idx <- min(which(xl>=observation_start_time))
      end_idx <- max(which(xl<=observation_end_time))
      
      if(observation_time>1){
        previous_no_observe_time <- previous_no_observe_time + clean_data$time_stamp[min(day_observation_idx)]-
          clean_data$time_stamp[min(day_observation_idx)-1]
      }
      
      lines(xl[start_idx:end_idx]+previous_no_observe_time,
            yl[start_idx:end_idx],
            cex=0.6, col=c,lwd=lwidth,lty=ltype)
      
      if(day_observation_time[observation_time]==1){
        previous_x <- xl[end_idx]+previous_no_observe_time
        previous_y <- yl[end_idx]
      }else{
        segments(x0=previous_x,
                 y0=previous_y,
                 x1=xl[start_idx]+previous_no_observe_time,
                 y1=yl[start_idx],
                 cex=0.6, col=c,lwd=0.6)
        previous_x <- xl[end_idx]+previous_no_observe_time
        previous_y <- yl[end_idx]
      }
    }
    previous_no_observe_time <- previous_no_observe_time + 24 - clean_data$time_stamp[max(day_observation_idx)]
  }  
}
