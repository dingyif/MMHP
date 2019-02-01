metricsStateSeparation <- function(start, end, type, rank){
  if(type == "total"){
    if(start == 0){
      m <- adjm[end,rank,rank]
    }else{
      m <- adjm[end,rank,rank] - adjm[start,rank,rank]
    }
  }else if(type == "utility"){
    if(start == 0){
      m <- utility_state_cumcount[end,rank,rank]
    }else{
      m <- utility_state_cumcount[end,rank,rank] - utility_state_cumcount[start,rank,rank]
    }
  }else{
    if(start == 0){
      m <- adjm[end,rank,rank] - utility_state_cumcount[end,rank,rank]
    }else{
      m <- adjm[end,rank,rank] - adjm[start,rank,rank] - (utility_state_cumcount[end,rank,rank] - utility_state_cumcount[start,rank,rank])
    }
  }
  
  return(list(metric_dci=dci(m),
              metric_ttri=ttri(m)[[2]],
              metric_i=iAndSI(m)[1]))
}
