handelExceptionMmhpsd <- function(event_time, max_iter = 6){
  parameter_list <- c(list(lamb0 = c(0.001, 1), 
                             nu0 = c(0.1, 0.6),
                             eta0 = c(0.04, 0.9),
                             Q0 = matrix(c(-0.08,8,0.08,-8),nrow=2,ncol=2),
                             pai0 = c(0.5,0.5)),
                      list(lamb0 = c(0.01, 30), 
                             nu0 = c(0.01, 0.03),
                             eta0 = c(0.4, 9),
                             Q0 = matrix(c(-0.5,0.4,0.5,-0.4),nrow=2,ncol=2),
                             pai0 = c(0.5,0.5)),
                      list(lamb0 = c(1, 30), 
                           nu0 = c(0.01, 0.03),
                           eta0 = c(0.4, 9),
                           Q0 = matrix(c(-0.5,0.4,0.5,-0.4),nrow=2,ncol=2),
                           pai0 = c(0.5,0.5)),
                      list(lamb0 = c(1, 30), 
                           nu0 = c(0.1, 0.03),
                           eta0 = c(0.4, 9),
                           Q0 = matrix(c(-0.05,0.4,0.05,-0.4),nrow=2,ncol=2),
                           pai0 = c(0.5,0.5)),
                      list(lamb0 = c(1, 3), 
                           nu0 = c(0.01, 0.03),
                           eta0 = c(0.4, 0.9),
                           Q0 = matrix(c(-0.05,0.4,0.05,-0.4),nrow=2,ncol=2),
                           pai0 = c(0.5,0.5)),
                      list(lamb0 = c(0.1, 3), 
                           nu0 = c(0.01, 0.03),
                           eta0 = c(0.4, 0.9),
                           Q0 = matrix(c(-0.05,0.4,0.05,-0.4),nrow=2,ncol=2),
                           pai0 = c(0.5,0.5)))
  return_result <- try(log("error"))
  i <- 1
  while(class(return_result) == "try-error" & i<max_iter){
    return_result <- try(mmhp(ti = event_time,
                              lamb0 = parameter_list[[(i-1)*5+1]], 
                              nu0 = parameter_list[[(i-1)*5+2]],
                              eta0 = parameter_list[[(i-1)*5+3]],
                              Q0 = parameter_list[[(i-1)*5+4]],
                              pai0 = parameter_list[[(i-1)*5+5]],
                              nlmp = TRUE, fortran = TRUE))
    i <- i+1
  }
  return(return_result)
}