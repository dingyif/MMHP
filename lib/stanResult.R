stanResult <- function(cohort_number,...){
  
  # check for additional function arguments
  if(length(list(...))){
    Lst <- list(...)
    if( !is.null(Lst$path) ){
      path <- Lst$path
    }else{
      path <- paste(getwd(),"/stan_result/",sep='')
    }
  }
  
  # clean the data 
  clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
  
  # if the model has been fitted, load the Rdata from the input path,
  # otherwise fit the model and save the result in the path
  if(file.exists(paste(path, cohort_names[current_cohort], "/mmhp_stan_result_", cohort_names[current_cohort], ".RData", sep=''))){
    load(paste(path, cohort_names[current_cohort], "/mmhp_stan_result_", cohort_names[current_cohort], ".RData",sep=''))
  }else{
    fit_mmhp <- stan("../lib/mmhp_real_data.stan",
                     data=list(max_Nm=max(clean_data$N_count),
                               N_til=sum(clean_data$N_count>=cut_off),
                               Nm=as.vector(clean_data$N_count[clean_data$N_count>=cut_off]),
                               time_matrix=clean_data$time_matrix,
                               max_interevent=clean_data$max_interevent),
                     iter=1000, chains=4)
    sim_mmhp_pair <- rstan::extract(fit_mmhp)
    dir.create(paste(path,cohort_names[current_cohort],sep=''), recursive = TRUE)
    save(sim_mmhp_pair, file = paste(path,cohort_names[current_cohort],"/mmhp_stan_result_",cohort_names[current_cohort],".RData",sep=''))
  }
  
  # reorganize the stan output into matrix (or matrix of list for latent process)
  ztv_mmhp_pair <- t(ifelse(apply(sim_mmhp_pair$zt_v[1001:2000,,],c(2,3),function(x) sum(x==1)) >500,rep(1,1000),rep(2,1000)))
  pztv_mmhp_pair <- t(apply(sim_mmhp_pair$zt_v[1001:2000,,],c(2,3),mean))
  lambda0_pair <- apply(sim_mmhp_pair$lambda0[1001:2000,],2,mean)
  lambda1_pair <- apply(sim_mmhp_pair$lambda0[1001:2000,],2,mean)
  alpha_pair <- apply(sim_mmhp_pair$alpha[1001:2000,],2,mean)
  beta_pair <- apply(sim_mmhp_pair$beta[1001:2000,],2,mean)
  q1_pair <- apply(sim_mmhp_pair$q1[1001:2000,],2,mean)
  q2_pair <- apply(sim_mmhp_pair$q2[1001:2000,],2,mean)
  
  lambda0_matrix <- matrix(0,ncol=mice_number,nrow=mice_number)
  lambda1_matrix <- matrix(0,ncol=mice_number,nrow=mice_number)
  alpha_matrix <- matrix(0,ncol=mice_number,nrow=mice_number)
  beta_matrix <- matrix(0,ncol=mice_number,nrow=mice_number)
  q1_matrix <- matrix(0,ncol=mice_number,nrow=mice_number)
  q2_matrix <- matrix(0,ncol=mice_number,nrow=mice_number)
  ztv_list_pair <- matrix(rep(list(),mice_number^2),mice_number,mice_number)
  pztv_list_pair <- matrix(rep(list(),mice_number^2),mice_number,mice_number)
  ztv_list_sample_pair <- array(rep(list(),mice_number^2*1000),dim=c(mice_number,mice_number,1000))
  
  for(i in 1:sum(clean_data$N_count>=cut_off)){
    lambda0_matrix[clean_data$I_fit[i],clean_data$J_fit[i]] <- lambda0_pair[i]
    lambda1_matrix[clean_data$I_fit[i],clean_data$J_fit[i]] <- lambda1_pair[i]
    alpha_matrix[clean_data$I_fit[i],clean_data$J_fit[i]] <- alpha_pair[i]
    beta_matrix[clean_data$I_fit[i],clean_data$J_fit[i]] <- beta_pair[i]
    q1_matrix[clean_data$I_fit[i],clean_data$J_fit[i]] <- q1_pair[i]
    q2_matrix[clean_data$I_fit[i],clean_data$J_fit[i]] <- q2_pair[i]
    ztv_list_pair[clean_data$I_fit[i],clean_data$J_fit[i]][[1]] <- 
      ztv_mmhp_pair[i,1:(clean_data$N_count[clean_data$I_fit[i],clean_data$J_fit[i]])]
    pztv_list_pair[clean_data$I_fit[i],clean_data$J_fit[i]][[1]] <- 
      pztv_mmhp_pair[i,1:(clean_data$N_count[clean_data$I_fit[i],clean_data$J_fit[i]])]
    for(k in 1:1000){
      ztv_list_sample_pair[clean_data$I_fit[i],clean_data$J_fit[i],k][[1]] <- 
        sim_mmhp_pair$zt_v[k+1000,1:(clean_data$N_count[clean_data$I_fit[i],clean_data$J_fit[i]]),i]
    }
  }
  
  return(list(lambda0_matrix = lambda0_matrix,
              lambda1_matrix = lambda1_matrix,
              alpha_matrix = alpha_matrix,
              beta_matrix = beta_matrix,
              q1_matrix = q1_matrix,
              q2_matrix = q2_matrix,
              ztv_list_pair = ztv_list_pair,
              pztv_list_pair = pztv_list_pair,
              ztv_list_sample_pair = ztv_list_sample_pair,
              mu_alpha = mean(sim_mmhp_pair$mu_alpha[1001:2000]),
              sigma_alpha = mean(sim_mmhp_pair$sigma_alpha[1001:2000]),
              mu_beta = mean(sim_mmhp_pair$mu_beta[1001:2000]),
              sigma_beta = mean(sim_mmhp_pair$sigma_beta[1001:2000])))
}