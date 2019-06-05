data{
  int<lower=1> N_til;//number of pairs => nrow(unique_pairs_df)
  int<lower=1> no_observations;//total number of observation windows -> max(return_df$observe.id),
  int<lower=1> max_Nm; //maximum of number of events for each pair each window => max(unlist(lapply(return_df$event.times,length))))
  int<lower=0,upper=max_Nm> Nm[N_til,no_observations]; //number of events for each pair => count_matrix
  vector[max_Nm+1] time_matrix[N_til,no_observations]; // include termination time in the last entry
  real<lower=0> max_interevent[N_til];
}
parameters{
  vector<lower=0>[N_til] lambda0; //baseline rate for each pair
  vector<lower=0>[N_til] c; //baseline rate for each pair
  vector<lower=0,upper=1>[N_til] w1; //CTMC transition rate
  vector<lower=0,upper=1>[N_til] w2; //CTMC transition rate
}
transformed parameters{
  vector<lower=0>[N_til] q1;
  vector<lower=0>[N_til] q2;
  q1 = (lambda0).*w1;
  q2 = (lambda0).*w2;
}
model{
  real integ; // Placeholder variable for calculating integrals
  row_vector[2] forward[max_Nm]; // Forward variables from forward-backward algorithm
  row_vector[2] forward_termination; 
  row_vector[2] probs_1[max_Nm+1]; // Probability vector for transition to state 1 (active state)
  row_vector[2] probs_2[max_Nm+1]; // Probability vector for transition to state 2 (inactive state)
  vector[max_Nm+1] interevent;
  real temp_lambda0;
  real temp_c;
  real temp_q1;
  real temp_q2;

  //priors
  c ~ lognormal(0,1);
  w1 ~ beta(0.5,0.5);
  w2 ~ beta(0.5,0.5);
  
  for(i in 1:N_til){ //for each pair
    lambda0[i] ~ normal(1/max_interevent[i],0.5);
    temp_lambda0 = lambda0[i];
    temp_c = c[i];
    temp_q1 = q1[i];
    temp_q2 = q2[i];
    for(j in 1:no_observations){ // for each observation period
      interevent = time_matrix[i,j];
      if(Nm[i,j]==0){ // there is no event occured in this period
        integ = interevent[1]*temp_lambda0;
        probs_1[1][1] = log(temp_q2/(temp_q1+temp_q2)+temp_q1/(temp_q1+temp_q2)*exp(-(temp_q1+temp_q2)*interevent[1])); //1->1
        probs_2[1][2] = log(temp_q1/(temp_q1+temp_q2)+temp_q2/(temp_q1+temp_q2)*exp(-(temp_q1+temp_q2)*interevent[1])); //2->2
        probs_1[1][2] = log1m_exp(probs_2[1][2]); //2->1
        probs_2[1][1] = log1m_exp(probs_1[1][1]); //1->2

        forward_termination[1] = log_sum_exp(probs_1[1])-integ*(1+temp_c);
        forward_termination[2] = log_sum_exp(probs_2[1])-integ;
        target += log_sum_exp(forward_termination);
      }else{ // there is event occured
        // ---- prepare for forward algorithm
        // --- log probability of Markov transition logP_ij(t)
        for(n in 1:(Nm[i,j]+1)){
           probs_1[n][1] = log(temp_q2/(temp_q1+temp_q2)+temp_q1/(temp_q1+temp_q2)*exp(-(temp_q1+temp_q2)*interevent[n])); //1->1
           probs_2[n][2] = log(temp_q1/(temp_q1+temp_q2)+temp_q2/(temp_q1+temp_q2)*exp(-(temp_q1+temp_q2)*interevent[n])); //2->2
           probs_1[n][2] = log1m_exp(probs_2[n][2]); //2->1
           probs_2[n][1] = log1m_exp(probs_1[n][1]); //1->2
        }
    
        //consider n = 1
        integ = interevent[1]*temp_lambda0;
        forward[1][1] = log_sum_exp(probs_1[1]) + log(temp_lambda0*(1+temp_c)) - integ*(1+temp_c); 
        forward[1][2] = log_sum_exp(probs_2[1]) + log(temp_lambda0) - integ; 
        
        if(Nm[i,j]>1){
          for(n in 2:Nm[i,j]){
            integ = interevent[n]*temp_lambda0;
            forward[n][1] = log_sum_exp(forward[n-1] + probs_1[n]) + log(temp_lambda0*(1+temp_c))- integ*(1+temp_c);
            forward[n][2] = log_sum_exp(forward[n-1] + probs_2[n]) + log(temp_lambda0) - integ;
          }
        }
        
        integ = interevent[Nm[i,j]+1]*temp_lambda0;
        forward_termination[1] = log_sum_exp(forward[Nm[i,j]] + probs_1[Nm[i,j]+1]) - integ*(1+temp_c);
        forward_termination[2] = log_sum_exp(forward[Nm[i,j]] + probs_2[Nm[i,j]+1]) - integ;
        
        target += log_sum_exp(forward_termination);
      }
    }
  }
}
