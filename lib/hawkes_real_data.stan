data{
  int<lower=1> N_til;//number of pairs => nrow(unique_pairs_df)
  int<lower=1> no_observations;//total number of observation windows -> max(return_df$observe.id),
  int<lower=1> max_Nm; //maximum of number of events for each pair each window => max(unlist(lapply(return_df$event.times,length))))
  int<lower=0,upper=max_Nm> Nm[N_til,no_observations]; //number of events for each pair => count_matrix
  vector[max_Nm+1] interevent_time_matrix[N_til,no_observations]; // include termination time in the last entry
  vector[max_Nm] event_matrix[N_til,no_observations]; // event times in each observation window
  real<lower=0> finishing_time[no_observations]; //for each pair, each observation window, what is the finishing time
}
parameters{
  vector<lower=0>[N_til] lambda1; //
  vector<lower=0>[N_til] alpha; //
  vector<lower=0>[N_til] beta; //
}
model{
  real r[max_Nm]; // record variable for Hawkes process
  vector[max_Nm+1] interevent;
  vector[max_Nm] event;
  
  lambda1 ~ lognormal(0,0.5);
  beta ~ normal(0,10);
  //tilde_beta ~ normal(0,20);

  for(i in 1:N_til){
    for(j in 1:no_observations){
      if(Nm[i,j]==0){ 
        target += -lambda1[i]*finishing_time[j];
      }else{
        interevent = interevent_time_matrix[i,j];
        target += -lambda1[i]*finishing_time[j] - alpha[i]/beta[i]*sum((1-exp(-beta[i]*(finishing_time[j]-segment(event_matrix[i,j],1,Nm[i,j])))));
        r[1] = 0; 
        target += log(lambda1[i]+alpha[i]*r[1]);
        if(Nm[i,j]>1){
          for(n in 2:Nm[i,j]){
            r[n] = exp(-beta[i]*interevent[n])*(r[n-1]+1); 
            target += log(lambda1[i]+alpha[i]*r[n]);
          }
        }
      }
    }
  }
}
