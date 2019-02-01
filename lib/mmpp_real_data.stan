data{
  int<lower=1> max_Nm; //maximum of number of events for each pair
  int<lower=1> N_til;//number of pairs who have more than two interactions
  int<lower=0,upper=max_Nm> Nm[N_til]; //number of events for each pair
  vector[max_Nm] time_matrix[N_til];
  real max_interevent[N_til];
}
parameters{
  //real mu; //mean for baseline rate
  //real mu_c;
  //real<lower=0> sigma;
  //real<lower=0> sigma_c;
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
  row_vector[2] probs_1[max_Nm]; // Probability vector for transition to state 1 (active state)
  row_vector[2] probs_2[max_Nm]; // Probability vector for transition to state 2 (inactive state)
  vector[max_Nm] interevent;

  //priors
  //mu ~ normal(1,5);
  //sigma ~ inv_gamma(2,4);
  c ~ lognormal(0,1);
  w1 ~ beta(0.5,0.5);
  w2 ~ beta(0.5,0.5);

  //consider n = 1
  for(i in 1:N_til){
    lambda0[i] ~ normal(1/max_interevent[i],0.1);
    interevent = time_matrix[i];
    integ = interevent[1]*lambda0[i];
    forward[1][1] = log((1+c[i])*lambda0[i]) - (1+c[i])*integ;
    forward[1][2] = log(lambda0[i]) - integ; //calculate forward variables, uniform initial distribution for latent state

    probs_1[1][1] = -q1[i]*interevent[1];
    probs_2[1][2] = -q2[i]*interevent[1];
    probs_1[1][2] = log1m_exp(probs_2[1][2]);
    probs_2[1][1] = log1m_exp(probs_1[1][1]);
      
    for(n in 2:Nm[i]){
      probs_1[n][1] = -q1[i]*interevent[n];
      probs_2[n][2] = -q2[i]*interevent[n];
      probs_1[n][2] = log1m_exp(probs_2[n][2]);
      probs_2[n][1] = log1m_exp(probs_1[n][1]);
      integ = interevent[n]*lambda0[i];
      forward[n][1] = log_sum_exp(forward[n-1] + probs_1[n]) + log(lambda0[i]*(1+c[i])) - integ*(1+c[i]);
      forward[n][2] = log_sum_exp(forward[n-1] + probs_2[n]) + log(lambda0[i]) - integ;
    }
    target += log_sum_exp(forward[Nm[i]]);
  }
}
generated quantities{
  //-------------Parameters for local decoding
  vector<lower=0,upper=1>[max_Nm] pzt[N_til]; // probability of event at active state for all pairs
  vector<lower=1,upper=2>[max_Nm] zt[N_til]; // local decoding
  row_vector[2] forward[max_Nm]; // forward and backward variables
  row_vector[2] backward[max_Nm];
  row_vector[2] probs_1[max_Nm]; // Probability vector for transition to state 1 (active state)
  row_vector[2] probs_2[max_Nm]; // Probability vector for transition to state 2 (inactive state)
  vector[max_Nm] interevent;
  real integral;
  int l;
  row_vector[2] intensities;
  row_vector[2] integrals;
  
  //-------------Parameters for global decoding (Viterbi algorithm)
  int<lower=1,upper=2> zt_v[max_Nm,N_til]; // global decoding
  real m[max_Nm,2]; // forward and backward variables
  int b[max_Nm,2];
  real logp1;
  real logp2;

  //------------consider all the processes
  for(i in 1:N_til){
    interevent = time_matrix[i];
    //-------------local decoding
    //-----calculate forward variables
    integral = interevent[1]*lambda0[i];
    forward[1][1] = log((1+c[i])*lambda0[i]) - integral*(1+c[i]);
    forward[1][2] = log(lambda0[i]) - integral; //calculate forward variables, uniform initial distribution for latent state
    m[1,1] = log(lambda0[i]*(1+c[i])) - integral*(1+c[i]);
    m[1,2] = log(lambda0[i]) - integral; 
    
    probs_1[1][1] = -q1[i]*interevent[1];
    probs_2[1][2] = -q2[i]*interevent[1];
    probs_1[1][2] = log1m_exp(probs_2[1][2]);
    probs_2[1][1] = log1m_exp(probs_1[1][1]);
    
    for(n in 2:Nm[i]){
      probs_1[n][1] = -q1[i]*interevent[n];
      probs_2[n][2] = -q2[i]*interevent[n];
      probs_1[n][2] = log1m_exp(probs_2[n][2]);
      probs_2[n][1] = log1m_exp(probs_1[n][1]);
      integral = interevent[n]*lambda0[i];

      //----local
      forward[n][1] = log_sum_exp(forward[n-1] + probs_1[n]) + log(lambda0[i]*(1+c[i])) - integral*(1+c[i]); 
      forward[n][2] = log_sum_exp(forward[n-1] + probs_2[n]) + log(lambda0[i]) - integral;

      //----global
      logp1 = m[n-1,1] + probs_1[n][1];
      logp2 = m[n-1,2] + probs_1[n][2];
      if(logp1>logp2){
        m[n,1] = logp1 + log(lambda0[i]*(1+c[i])) - integral*(1+c[i]);
        b[n,1] = 1;
      }else{
          m[n,1] = logp2 + log(lambda0[i]*(1+c[i])) - integral*(1+c[i]);
          b[n,1] = 2;
      }
      logp1 = m[n-1,1] + probs_2[n][1];
      logp2 = m[n-1,2] + probs_2[n][2];
      if(logp1>logp2){
        m[n,2] = logp1 + log(lambda0[i]) - integral;
        b[n,2] = 1;
      }else{
        m[n,2] = logp2 + log(lambda0[i]) - integral;
        b[n,2] = 2;
      }
    }
    
    //-----calculate backward variables
    backward[Nm[i]][1] = 0;
    backward[Nm[i]][2] = 0;
    intensities[1] = lambda0[i]*(1+c[i]);
    intensities[2] = lambda0[i];
    for(n in 2:Nm[i]){
      l = Nm[i]-n+1;
      integrals[2] = lambda0[i]*interevent[l+1];
      integrals[1] = integrals[2]*(1+c[i]);
      backward[l][1] = log_sum_exp(backward[l+1] + probs_1[l+1] + log(intensities) - integrals);
      backward[l][2] = log_sum_exp(backward[l+1] + probs_2[l+1] + log(intensities) - integrals);
    }
    //----infer the latent state for local decoding
    for(n in 1:Nm[i]){
      pzt[i][n] = inv_logit(forward[n][1] + backward[n][1] - forward[n][2] - backward[n][2]);
      if(forward[n][2] + backward[n][2] > forward[n][1] + backward[n][1]){
        zt[i][n] = 2;
      }else{
        zt[i][n] = 1;
      }
    }
     if(Nm[i]!=max_Nm){
       for(n in (Nm[i]+1):max_Nm){
         pzt[i][n] = 0;
         zt[i][n] = 2;
       }
     }
    
    //----calulate zt_v for global decoding
    if(m[Nm[i],1]>m[Nm[i],2]){
      zt_v[Nm[i],i] = 1;
    }else{
      zt_v[Nm[i],i] = 2;
    }
  
    for(n in 1:(Nm[i]-1)){
      zt_v[Nm[i]-n,i] = b[Nm[i]-n+1,zt_v[Nm[i]-n+1,i]];
    }
    
    if(Nm[i]!=max_Nm){
      for(n in (Nm[i]+1):max_Nm){
        zt_v[n,i] = 2;
      }
    }
  }
}
