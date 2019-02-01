data{
  int<lower=1> max_Nm; //maximum of number of events for each pair
  int<lower=1> N_til;//number of pairs who have more than "cut_off" interactions
  int<lower=0,upper=max_Nm> Nm[N_til]; //number of events for each pair
  vector[max_Nm] time_matrix[N_til];
  real max_interevent[N_til];
}
parameters{
  real mu_alpha;
  real mu_beta;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  vector<lower=0>[N_til] lambda0; //baseline rate for each pair
  vector<lower=0,upper=1>[N_til] w_q2; //CTMC transition rate
  real<lower=0> alpha[N_til];
  real<lower=0> beta[N_til];
  vector<lower=0,upper=1>[N_til] w;
}
transformed parameters{
  vector<lower=0>[N_til] q2;
  vector<lower=0>[N_til] q1;
  q2 = (lambda0).*w_q2;
  q1 = (q2).*w;
}
model{
  real integ; // Placeholder variable for calculating integrals
  row_vector[2] forward[max_Nm]; // Forward variables from forward-backward algorithm
  row_vector[2] probs_1[max_Nm]; // Probability vector for transition to state 1 (active state)
  row_vector[2] probs_2[max_Nm]; // Probability vector for transition to state 2 (inactive state)
  real r[max_Nm]; // record variable for Hawkes process
  vector[max_Nm] interevent;

  //priors
  alpha ~ lognormal(mu_alpha,sigma_alpha);
  beta ~ lognormal(mu_beta,sigma_beta);
  w_q2 ~ beta(0.5,0.5);
  w ~ beta(0.5,0.1);
  
  for(i in 1:N_til){
    lambda0[i] ~ normal(1/max_interevent[i],0.1);
    interevent = time_matrix[i];

    //consider n = 1
    integ = interevent[1]*lambda0[i];
    forward[1][1] = log(lambda0[i]) - integ;
    forward[1][2] = log(lambda0[i]) - integ; //calculate forward variables, uniform initial distribution for latent state
    r[1] = 0; 
    probs_1[1][1] = log(q2[i]/(q1[i]+q2[i])+q1[i]/(q1[i]+q2[i])*exp(-(q1[i]+q2[i])*interevent[1])); //1->1
    probs_2[1][2] = log(q1[i]/(q1[i]+q2[i])+q2[i]/(q1[i]+q2[i])*exp(-(q1[i]+q2[i])*interevent[1])); //2->2
    probs_1[1][2] = log1m_exp(probs_2[1][2]); //2->1
    probs_2[1][1] = log1m_exp(probs_1[1][1]); //1->2
      
    for(n in 2:Nm[i]){
      probs_1[n][1] = log(q2[i]/(q1[i]+q2[i])+q1[i]/(q1[i]+q2[i])*exp(-(q1[i]+q2[i])*interevent[n])); //1->1
      probs_2[n][2] = log(q1[i]/(q1[i]+q2[i])+q2[i]/(q1[i]+q2[i])*exp(-(q1[i]+q2[i])*interevent[n])); //2->2
      probs_1[n][2] = log1m_exp(probs_2[n][2]); //2->1
      probs_2[n][1] = log1m_exp(probs_1[n][1]); //1->2
      integ = interevent[n]*lambda0[i];
      r[n] = exp(-beta[i]*interevent[n])*(r[n-1]+1);
      forward[n][1] = log_sum_exp(forward[n-1] + probs_1[n]) + log(lambda0[i]+alpha[i]*exp(-beta[i]*interevent[n])*(r[n-1]+1)) - integ + alpha[i]/beta[i]*(r[n]-r[n-1]-1);
      forward[n][2] = log_sum_exp(forward[n-1] + probs_2[n]) + log(lambda0[i]) - integ;
    }
    target += log_sum_exp(forward[Nm[i]]);
  }
}
generated quantities{
  //-------------Parameters
  row_vector[2] probs_1[max_Nm]; // Probability vector for transition to state 1 (active state)
  row_vector[2] probs_2[max_Nm]; // Probability vector for transition to state 2 (inactive state)
  real r[max_Nm];
  vector[max_Nm] interevent;
  real integral;
  int l;
  real temp;
  
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
    r[1] = 0; 
    m[1,1] = log(lambda0[i]) - integral;
    m[1,2] = log(lambda0[i]) - integral; 
    
    probs_1[1][1] = log(q2[i]/(q1[i]+q2[i])+q1[i]/(q1[i]+q2[i])*exp(-(q1[i]+q2[i])*interevent[1])); //1->1
    probs_2[1][2] = log(q1[i]/(q1[i]+q2[i])+q2[i]/(q1[i]+q2[i])*exp(-(q1[i]+q2[i])*interevent[1])); //2->2
    probs_1[1][2] = log1m_exp(probs_2[1][2]); //2->1
    probs_2[1][1] = log1m_exp(probs_1[1][1]); //1->2
    
    for(n in 2:Nm[i]){
      probs_1[n][1] = log(q2[i]/(q1[i]+q2[i])+q1[i]/(q1[i]+q2[i])*exp(-(q1[i]+q2[i])*interevent[n])); //1->1
      probs_2[n][2] = log(q1[i]/(q1[i]+q2[i])+q2[i]/(q1[i]+q2[i])*exp(-(q1[i]+q2[i])*interevent[n])); //2->2
      probs_1[n][2] = log1m_exp(probs_2[n][2]); //2->1
      probs_2[n][1] = log1m_exp(probs_1[n][1]); //1->2
      integral = interevent[n]*lambda0[i];

      r[n] = exp(-beta[i]*interevent[n])*(r[n-1]+1);
      temp = log(lambda0[i]+alpha[i]*exp(-beta[i]*interevent[n])*(r[n-1]+1)) + alpha[i]/beta[i]*(r[n]-r[n-1]-1);
      
      //----global
      logp1 = m[n-1,1] + probs_1[n][1];
      logp2 = m[n-1,2] + probs_1[n][2];
      if(logp1>logp2){
        m[n,1] = logp1 + temp - integral;
        b[n,1] = 1;
      }else{
          m[n,1] = logp2 + temp - integral;
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
