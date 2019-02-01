data{
	int<lower=0> N; // number of events
	real<lower=0> interevent[N]; // interevent time
	real<lower=0> event[N]; // interevent time
}
parameters{
	real<lower=0> lambda1; //baseline rate
	real<lower=0,upper=1> w_lambda; 
	//real<lower=0,upper=1> w_q1; //CTMC transition rate
	//real<lower=0,upper=1> w_q2;
	real<lower=0> q2;
  real<lower=0> q1;
	real<lower=0> alpha;
	real<lower=0> beta;
}
transformed parameters{
  real<lower=0> lambda0; //baseline rate
  //real<lower=0> q2;
  //real<lower=0> q1;
  row_vector[2] probs_1[N]; // Probability vector for transition to state 1 (active state)
  row_vector[2] probs_2[N]; // Probability vector for transition to state 2 (inactive state)
  lambda0 = lambda1*w_lambda;
  //q1 = lambda0*w_q1;
  //q2 = q1*w_q2;
   
   for(n in 1:N){
     probs_1[n][1] = log(q2/(q1+q2)+q1/(q1+q2)*exp(-(q1+q2)*interevent[n])); //1->1
     probs_2[n][2] = log(q1/(q1+q2)+q2/(q1+q2)*exp(-(q1+q2)*interevent[n])); //2->2
     probs_1[n][2] = log1m_exp(probs_2[n][2]); //2->1
     probs_2[n][1] = log1m_exp(probs_1[n][1]); //1->2
  }
}
model{
	row_vector[2] forward[N]; // Forward variables from forward-backward algorithm
	row_vector[2] integ[2]; // integral vector
  real r[N]; // record variable for Hawkes process
  real A; //record variable for integration
  
	//priors
	lambda1 ~ lognormal(0,1);
	alpha ~ normal(0,5);
	beta ~ lognormal(0,0.5);
	q1 ~ lognormal(-1,1);
	q2 ~ lognormal(-1,1);

	//consider n = 1
	forward[1][1] = log(lambda1) - interevent[1]*lambda1;
	forward[1][2] = log(lambda0) - interevent[1]*lambda0; //calculate forward variables, uniform initial distribution for latent state
	r[1] = 0; 
	
	for(n in 2:N){
	  r[n] = exp(-beta*interevent[n])*(r[n-1]+1);
	  forward[n][1] = log_sum_exp(forward[n-1] + probs_1[n]) + log(lambda1+alpha*exp(-beta*interevent[n])*(r[n-1]+1)) - interevent[n]*lambda1 + alpha/beta*(r[n]-r[n-1]-1);
	  forward[n][2] = log_sum_exp(forward[n-1] + probs_2[n]) + log(lambda0) - interevent[n]*lambda0;
	}
	target += log_sum_exp(forward[N]);
}

generated quantities{
	real r[N];
  
  //-------------Parameters for global decoding (Viterbi algorithm)
  int<lower=1,upper=2> zt_v[N]; // global decoding
  real m[N,2]; // forward and backward variables
  int b[N,2];
  real logp1;
  real logp2;
  
  m[1,1] = log(lambda1) - interevent[1]*lambda1;
  m[1,2] = log(lambda0) - interevent[1]*lambda0; 
	r[1] = 0;
	
	for(n in 2:N){
	  r[n] = exp(-beta*interevent[n])*(r[n-1]+1);
	  
	  //----global
	  logp1 = m[n-1,1] + probs_1[n][1];
	  logp2 = m[n-1,2] + probs_1[n][2];
	  if(logp1>logp2){ //1->1
	    m[n,1] = logp1 + log(lambda1+alpha*exp(-beta*interevent[n])*(r[n-1]+1)) - interevent[n]*lambda1;
	    b[n,1] = 1;
	  }else{ //2->1
        m[n,1] = logp2 + log(lambda1+alpha*exp(-beta*interevent[n])*(r[n-1]+1)) - interevent[n]*lambda1;
        b[n,1] = 2;
	  }
	  logp1 = m[n-1,1] + probs_2[n][1];
	  logp2 = m[n-1,2] + probs_2[n][2];
	  if(logp1>logp2){ //1->2
	    m[n,2] = logp1 + log(lambda0) - interevent[n]*lambda0;
	    b[n,2] = 1;
	  }else{ //2->2
        m[n,2] = logp2 + log(lambda0) - interevent[n]*lambda0;
        b[n,2] = 2;
	  }
	}
  
  //calulate zt_v for global decoding
  if(m[N,1]>m[N,2]){
	  zt_v[N] = 1;
	}else{
	  zt_v[N] = 2;
	}
	
	for(n in 1:(N-1)){
	  zt_v[N-n] = b[N-n+1,zt_v[N-n+1]];
	}
}
