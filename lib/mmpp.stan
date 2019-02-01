data{
	int<lower=0> N;
	real<lower=0> interevent[N];
}
parameters{
	real<lower=0> lambda0; //baseline rate
	real<lower=0> c; //increment ratio over baseline rate
  real<lower=0,upper=1> w1;
  real<lower=0,upper=1> w2;
  //real<lower=0,upper=1> w_q2;
  //real<lower=0> q1; //CTMC transition rate
	//real<lower=0> q2;
}
transformed parameters{
	row_vector[2] probs_1[N]; // Probability vector for transition to state 1 (active state)
	row_vector[2] probs_2[N]; // Probability vector for transition to state 2 (inactive state)
  real<lower=0> q1; //CTMC transition rate
	real<lower=0> q2;
	
	q1 = lambda0*w1;
	q2 = lambda0*w2;
	
	for(n in 1:N){
	  probs_1[n][1] = -q1*interevent[n];
    probs_2[n][2] = -q2*interevent[n];
    probs_1[n][2] = log1m_exp(probs_2[n][2]);
    probs_2[n][1] = log1m_exp(probs_1[n][1]);
	}
}
model{
	real integ; // Placeholder variable for calculating integrals
	row_vector[2] forward[N]; // Forward variables from forward-backward algorithm

	//priors
	lambda0 ~ normal(1,0.01);
	c ~ lognormal(0,3);
	w1 ~ beta(0.5,0.5);
  w2 ~ beta(0.5,0.5);
  //q1 ~ lognormal(0,1);
  //q2 ~ lognormal(0,1);

	//consider n = 1
	integ = interevent[1]*lambda0;
	forward[1][1] = log((1+c)*lambda0) - (1+c)*integ;
	forward[1][2] = log(lambda0) - integ; //calculate forward variables, uniform initial distribution for latent state

	for(n in 2:N){
	  integ = interevent[n]*lambda0;
	  forward[n][1] = log_sum_exp(forward[n-1] + probs_1[n]) + log((1+c)*lambda0) - (1+c)*integ;
	  forward[n][2] = log_sum_exp(forward[n-1] + probs_2[n]) + log(lambda0) - integ;
	}

	target += log_sum_exp(forward[N]);
}
generated quantities{
  //-------------Parameters for local decoding
  real<lower=0,upper=1> pzt[N]; // probability of event at active state
  int<lower=1,upper=2> zt[N]; // local decoding
  row_vector[2] forward[N]; // forward and backward variables
  row_vector[2] backward[N];
  real integral;
  int l;
  row_vector[2] intensities;
  row_vector[2] integrals;

  //-------------Parameters for global decoding (Viterbi algorithm)
  int<lower=1,upper=2> zt_v[N]; // global decoding
  real m[N,2]; // forward and backward variables
  int b[N,2];
  real logp1;
  real logp2;

  //-------------local decoding
  //calculate forward variables
  integral = interevent[1]*lambda0;
  forward[1][1] = log((1+c)*lambda0) - (1+c)*integral;
	forward[1][2] = log(lambda0) - integral;
	m[1,1] = log((1+c)*lambda0) - (1+c)*integral;
  m[1,2] = log(lambda0) - integral;

	for(n in 2:N){
	  integral = interevent[n]*lambda0;
	  //----local
	  forward[n][1] = log_sum_exp(forward[n-1] + probs_1[n]) + log((1+c)*lambda0) - (1+c)*integral;
	  forward[n][2] = log_sum_exp(forward[n-1] + probs_2[n]) + log(lambda0) - integral;
	  //----global
	  logp1 = m[n-1,1] + probs_1[n][1];
	  logp2 = m[n-1,2] + probs_1[n][2];
	  if(logp1>logp2){
	    m[n,1] = logp1 + log((1+c)*lambda0) - (1+c)*integral;
	    b[n,1] = 1;
	  }else{
        m[n,1] = logp2 + log((1+c)*lambda0) - (1+c)*integral;
        b[n,1] = 2;
	  }
	  logp1 = m[n-1,1] + probs_2[n][1];
	  logp2 = m[n-1,2] + probs_2[n][2];
	  if(logp1>logp2){
	    m[n,2] = logp1 + log(lambda0) - integral;
	    b[n,2] = 1;
	  }else{
        m[n,2] = logp2 + log(lambda0) - integral;
        b[n,2] = 2;
	  }
	}

	//calculate backward variables for local decoding
	backward[N][1] = 0;
	backward[N][2] = 0;
	intensities[1] = (1+c)*lambda0;
	intensities[2] = lambda0;
	for(n in 2:N){
	  l = N-n+1;
	  integrals[2] = lambda0*interevent[l+1];
	  integrals[1] = (1+c)*integrals[2];
	  backward[l][1] = log_sum_exp(backward[l+1] + probs_1[l+1] + log(intensities) - integrals);
	  backward[l][2] = log_sum_exp(backward[l+1] + probs_2[l+1] + log(intensities) - integrals);
	}
  //infer the latent state
  for(n in 1:N){
    pzt[n] = inv_logit(forward[n][1] + backward[n][1] - forward[n][2] - backward[n][2]);
    if(forward[n][2] + backward[n][2] > forward[n][1] + backward[n][1]){
      zt[n] = 2;
    }else{
      zt[n] = 1;
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

