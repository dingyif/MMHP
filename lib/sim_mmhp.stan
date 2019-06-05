data{
	int<lower=0> N; // number of events
	real<lower=0> interevent[N]; // interevent time
	real<lower=0> event[N]; // event time
	real<lower=0,upper=1> delta_1; // P(initial state = 1)
}
transformed data{
  row_vector[2] log_delta;
  log_delta[1] = log(delta_1);
  log_delta[2] = log(1-delta_1);
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
  row_vector[2] int_1[N]; // Integration of lambda when state transit to 1 (active state)
  row_vector[2] int_2[N]; // Integration of lambda when state transit to 2 (inactive state)
  real R[N];
  real K0;
  real K1;
  real K2;
  real K3;
  real K4;
  real K5;
  
  lambda0 = lambda1*w_lambda;
   
   // --- log probability of Markov transition logP_ij(t)
   for(n in 1:N){
     probs_1[n][1] = log(q2/(q1+q2)+q1/(q1+q2)*exp(-(q1+q2)*interevent[n])); //1->1
     probs_2[n][2] = log(q1/(q1+q2)+q2/(q1+q2)*exp(-(q1+q2)*interevent[n])); //2->2
     probs_1[n][2] = log1m_exp(probs_2[n][2]); //2->1
     probs_2[n][1] = log1m_exp(probs_1[n][1]); //1->2
  }
  
  // --- R for Hawkes
  R[1] = 0;
  for(n in 2:N){
	  R[n] = exp(-beta*interevent[n])*(R[n-1]+1);
  }
  
  // Integration of lambda
  for(n in 1:N){
    K0 = exp(-(q1+q2)*interevent[n]);
    K1 = (1-exp(-(q1+q2)*interevent[n]))/(q1+q2);
    K2 = (1-exp(-(q1+q2)*interevent[n]))/(q1+q2);
    K3 = R[n]*(exp(beta*interevent[n])-1)/beta;
    K4 = R[n]*(1-exp(-(beta+q1+q2)*interevent[n]))*exp(beta*interevent[n])/(beta+q1+q2);
    K5 = R[n]*(1-exp(-(q1+q2-beta)*interevent[n]))/(q1+q2-beta);
    int_1[n][1] = ((q2^2*lambda1+q2*q1*lambda0)*interevent[n] +
                   (q1^2*lambda1+q2*q1*lambda0)*K0*interevent[n] +
                   (lambda1-lambda0)*q2*q1*K1 + (lambda1-lambda0)*q2*q1*K2 +
                   alpha*K3*(q2^2+q1^2*K0) +
                   alpha*q1*q2*K4 + alpha*q1*q2*K5)/(q1+q2)^2/exp(probs_1[n][1]); //1->1
    int_1[n][2] = ((q2^2*lambda1+lambda0*q1*q2)*interevent[n] -
                   (lambda1*q1*q2+lambda0*q2^2)*K0*interevent[n] +
                   (lambda0-lambda1)*q2^2*K1 + (lambda1-lambda0)*q1*q2*K2 +
                   alpha*q2*K3*(q2-q1*K0) -
                   alpha*q2^2*K4 + alpha*q1*q2*K5)/(q1+q2)^2/exp(probs_1[n][2]); //2->1
    int_2[n][1] = ((q1*q2*lambda1+q1^2*lambda0)*interevent[n] -
                   (q1^2*lambda1+q1*q2*lambda0)*K0*interevent[n] +
                   (lambda1-lambda0)*q1^2*K1 + q1*q2*(lambda0-lambda1)*K2 +
                   alpha*q1*K3*(q2-q1*K0) +
                   alpha*q1^2*K4 - alpha*q2*q1*K5)/(q1+q2)^2/exp(probs_2[n][1]); //1->2
    int_2[n][2] = ((q1*q2*lambda1+lambda0*q1^2)*interevent[n] +
                   (q1*q2*lambda1+lambda0*q2^2)*K0*interevent[n] +
                   (lambda0-lambda1)*q1*q2*K1 + (lambda0-lambda1)*q1*q2*K2 +
                   alpha*q1*q2*K3*(1+K0) -
                   alpha*q1*q2*K4 - alpha*q1*q2*K5)/(q1+q2)^2/exp(probs_2[n][2]); //2->2
  }
}

model{
	row_vector[2] forward[N]; // Forward variables from forward-backward algorithm
  
	//priors
	lambda1 ~ lognormal(0,1);
	alpha ~ normal(0,5);
	beta ~ lognormal(0,0.5);
	q1 ~ lognormal(-1,1);
	q2 ~ lognormal(-1,1);

	//consider n = 1
	forward[1][1] = log(lambda1) + log_sum_exp(probs_1[1]-int_1[1]+log_delta); 
	forward[1][2] = log(lambda0) + log_sum_exp(probs_2[1]-int_2[1]+log_delta); 
	
	for(n in 2:N){
	  forward[n][1] = log_sum_exp(forward[n-1] + probs_1[n] - int_1[n]) + log(lambda1+alpha*R[n]);
	  forward[n][2] = log_sum_exp(forward[n-1] + probs_2[n] - int_2[n]) + log(lambda0);
	}
	target += log_sum_exp(forward[N]);
}
