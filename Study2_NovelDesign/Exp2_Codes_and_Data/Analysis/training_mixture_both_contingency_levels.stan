// Stan model - Study 2. RT: 2 Gaussians mixture model, with 2 contingency conditions
data {
  // Input: data shape 
  int N_subjects; // n=20 or n = 58
  int N_trials; // maximal number of Go trials
  int N_trials_valid[N_subjects]; // actual number of trials of each participant
  // Input: Dependent (RT-effective) and independent (Run and contingency) data
  real RT[N_trials, N_subjects]; // dependent variable: RT (RT - Cue)
  real Run[N_trials, N_subjects]; // training run [0,1] indicating [1st, 20th] 
  real Contingency[N_trials, N_subjects]; // 50% (0) or 100% (1) contingency condition
  real theta_b0; // b0 = -3.1, mixture proportion defined at 1st run (set to 0.1%)
  real Cue; // fixed at 850ms
}
parameters {
  // Group-level parameters of the two Gaussians       
  real<lower=0, upper=Cue-100> mu1_fix; // early Gaussian Mean (upper limit at 100ms before cue-onset)
  real<lower=Cue> mu2_fix; // late Gaussian Mean (lower limit at cue-onset)
  real<lower=0> sigma_e[2]; // SD of the two Gaussians
  real<lower=0> mu1_sd[1]; // SD of individualized early Gaussian mean parameter 
  // Group-level (fixed-effect) parameters: theta-slope individualized learning parameter 
  real theta_fix[2]; // Means of individualized theta-slope parameters (for each contingency condition)
  real<lower=0> theta_sd[2]; // SDs of individualized theta-slope parameter
  // Participant-level (random-effect) parameter: individualized learning parameter and anticipatory mean
  real <upper=Cue+100> mu1_rand[N_subjects]; // individualized parameter for early RT Gaussian's mean
  real theta_b1_rand[N_subjects]; // individualized parameter to determine mixture proportion (100% contingency)
  real theta_b2_rand[N_subjects]; // individualized parameter to determine mixture proportion (50% contingency)
  
}
model {
  //priors
  mu1_fix ~ normal(500,500);
  mu2_fix ~ normal(1000,500);
  sigma_e ~ normal(0,1000);
  mu1_sd ~ normal(0,150);
  theta_fix ~ normal(0,0.7);
  theta_sd ~ normal(0.7,0.7);
  //likelihood
  for (s in 1:N_subjects){ // go over participants
  // participant-level [random-effect] from group-level [fix-effect] parameters
  mu1_rand[s] ~ normal(mu1_fix,mu1_sd); // individualized early RT mean
  theta_b1_rand[s] ~ normal(theta_fix[1],theta_sd[1]); // 100% contingency
  theta_b2_rand[s] ~ normal(theta_fix[2],theta_sd[2]); // 50% contingency
  for (t in 1:N_trials_valid[s]){ // go over participantג€™s trials
  /* define a linear function of training run with individualized theta-slope use the
  inverse of the Normal distribution CDF to go from [-inf,inf] to [0,1] range */
  // maximize likelihood of parameters' fit to RT data
  target += log_mix(Phi_approx(theta_b0 + (Contingency[t,s]*theta_b1_rand[s] + (1-Contingency[t,s])*theta_b2_rand[s]) * Run[t,s]),
  normal_lpdf(RT[t,s] | mu1_rand[s] , sigma_e[1]), // anticipatory RT
  normal_lpdf(RT[t,s] | mu2_fix, sigma_e[2])); // cue-dependent RT
  }
  }
}
