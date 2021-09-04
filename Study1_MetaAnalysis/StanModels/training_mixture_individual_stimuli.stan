
// Stan model - Study 2. RT: 2 Gaussians mixture model, with stimulus-level learning parameters
data {
// Input: data shape 
  int N_subjects; // Number of participants
	int N_trials; // maximal number of Go trials
	int N_stim; // number of unique Go stimuli-participants
  int N_trials_valid[N_subjects]; // actual number of trials of each participant
  // Input: Dependent (RT-effective) and independent (Run and contingency) data
  real RT[N_trials, N_subjects]; // dependent variable: effective RT (RT - Cue)
 	real Run[N_trials, N_subjects]; // training run [0,1] indicating [1st, 20th] 
  int Stim[N_trials, N_subjects]; // stimulus indicator
  real theta_b0; // b0 = -3.1, mixture proportion defined at 1st run (set to 0.1%)
}
parameters {
  // Group-level parameters of the two Gaussians       
  real <upper=100> mu1_fix; // early Gaussian Mean (upper limit: 100ms after cue)
  real <lower=200> mu2_fix; // late Gaussian Mean (lower limit: 200ms after cue)
	real <lower=0> sigma_e[2]; // SD of the two Gaussians
	// Group-level (fixed-effect) parameters: theta-slope individualized learning parameter 
	real theta_slope_fix; // individualized theta-slope parameter Means
	real <lower=0> theta_stim_sd; // SD of stimulus-level theta-slope parameter
	// Stimulus-level (random-effect) parameters: individualized learning parameters
	real theta_slope_stim [N_stim]; // stimuli individualized parameter 
}
model {
  //priors
  mu1_fix ~ normal(500,500);
  mu2_fix ~ normal(1000,500);
  sigma_e ~ normal(0,1000);
  theta_slope_fix ~ normal(0,0.7);
  // theta_sub_sd ~ cauchy(0,1);
	theta_stim_sd ~ cauchy(0,1);
  // likelihood
  // Stimulus-level [random-effect] from group-level [fix-effect] parameters
  theta_slope_stim ~ normal(theta_slope_fix, theta_stim_sd);
  for (i in 1:N_subjects){ // go over participants
  	for (t in 1:N_trials_valid[i]){ // go over participant’s trials
	  // identify stimulus index
    int stim_i;
    stim_i = Stim[t,i];
  	/* define a linear function of training run with individualized theta-slope use the
  	inverse of the Normal distribution CDF to go from [-inf,inf] to [0,1] range */
  	// maximize likelihood of parameters' fit to RT data
  	target += log_mix(Phi_approx(theta_b0 + theta_slope_stim[stim_i] * Run[t,i]), // mixture proportion
  	normal_lpdf(RT[t,i] | mu1_fix, sigma_e[1]), // anticipatory RT
  	normal_lpdf(RT[t,i] | mu2_fix, sigma_e[2])); // cue-dependent RT
  	}
  	}
}

