// Stan model mixture gaus model - training with random theta
      data {
      int N_subjects;
      int N_trials;
      int N_trials_valid[N_subjects];
      real Cue;
      real theta_b0;
      real RT[N_trials, N_subjects];
      real Run[N_trials, N_subjects];
      real Contingency[N_trials, N_subjects];
      }
      
      parameters {
      // Fixed effect: mu1(Anticipatory), mu2(Cue dependent), theta1(100% contingency), theta2(50% contingency), sigma_e1, sigma_e2
      real<lower=0, upper=Cue+100> mu1_fix; // Anticipatory mean
      real<lower=Cue> mu2_fix; // Cue dependent mean
      real theta_fix[2]; // theta0 + slope (for run)
      real<lower=0> sigma_e[2]; // error terms
      
      // SD for random effect: mu1, theta1, theta2
      real<lower=0> mu1_sd[1]; // SD of the random effects around the fixed mean 
      real<lower=0> theta_sd[2]; // SD of the random effects around the fixed thetas
      
      // Random effect parameter: mu1, theta1, theta2
      real <upper=Cue+100> mu1_rand[N_subjects]; // random intercept for early onest Gaussian's mean
      // real theta_b0_rand[N_subjects]; // random intercept for theta - fix at -2.5
      real theta_b1_rand[N_subjects]; // random intercept for theta
      real theta_b2_rand[N_subjects]; // random intercept for theta

      }
      
      
      model {
      //priors
      mu1_fix ~ normal(500,500);
      mu2_fix ~ normal(1000,500);
      theta_fix ~ normal(0,0.7);
      sigma_e ~ normal(0,1000);
      theta_sd ~ normal(0.7,0.7);
      mu1_sd ~ normal(0,150);
      
      
      //likelihood
      for (s in 1:N_subjects){
      mu1_rand[s] ~ normal(mu1_fix,mu1_sd);
      // theta_b0_rand[s] ~ normal(theta_fix[1],theta_sd[1]); fix at -2.5
      theta_b1_rand[s] ~ normal(theta_fix[1],theta_sd[1]); // 100% contingency
      theta_b2_rand[s] ~ normal(theta_fix[2],theta_sd[2]); // 50% contingency
      
      for (t in 1:N_trials_valid[s]){
      target += log_mix(Phi_approx(theta_b0 + (Contingency[t,s]*theta_b1_rand[s] + (1-Contingency[t,s])*theta_b2_rand[s]) * Run[t,s]),
      normal_lpdf(RT[t,s] | mu1_rand[s] , sigma_e[1]),
      normal_lpdf(RT[t,s] | mu2_fix, sigma_e[2]));
      }
      }
      }
