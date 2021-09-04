// Stan model mixture gaus model - training with random theta
      data {
      int N_subjects;
      int N_trials;
      int N_trials_valid[N_subjects];
      real RT[N_trials, N_subjects];
      real Run[N_trials, N_subjects];
      real theta_b0;
      }
      
      parameters {
      real <upper=0> mu_fix_1; // fixed intercept - 2 distributions' means
      real <lower=0> mu_fix_2; // fixed intercept - 2 distributions' means
      real theta_fix; // fixed intercept for theta + slope
      // real theta_b1_fix[1]; // fixed slope for run
      real<lower=0> theta_sd; // SD of the random effects around the fixed effect
      real<lower=0> sigma_e[2];
      real theta_b1[N_subjects]; // random intercept for theta
      }
      
      model {
      //priors
      mu_fix_1 ~ normal(-500,500);
      mu_fix_2 ~ normal(500,500);
      theta_fix ~ normal(0,1);
      theta_sd ~ normal(0,0.7);
      sigma_e ~ normal(0,1000);
      
      
      //likelihood
      for (s in 1:N_subjects){
      real theta;
      theta_b1[s] ~ normal(theta_fix,theta_sd);
      
      for (t in 1:N_trials_valid[s]){
      theta = Phi_approx(theta_b0 + (theta_b1[s] * Run[t,s]));
      target += log_mix(theta,
      normal_lpdf(RT[t,s] | mu_fix_1 , sigma_e[1]),
      normal_lpdf(RT[t,s] | mu_fix_2 , sigma_e[2]));
      }
      }
      }
      
