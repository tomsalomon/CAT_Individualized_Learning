
library(lme4)
library(rstudioapi)    
library(rstan)
library(reshape2)
library(lm.beta)
library(lmerTest)
library(dplyr)
library(tidyr)
library(gganimate)

rm(list=ls())
options(mc.cores = parallel::detectCores())

# Fix terrible plotting on Windows
trace(grDevices:::png, quote({if (missing(type) && missing(antialias)) {
  type <- "cairo-png"
  antialias <- "subpixel"}}), print = FALSE)

# Get current path
pwd = dirname(rstudioapi::getActiveDocumentContext()$path)
pathSplit=strsplit(pwd, "/")
pathSplit=pathSplit[[1]]
main_path=paste0(pathSplit[1:(length(pathSplit)-1)],"/",collapse="")
analysis_path=paste0(pwd,"/")
path=paste0(main_path,"/Output/")
models_path = paste0(pwd,'/StanModels/')
setwd(pwd)

# Valid participants (n = 59):
subjects = c(101:115,117:118,121:128,131:162,164:165)

# Excluded participants:
# 116, 130 - Transitivity
# 163 - FA (max = 37.59%)
# 119, 120 - Miss (max = 38.75%)
# 129 - Technical issues

# Load data ----
# CAT
filelist=c()
for (s in subjects){
  filelist=c(filelist,Sys.glob(paste(path, "BM_",s,"_CAT_*.txt",sep="")))
}
CAT_data=c()
for (f in filelist){
  CAT_data=rbind(CAT_data,read.table(f,header=T,na.strings=c(999,999000)))
}
# Probe
filelist=c()
for (s in subjects){
  filelist=c(filelist,Sys.glob(paste(path, "BM_",s,"_Probe_*.txt",sep="")))
}
Probe_data=c()
for (f in filelist){
  Probe_data=rbind(Probe_data,read.table(f,header=T,na.strings=c(999,999000)))
}
Probe_data$sub_i = as.numeric(factor(Probe_data$subjectID))

CAT_data$sub_i = as.numeric(factor(CAT_data$Subject))
CAT_data$Contingency2 = factor(CAT_data$Contingency, labels =  c("0% Contingency", "50% Contingency","100% Contingency"))
CAT_data$Run2 = as.factor(CAT_data$Run)
CAT_data$Response = 1 - is.na(CAT_data$RT)
CAT_data$RT_eff = CAT_data$RT - CAT_data$CueOnset
# Dataframe of Go trials with response
CAT_data_GO = subset(CAT_data, Go==1 & Response == 1 & RT>=250)

# RT distribution plot ----
ggplot(data = CAT_data_GO) +
  facet_grid(. ~ Contingency2) +
  geom_density(aes(x = RT, y = ..density.., group=Run2, color = Run)) +
  scale_color_gradient2(low = "firebrick1", mid = "yellow3", high = "green4", midpoint = 10) +
  labs(color = "Training\nRun", x = "Reaction time (ms)") +
  theme_bw() +
  annotate(geom="text", x=400, y=.006, label="Cue Onset", size = 3) +
  geom_vline(xintercept = 850, linetype =2) +
  scale_y_continuous(limits =c(0,.007), expand = c(0,0)) + 
  scale_x_continuous(limits =c(0,1500), expand = c(0,0)) + 
  #theme(legend.title = element_text(size = 3), legend.text = element_text(size = 5))+ 
  theme(legend.key.size=unit(.4,"cm"),panel.spacing = unit(2, "lines")) 

# RT distribution plot ----
ggplot(data = CAT_data_GO) +
  facet_grid(. ~ Contingency2) +
  geom_density(aes(x = RT, y = ..density.., color = Run2)) +
  labs(color = "Training\nRun", x = "Reaction time (ms)") +
  theme_bw() +
  annotate(geom="text", x=400, y=.006, label="Cue Onset", size = 3) +
  geom_vline(xintercept = 850, linetype =2) +
  scale_y_continuous(limits =c(0,.007), expand = c(0,0)) + 
  scale_x_continuous(limits =c(0,1500), expand = c(0,0)) +   
  theme(legend.key.size=unit(.4,"cm"),panel.spacing = unit(2, "lines")) 


# Mean RT linear trend
ggplot(data = CAT_data_GO, aes(x = Run, y = RT, color = Contingency2)) +
  geom_point(position = 'jitter', alpha = 0.3) +
  geom_smooth(method = 'lm')

mean_RT = CAT_data_GO %>% group_by(Run,Contingency, Contingency2) %>% summarise(MeanRT = mean(RT))
ggplot(data = mean_RT, aes(x = Run, y = MeanRT, fill = Contingency2)) +
  geom_bar(stat = 'identity', position = 'dodge') 

ggplot(data = mean_RT, aes(x = Run, y = MeanRT - 850, fill = Contingency2)) +
  geom_bar(stat = 'identity', position = 'dodge') 

# Logistic regression plot
CAT_data_GO = CAT_data_GO %>% mutate(anticipatoty = RT <= (CueOnset + 100))
ggplot(CAT_data_GO, aes(x=Run, y=as.numeric(anticipatoty), color = Contingency2)) + geom_jitter(width = 0.5, height = 0.01, alpha=0.1) + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)

ggplot(CAT_data_GO, aes(x=Run, y=as.numeric(anticipatoty), color = Contingency2)) + 
  facet_wrap(. ~ Subject) + 
  geom_jitter(width = 0.5, height = 0.01, alpha=0.1) + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)
# RT animation
anim = ggplot(data = CAT_data_GO) +
  geom_histogram(aes(x = RT, y = ..density.., fill = Contingency2),
                 color ="black", alpha = 0.4, bins = 50, position = "identity") + 
  geom_density(aes(x = RT, y = ..density.., 
                   group = Contingency2), alpha = .2)+
  labs(color = "Training\nRun", x = "Reaction time (ms)") +
  theme_bw() +
  annotate(geom="text", x=400, y=.006, label="Cue Onset", size = 3) +
  geom_vline(xintercept = 850, linetype =2) +
  scale_y_continuous(limits =c(0,.008), expand = c(0,0)) + 
  scale_x_continuous(limits =c(0,1500), expand = c(0,0)) + 
  theme(legend.title = element_text(size = 3), legend.text = element_text(size = 5))+ 
  theme(legend.key.size=unit(.4,"cm"),panel.spacing = unit(2, "lines")) + 
  transition_states(Run2,
                    transition_length = 2,
                    state_length = 2) + 
  ggtitle('Run number: {closest_state} of 20')
animate(anim, type = "cairo")
anim_save(animation = anim, filename = "Training_RT_animation.gif",type = 'cairo')

# STAN model - both contingency levels: data organization ====
stan_file = "training_mixture_both_contingency_levels.stan"

subs2test = 1:1000 # set to smaller number to examine less participants for speed
dat = subset(CAT_data_GO, sub_i %in% subs2test)
N_trials_per_sub = dat %>% group_by(sub_i) %>% summarise(n()) %>% pull
N_trials_max = max(N_trials_per_sub)
N_runs = max(dat$Run)
N_subjects = nlevels(factor(dat$sub_i))

RT_table = matrix(data = 999000,nrow = N_trials_max,ncol = N_subjects)
Run_table = RT_table
Contingency_table = RT_table
N_trials_valid = c()
for (n in 1:N_subjects){
  # find values
  RT_tmp = dat$RT[dat$sub_i == n]
  Run_tmp = dat$Run[dat$sub_i == n]
  Contingency_tmp = dat$Contingency[dat$sub_i == n]
  # scale values
  Run_tmp = (Run_tmp-min(Run_tmp))/max(Run_tmp-min(Run_tmp))
  Contingency_tmp = (Contingency_tmp-min(Contingency_tmp))/max(Contingency_tmp-min(Contingency_tmp))
  # organize together
  N_trials_valid[n] = sum(dat$sub_i == n)
  RT_table[1:N_trials_valid[n],n] = RT_tmp
  Run_table[1:N_trials_valid[n],n] = Run_tmp 
  Contingency_table[1:N_trials_valid[n],n] = Contingency_tmp 
}

stan_data = list(
  Run = Run_table,
  Cue = mean(dat$CueOnset),
  theta_b0 = -3.1, # fix the proportion of early onset at pnorm(theta_b0) at the first run
  RT = RT_table,
  Contingency = Contingency_table,
  N_trials = N_trials_max,
  N_trials_valid = N_trials_valid,
  N_runs = N_runs,
  N_subjects = N_subjects
)

# STAN model - both contingency levels: write model ----

write("// Stan model mixture gaus model - training with random theta
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
      real<lower=0, upper=Cue-100> mu1_fix; // Anticipatory mean
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
      }", 
      stan_file)

# STAN model - both contingency levels: run ----

#fit_stan2 <- stan(file = stan_file, data = stan_data, iter = 2000, chains = 4,seed = 1, control = list(adapt_delta =.8))
#save(list = c("fit_stan2", "stan_data"), file = paste0(analysis_path,"fit_stan2.Rda"))
load(file = paste0(analysis_path,"fit_stan2.Rda"))
print(fit_stan2)
traceplot(fit_stan2, inc_warmup = FALSE)

# STAN model - both contingency levels: simulate posterior data ----
params = c()
for (chain in 1:4){
  params = rbind(params,as.data.frame(rstan::extract(fit_stan2, permuted=FALSE, inc_warmup = FALSE)[,chain,]))
}
chains_length <- length(params$lp__)
N_subjects = stan_data$N_subjects
N_trials = stan_data$N_trials

simulations = data.frame(Sub = rep(1:N_subjects,rep(N_trials,N_subjects)),
                         Trial = rep(1:N_trials,N_subjects),
                         Run = as.vector(stan_data$Run),
                         Cue = as.vector(stan_data$Cue), 
                         Contingency = as.vector(stan_data$Contingency), 
                         RT = as.vector(stan_data$RT)
)

simulations = simulations[simulations$RT<=10000,]
n_sims = nrow(simulations)
simulations$sim_i = sample(1:chains_length,n_sims,replace = TRUE)
simulations$theta_b0 = stan_data$theta_b0
simulations$theta_b1 = NA
for (sub_i in 1:N_subjects){
  theta_b1_varname = paste0("theta_b1_rand[",sub_i,"]")
  theta_b2_varname = paste0("theta_b2_rand[",sub_i,"]")
  mu1_b0_varname = paste0("mu1_rand[",sub_i,"]")
  mu2_b0_varname = "mu2_fix"
  
  simulations$theta_b1[simulations$Sub == sub_i] = params[simulations$sim_i[simulations$Sub == sub_i],theta_b1_varname]
  simulations$theta_b2[simulations$Sub == sub_i] = params[simulations$sim_i[simulations$Sub == sub_i],theta_b2_varname]
  simulations$mu_1[simulations$Sub == sub_i] = params[simulations$sim_i[simulations$Sub == sub_i],mu1_b0_varname]
  simulations$mu_2[simulations$Sub == sub_i] = params[simulations$sim_i[simulations$Sub == sub_i],mu2_b0_varname]
  
  cat("\rprogress:",round(100*sub_i/N_subjects),"%")
  flush.console() 
}
simulations$Contingency2 = as.factor(simulations$Contingency)
levels(simulations$Contingency2) = c("50%", "100%")
simulations$theta = simulations$theta_b0 + (simulations$theta_b1*simulations$Contingency + (1-simulations$Contingency)*simulations$theta_b2)*simulations$Run
simulations$p = pnorm(simulations$theta)
simulations$mu_2 = params$mu2_fix[simulations$sim_i]
simulations$sigma_1 = params$`sigma_e[1]`[simulations$sim_i]
simulations$sigma_2 = params$`sigma_e[2]`[simulations$sim_i]
simulations$selection = rbinom(n_sims,1,simulations$p)
simulations$RT_simulation = rnorm(n_sims, simulations$mu_1 ,simulations$sigma_1) * simulations$selection +
  rnorm(n_sims, simulations$mu_2 ,simulations$sigma_2) * (1-simulations$selection)
hist(simulations$RT_simulation - simulations$Cue,100)

simulations$RT_simulation_eff = simulations$RT_simulation - simulations$Cue
simulations$RT_eff = simulations$RT - simulations$Cue
simulations$Run2 = as.factor(simulations$Run*19 + 1)
simulations$selection_fac= as.factor(simulations$selection)

# Actual Data
ggplot(data = simulations) +
  geom_density(aes(x =  RT,y = ..density.., color = Run2)) + 
  facet_grid(. ~ Contingency2)
# Posterior Distribution
ggplot(data = simulations) +
  geom_density(aes(x =  RT_simulation,y = ..density.., color = Run2)) + 
  facet_grid(. ~ Contingency2)

ggplot(data = simulations) + geom_density(aes(x = RT, y = ..density..,color = 'Original Data')) + 
  geom_density(aes(x = RT_simulation, y = ..density..,color = 'Posterior Simulation')) + 
  facet_grid(. ~ Contingency2)

ggplot(data = simulations) + geom_density(aes(x = RT_eff, y = ..density..,color = 'Original Data')) + 
  geom_density(aes(x = RT_simulation_eff, y = ..density..,color = 'Posterior Simulation')) + 
  scale_x_continuous(name = "RT - from cue onset") + 
  facet_grid(. ~ Contingency2)


# Simulated posterior - split by Gaussians: count
ggplot() +
  geom_density(data = subset(simulations, selection ==0), aes(x =  RT_simulation,y = ..count.., color = Run2)) + 
  geom_density(data = subset(simulations, selection ==1), aes(x =  RT_simulation,y = ..count.., color = Run2))

ggplot() +
  geom_density(data = subset(simulations, selection ==0), aes(x =  RT_eff,y = ..count.., color = Run2)) 
ggplot() +
  geom_density(data = subset(simulations, selection ==1), aes(x =  RT_eff,y = ..count.., color = Run2))

# STAN model - both contingency levels: correlation with choice ----
params %>% select(contains(c("mu1_rand"))) %>% summarise_all(mean) %>%
  gather(Sub, Mu1, contains("mu1")) 

gather(Sub, Theta1, contains("theta_b1")) %>%
  gather(Sub, Theta2, contains("theta_b2")) 


model_summary = as.data.frame(summary(fit_stan2)$summary) 
theta1_mean = model_summary$mean[grep("theta_b1_rand",row.names(model_summary))]
theta2_mean = model_summary$mean[grep("theta_b2_rand",row.names(model_summary))]
mu1_mean = model_summary$mean[grep("mu1_rand",row.names(model_summary))]

probe_mean = tapply(Probe_data$Outcome, list(Probe_data$sub_i, Probe_data$Contingency), mean, na.rm = T)
DM_dat = melt(probe_mean, varnames = c("sub","Contingency"), value.name = "probe_effect")
DM_dat$theta = c(theta2_mean,theta1_mean)
DM_dat$Contingency2 = as.factor(DM_dat$Contingency)
levels(DM_dat$Contingency2) = c("50% Contingency", "100% Contingency")
# new - tidyr version
DM_dat = Probe_data %>% group_by(Contingency, sub_i) %>%
  summarise(probe_effect = mean(Outcome, na.rm=TRUE)) %>% group_by() %>%
  mutate(theta = c(theta2_mean,theta1_mean),
         Contingency2 = factor(Contingency, labels = c("50% Contingency", "100% Contingency"))) %>%
  as.data.frame()

ggplot(data = DM_dat, aes(x = theta, y = probe_effect, color = Contingency2)) +
  labs(x = expression(paste(theta," slope")), y = 'proportion of trials Go stimuli were chosen') +
  geom_smooth(method = 'lm',se = TRUE, fullrange=TRUE) +
  #geom_smooth(method = 'glm', method.args = list(family = "binomial"),se = TRUE, fullrange=TRUE) +
  geom_point() +
  theme_bw() + 
  geom_abline(slope = 0, intercept = .5, linetype = 2) +
  facet_grid(. ~ Contingency2, scales = "free_x") + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1), expand = c(0,0)) + 
  theme(legend.position = "none")



DM_dat %>% group_by(Contingency) %>% summarise(mean_theta = mean(theta), sd_theta = sd(theta))

summary(lm.beta(lm(probe_effect ~ theta, data = subset(DM_dat,Contingency==0.5))))
summary(lm.beta(lm(probe_effect ~ theta, data = subset(DM_dat,Contingency==1))))
summary(lmer(theta ~ 1 + Contingency2 + (1|sub), data = DM_dat))

# STAN model - both contingency levels and correlation: data organization ====
stan_file = "training_mixture_both_contingency_levels_regression.stan"

subs2test = 1:100 # set to examine less participants for speed
dat = subset(CAT_data_GO, sub_i %in% subs2test)
N_trials_per_sub = as.vector(tapply(dat$RT, list(dat$sub_i), length))
N_trials_max = max(N_trials_per_sub)
N_runs = max(dat$Run)
N_subjects = nlevels(factor(dat$sub_i))

RT_table = matrix(data = 999000,nrow = N_trials_max,ncol = N_subjects)
Run_table = RT_table
Contingency_table = RT_table
N_trials_valid = c()
Prop_chose_go_100 = c()
Prop_chose_go_50 = c()

for (n in 1:N_subjects){
  # find values
  RT_tmp = dat$RT[dat$sub_i == n]
  Run_tmp = dat$Run[dat$sub_i == n]
  Contingency_tmp = dat$Contingency[dat$sub_i == n]
  # scale values
  Run_tmp = (Run_tmp-min(Run_tmp))/max(Run_tmp-min(Run_tmp))
  Contingency_tmp = (Contingency_tmp-min(Contingency_tmp))/max(Contingency_tmp-min(Contingency_tmp))
  # organize together
  N_trials_valid[n] = sum(dat$sub_i == n)
  RT_table[1:N_trials_valid[n],n] = RT_tmp
  Run_table[1:N_trials_valid[n],n] = Run_tmp 
  Contingency_table[1:N_trials_valid[n],n] = Contingency_tmp 
  # Proportion of trials Go was chosen
  Prop_chose_go_100[n] = mean(Probe_data$Outcome[Probe_data$sub_i==n & Probe_data$Contingency==1],na.rm=T)
  Prop_chose_go_50[n] = mean(Probe_data$Outcome[Probe_data$sub_i==n & Probe_data$Contingency==.5],na.rm=T)
}

stan_data = list(
  Run = Run_table,
  Cue = mean(dat$CueOnset),
  theta_b0 = -3.1, # fix the proportion of early onset at pnorm(theta_b0) at the first run
  RT = RT_table,
  Contingency = Contingency_table,
  N_trials = N_trials_max,
  N_trials_valid = N_trials_valid,
  N_runs = N_runs,
  N_subjects = N_subjects,
  Prop_chose_go_100 = Prop_chose_go_100,
  Prop_chose_go_50 = Prop_chose_go_50
)

# STAN model - both contingency levels and correlation: write model ----
write("// Stan model mixture gaus model - training with random theta and correlation
      data {
      int N_subjects;
      int N_trials;
      int N_trials_valid[N_subjects];
      real Cue;
      real theta_b0;
      real RT[N_trials, N_subjects];
      real Run[N_trials, N_subjects];
      real Contingency[N_trials, N_subjects];
      real Prop_chose_go_100[N_subjects];
      real Prop_chose_go_50[N_subjects];
      }
      
      parameters {
      // Fixed effect: mu1(Anticipatory), mu2(Cue dependent), theta1(100% contingency), theta2(50% contingency), sigma_e1, sigma_e2
      real<lower=0, upper=Cue-100> mu1_fix; // Anticipatory mean
      real<lower=Cue> mu2_fix; // Cue dependent mean
      real theta_fix[2]; // theta0 + slope (for run)
      real<lower=0> sigma_e[4]; // error terms
      real b0[2];
      real b1[2];
      
      // SD for random effect: mu1, theta1, theta2
      real<lower=0> mu1_sd[1]; // SD of the random effects around the fixed mean 
      real<lower=0> theta_sd[2]; // SD of the random effects around the fixed thetas
      
      // Random effect parameter: mu1, theta1, theta2
      real <upper=Cue-100> mu1_rand[N_subjects]; // random intercept for early onest Gaussian's mean
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
      target += normal_lpdf(Prop_chose_go_100[s] | b0[1] + b1[1]*theta_b1_rand[s] , sigma_e[3]);
      target += normal_lpdf(Prop_chose_go_50[s] | b0[2] + b1[2]*theta_b2_rand[s] , sigma_e[4]);

      for (t in 1:N_trials_valid[s]){
      target += log_mix(Phi_approx(theta_b0 + (Contingency[t,s]*theta_b1_rand[s] + (1-Contingency[t,s])*theta_b2_rand[s]) * Run[t,s]),
      normal_lpdf(RT[t,s] | mu1_rand[s] , sigma_e[1]),
      normal_lpdf(RT[t,s] | mu2_fix, sigma_e[2]));
      }
      }
      }", 
      stan_file)

# STAN model - both contingency levels: run ----

fit_stan_corr <- stan(file = stan_file, data = stan_data, iter = 2000, chains = 4,seed = 1, control = list(adapt_delta =.8))
# save(list("fit_stan_corr","stan_data"), file = paste0(analysis_path,"fit_stan_corr.Rda"))
load(file = paste0(analysis_path,"fit_stan_corr.Rda"))
print(fit_stan_corr)
traceplot(fit_stan2, inc_warmup = FALSE)

# STAN model - both contingency levels: simulate posterior data [NOT WORKING YET!]----
params = c()
for (chain in 1:4){
  params = rbind(params,as.data.frame(rstan::extract(fit_stan_corr, permuted=FALSE, inc_warmup = FALSE)[,chain,]))
}
chains_length <- length(params$lp__)
N_subjects = stan_data$N_subjects
N_trials = stan_data$N_trials

simulations = data.frame(Sub = rep(1:N_subjects,rep(N_trials,N_subjects)),
                         Trial = rep(1:N_trials,N_subjects),
                         Run = as.vector(stan_data$Run),
                         Cue = as.vector(stan_data$Cue), 
                         Contingency = as.vector(stan_data$Contingency), 
                         RT = as.vector(stan_data$RT)
)

simulations = simulations[simulations$RT<=10000,]
n_sims = nrow(simulations)
simulations$sim_i = sample(1:chains_length,n_sims,replace = TRUE)
simulations$theta_b0 = stan_data$theta_b0
simulations$theta_b1 = NA
for (sub_i in 1:N_subjects){
  theta_b1_varname = paste0("theta_b1_rand[",sub_i,"]")
  theta_b2_varname = paste0("theta_b2_rand[",sub_i,"]")
  mu1_b0_varname = paste0("mu1_rand[",sub_i,"]")
  mu2_b0_varname = "mu2_fix"
  
  simulations$theta_b1[simulations$Sub == sub_i] = params[simulations$sim_i[simulations$Sub == sub_i],theta_b1_varname]
  simulations$theta_b2[simulations$Sub == sub_i] = params[simulations$sim_i[simulations$Sub == sub_i],theta_b2_varname]
  simulations$mu_1[simulations$Sub == sub_i] = params[simulations$sim_i[simulations$Sub == sub_i],mu1_b0_varname]
  simulations$mu_2[simulations$Sub == sub_i] = params[simulations$sim_i[simulations$Sub == sub_i],mu2_b0_varname]
  
  cat("\rprogress:",round(100*sub_i/N_subjects),"%")
  flush.console() 
}
simulations$Contingency2 = as.factor(simulations$Contingency)
levels(simulations$Contingency2) = c("50%", "100%")
simulations$theta = simulations$theta_b0 + (simulations$theta_b1*simulations$Contingency + (1-simulations$Contingency)*simulations$theta_b2)*simulations$Run
simulations$p = pnorm(simulations$theta)
simulations$mu_2 = params$mu2_fix[simulations$sim_i]
simulations$sigma_1 = params$`sigma_e[1]`[simulations$sim_i]
simulations$sigma_2 = params$`sigma_e[2]`[simulations$sim_i]
simulations$selection = rbinom(n_sims,1,simulations$p)
simulations$RT_simulation = rnorm(n_sims, simulations$mu_1 ,simulations$sigma_1) * simulations$selection +
  rnorm(n_sims, simulations$mu_2 ,simulations$sigma_2) * (1-simulations$selection)
hist(simulations$RT_simulation - simulations$Cue,100)

simulations$RT_simulation_eff = simulations$RT_simulation - simulations$Cue
simulations$RT_eff = simulations$RT - simulations$Cue
simulations$Run2 = as.factor(simulations$Run*19 + 1)
simulations$selection_fac= as.factor(simulations$selection)

# Actual Data
ggplot(data = simulations) +
  geom_density(aes(x =  RT,y = ..density.., color = Run2)) + 
  facet_grid(. ~ Contingency2)
# Posterior Distribution
ggplot(data = simulations) +
  geom_density(aes(x =  RT_simulation,y = ..density.., color = Run2)) + 
  facet_grid(. ~ Contingency2)

ggplot(data = simulations) + geom_density(aes(x = RT, y = ..density..,color = 'Original Data')) + 
  geom_density(aes(x = RT_simulation, y = ..density..,color = 'Posterior Simulation'))

ggplot(data = simulations) + geom_density(aes(x = RT, y = ..density..,color = 'Original Data')) + 
  geom_density(aes(x = RT_simulation, y = ..density..,color = 'Posterior Simulation')) + 
  facet_grid(. ~ Contingency2)

ggplot(data = simulations) + geom_density(aes(x = RT_eff, y = ..density..,color = 'Original Data')) + 
  geom_density(aes(x = RT_simulation_eff, y = ..density..,color = 'Posterior Simulation')) + 
  scale_x_continuous(name = "RT - from cue onset") + 
  facet_grid(. ~ Contingency2)

# STAN model - both contingency levels: correlation with choice ----

model_summary = as.data.frame(summary(fit_stan_corr)$summary)
model_summary[grep("b1\\[",row.names(model_summary)),]

# STAN model - individual items (NOT WORKING YET!): data organization ====
stan_file = "training_mixture_individual_items.stan"

subs2test = 1:100 # set to smaller number to examine less participants for speed
dat = subset(CAT_data_GO, sub_i %in% subs2test)
N_trials_per_sub = dat %>% group_by(sub_i) %>% summarise(n()) %>% pull
N_trials_max = max(N_trials_per_sub)
N_runs = max(dat$Run)
N_subjects = nlevels(factor(dat$sub_i))
N_stim = dat %>% filter(sub_i ==1, Contingency==1) %>% pull(StimNum) %>% as.factor() %>% nlevels()
RT_table = matrix(data = 999000,nrow = N_trials_max,ncol = N_subjects)
Run_table = RT_table
Contingency_table = RT_table
Stim_table = RT_table
N_trials_valid = c()
for (n in 1:N_subjects){
  # find values
  RT_tmp = dat$RT[dat$sub_i == n]
  Run_tmp = dat$Run[dat$sub_i == n]
  Contingency_tmp = dat$Contingency[dat$sub_i == n]
  Stim_tmp = dat$StimNum[dat$sub_i == n]
  
  # scale values
  Run_tmp = (Run_tmp-min(Run_tmp))/max(Run_tmp-min(Run_tmp))
  Contingency_tmp = (Contingency_tmp-min(Contingency_tmp))/max(Contingency_tmp-min(Contingency_tmp))
  Stim_tmp[Contingency_tmp==0] = as.numeric(as.factor(Stim_tmp[Contingency_tmp==0]))
  Stim_tmp[Contingency_tmp==1] = as.numeric(as.factor(Stim_tmp[Contingency_tmp==1]))
  # organize together
  N_trials_valid[n] = sum(dat$sub_i == n)
  RT_table[1:N_trials_valid[n],n] = RT_tmp
  Run_table[1:N_trials_valid[n],n] = Run_tmp 
  Contingency_table[1:N_trials_valid[n],n] = Contingency_tmp 
  Stim_table[1:N_trials_valid[n],n] = Stim_tmp 
}

stan_data = list(
  Run = Run_table,
  Cue = mean(dat$CueOnset),
  theta_b0 = -3.1, # fix the proportion of early onset at pnorm(theta_b0) at the first run
  RT = RT_table,
  Contingency = Contingency_table,
  Stim = Stim_table,
  N_trials = N_trials_max,
  N_trials_valid = N_trials_valid,
  N_runs = N_runs,
  N_subjects = N_subjects,
  N_stim = N_stim
)
