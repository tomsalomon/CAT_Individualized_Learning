# Initialization ----

# Load required R packages ----
library(rstudioapi) # to define current working directory
library(ggplot2)
library(lme4)
library(lmerTest)
library(rstan)
library(lm.beta)
library(reshape2)
library(tidyverse)
library(gifski)
library(kableExtra)
# library(janitor)
library(gganimate) # If you want to save animation plots
library(magick)# If you want to save animation plots
library(ggpubr)
library(bayesplot)

# Define R environment ----
# clear workspace
rm(list=ls())
# Define the current script location as the working directory
pwd = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(pwd)
models_path = paste0(pwd,'/StanModels')
figures_path = paste0(pwd,'/Figures/')

# Load data from appropriate dir location (change if not saved in the same directory as the plotting script)
data_path = pwd # Change according to the path where the data is saved
load(paste0(data_path,"/probe_data.Rda"))
load(paste0(data_path,"/training_data.Rda"))
exp2test = unique(training_data$Experiment) # change for debugging
options(mc.cores = parallel::detectCores()) # parallel computing for Rstan

# Custom functions ----
# Standard error function
sem = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }
# count_unique function
count_unique = function(x) {length(unique(x))}
# extract_stan function
extract_stan = function(x) {
  params = c()
  n_chains = dim(rstan::extract(x, permuted=FALSE))[2]
  for (chain in 1:n_chains){
    params = rbind(params,as.data.frame(rstan::extract(x, permuted=FALSE)[,chain,]))
  }
  return(params)
}
expit = function(u){exp(u)/(1+exp(u))}

# Descriptive statistics ----
# count number of subjects and training runs per experiment
training_data %>% group_by(Experiment, ExperimentName) %>% 
  summarise(n = count_unique(subjectID), Runs = max(runNum)) %>% 
  kbl() %>% kable_minimal()

# Filter only Go trials within training, with sensible RT
training_data2 = subset(training_data, !is.na(RT_eff) & RT >= 100 & RT <=1500) # Only Go trials within training

# Missing data
training_data %>% filter(WasCue==1) %>%
  summarise(n = n(),
            missed = mean(is.na(RT) | RT >1500, na.rm=TRUE), 
            early_response = mean(RT<100, na.rm=TRUE), 
            excluded = mean(is.na(RT) | RT >1500 | RT<100, na.rm=TRUE),
            n_excluded = n*excluded) %>%
  kbl() %>% kable_minimal()

exclusion_total = training_data %>% filter(WasCue==1) %>%
  summarise(n = n(),
            missed = mean(is.na(RT) | RT >1500, na.rm=TRUE), 
            early_response = mean(RT<100, na.rm=TRUE), 
            excluded = mean(is.na(RT) | RT >1500 | RT<100, na.rm=TRUE),
            n_excluded = n*excluded) %>%
  mutate(ExperimentName = "Total")

training_data %>% filter(WasCue==1) %>% group_by(ExperimentName)%>%
  summarise(n = n(),
            missed = mean(is.na(RT) | RT >1500, na.rm=TRUE), 
            early_response = mean(RT<100, na.rm=TRUE), 
            excluded = mean(is.na(RT) | RT >1500 | RT<100, na.rm=TRUE),
            n_excluded = n*excluded) %>%
  rbind(exclusion_total) %>%
  mutate(missed = paste0(round(missed*100,2),"%"),
         early_response = paste0(round(early_response*100,2),"%"),
         excluded = paste0(round(excluded*100,2),"%")) %>%
  kbl() %>% kable_minimal() %>% row_spec(max(exp2test)+1, bold = T)

RT_summary_by_run = training_data2 %>% group_by(runNum) %>% 
  summarise(RT_eff_mean = mean(RT_eff),
            RT_eff_sd = sd(RT_eff), 
            RT_eff_1_per = quantile(RT_eff,0.01),
            RT_eff_2.5_per = quantile(RT_eff,0.025),
            RT_eff_50_per = quantile(RT_eff,0.5),
            RT_eff_97.5_per = quantile(RT_eff,0.975),
            RT_eff_99_per = quantile(RT_eff,0.99),
            prop_anticipatory = mean(RT_eff<145)) # use the 1% in Run 1 (145ms) as the threshold for anticipatory response
RT_summary_by_run %>% kbl %>% kable_minimal()

RT_summary_by_run %>% filter(runNum ==1 | runNum==20) %>% 
  pivot_longer(cols = RT_eff_mean:prop_anticipatory) %>% 
  # gather(key = runNum, value = value) %>%
  mutate(runNum = paste0("Run_",runNum)) %>%
  spread(key = runNum, value = value) %>%
  kbl %>% kable_minimal()

training_data2$anticipatory = training_data2$RT_eff<RT_summary_by_run$RT_eff_1_per[1]
training_data2 %>% group_by(anticipatory, runNum) %>% 
  summarise(n = n(),
            RT_eff_mean = mean(RT_eff),
            RT_eff_sd = sd(RT_eff), 
            RT_eff_1_per = quantile(RT_eff,0.01),
            RT_eff_2.5_per = quantile(RT_eff,0.025),
            RT_eff_50_per = quantile(RT_eff,0.5),
            RT_eff_97.5_per = quantile(RT_eff,0.975),
            RT_eff_99_per = quantile(RT_eff,0.99)) %>% 
  kbl %>% kable_minimal()

# Animation plot ----
p = ggplot(data = subset(training_data2)) +
  theme_bw() + 
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) + # trim space from x-axis
  geom_histogram(aes(x=RT, y = ..density.., fill = runNum2, group = 1L), bins = 100) +
  geom_density(aes(x = RT, y = ..density.., color = runNum2, group = 1L), color = "black") + 
  geom_vline(xintercept = 1000, linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2) + 
  annotate(geom = "text", x = 250, y = .0065,label = "Stimulus Onset") + 
  annotate(geom = "text", x = 1250, y = .0065,label = "Stimulus Offset") +
  labs(color = "Training\nRun",fill = "Training\nRun", x = "RT (from trial onset)") + 
  theme(text = element_text(size=12), legend.text =  element_text(size=8),legend.key.size = unit(0.1,"inch"))

plot(p)
anim <- p +
  transition_states(runNum2,
                    transition_length = 2,
                    state_length = 2) + 
  ggtitle('Run number: {closest_state} of 20')
# will take a while to save as a GIF animation - uncomment to run again
# animate(anim,  height = 500, width = 700, units = "in", res = 300, renderer = magick_renderer())
# anim_save(filename = paste0(figures_path,"Training_RT_animation.gif"),animation = last_animation(), type = 'cairo') 

# RT plot (time locked to trial onset) ----
RT_plot = ggplot(data = subset(training_data2), aes(x=RT, color = runNum2)) +
  theme_bw() + 
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) + # trim space from x-axis
  geom_density(size = 0.3) + 
  geom_vline(xintercept = 1000, linetype = 2, size = 0.3) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) + 
  annotate(geom = "text", x = 180, y = .0075,label = "Stimulus Onset", size = 2) + 
  annotate(geom = "text", x = 1180, y = .0075,label = "Stimulus Offset", size = 2) +
  labs(color = "Training\nRun", x = "RT (from trial onset)") + 
  ggtitle("") +
  theme(text = element_text(size=8), legend.text =  element_text(size=8),legend.key.size = unit(0.1,"inch"))
# RT_plot
# dev.copy2pdf(file=paste0(figures_path,"MetaAnalysis_RT_density_locked_to_trial_onset.pdf"), width = 5, height = 3)

# RT plot (time locked to Cue) ----
RT_effective_plot = ggplot(data = subset(training_data2), aes(x=RT_eff, color = runNum2)) +
  theme_bw() + 
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) + # trim space from x-axis
  geom_density(size = 0.3) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) + 
  annotate(geom = "text", x = 180, y = .0065,label = "Cue Onset", size = 2) +  
  labs(color = "Training\nRun", x = "RT effective (from Cue onset)") + 
  ggtitle("") +
  theme(text = element_text(size=8), legend.text =  element_text(size=8),legend.key.size = unit(0.1,"inch"))  
# RT_effective_plot
# dev.copy2pdf(file=paste0(figures_path,"MetaAnalysis_RT_density_locked_to_cue_onset.pdf"), width = 5, height = 3)

# RT plot (time locked to Cue) - split to late and early runs
n_run_groups = 4
run_groups_width = 20/n_run_groups
RT_effective_grouped_plot = training_data2 %>% 
  # filter(runNum>=18 | runNum<=3) %>%
  #mutate(late_run = factor(runNum>=17,labels = c('Early runs','Late runs'))) %>%
  mutate(run_Q = cut(runNum, breaks=n_run_groups, include.lowest = TRUE, labels = FALSE),
         late_run = paste0("Q",run_Q, " (Runs ",(run_Q-1)*run_groups_width+1," - ",(run_Q)*run_groups_width,")")) %>%
  ggplot(aes(x=RT_eff, color = runNum2)) +
  geom_density(size = 0.3) + 
  facet_grid(. ~ late_run) + 
  theme_bw() + 
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) + # trim space from x-axis
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) + 
  annotate(geom = "text", x = 300, y = .0065,label = "Cue Onset", size = 2) + 
  labs(color = "Training\nRun", x = "RT effective (from Cue onset)") + 
  ggtitle("") +
  theme(text = element_text(size=8), legend.text =  element_text(size=8),legend.key.size = unit(0.1,"inch"))  
# RT_effective_grouped_plot
# dev.copy2pdf(file=paste0(figures_path,"MetaAnalysis_RT_density_locked_to_cue_onset_run_groups.pdf"), width = 8, height = 3)

# Plot 3 panels into one figure
ggarrange(ggarrange(RT_plot, RT_effective_plot,
                    ncol = 2,
                    common.legend = TRUE,
                    font.label = list(size = 10),
                    labels = c("a.", "b."),
                    legend = "none"), 
          RT_effective_grouped_plot,
          font.label = list(size = 10),
          labels = c("","c."),
          common.legend = TRUE, nrow = 2, legend = "right")
dev.copy2pdf(file=paste0(figures_path,"MetaAnalysis_RT_density_all.pdf"), width = 7, height = 4)

# Plot probe results ----
means = probe_data %>% subset(PairType<=1) %>% 
  group_by(subjectID2, Experiment, ExperimentName) %>% 
  summarise(CAT_effect_sub=mean(Outcome,na.rm=T))%>% 
  group_by(Experiment, ExperimentName) %>%
  summarise(CAT_effect=mean(CAT_effect_sub), SE = sem(CAT_effect_sub)) %>% 
  mutate(ymin = CAT_effect - SE, ymax = CAT_effect + SE) %>%
  group_by() %>%
  mutate(Stimulus = ifelse(Experiment %in% c(6, 10:14), "Faces",
                           ifelse(Experiment %in% c(1,7,15),"Fractals",
                                  ifelse(Experiment %in% c(2,8),"Positive affective",
                                         ifelse(Experiment %in% c(3,9),"Negative affective",
                                                "Snacks")))))

ggplot(data = means, aes(x = Experiment, y = CAT_effect, fill = Stimulus)) + 
  geom_bar(stat = 'identity', color = 'black') +
  scale_y_continuous(limits = c(0,1), breaks=seq(0, 1, 0.1),expand = c(0,0),  name = "Proportion of trials high-value Go stimuli were chosen") +
  scale_x_continuous(limits = c(0,nrow(means)+1), breaks=seq(1, nrow(means), 1),expand = c(0,0),  name = "Experiment") +
  geom_errorbar(data = means, aes(ymin=ymin, ymax=ymax,x = Experiment), width=.1) + 
  theme_bw() + 
  theme(legend.position="top") + 
  geom_abline(intercept = (0.5),slope=0,linetype =2)  # chance level 50% reference line
dev.copy2pdf(file=paste0(figures_path,"MetaAnalysis_Probe_results_all_experiments.pdf"), width = 7, height = 5)

# Stan model: effective RT - FINAL MODEL ----
stan_file = paste0(models_path,"/training_mixture_RT_eff.stan")

dat = subset(training_data2, Experiment %in% exp2test)
dat$sub_i2 = as.numeric(factor(dat$subjectID2))
N_trials_per_sub = dat %>% group_by(sub_i2) %>% summarise(N=n()) %>% pull(N) 
N_trials = max(N_trials_per_sub)
N_runs = max(dat$runNum)
N_subjects = nlevels(factor(dat$sub_i2))

RT_table = matrix(data = 0,nrow = N_trials,ncol = N_subjects)
Cue_table = RT_table
Run_table = RT_table
N_trials_valid = c()
for (n in 1:N_subjects){
  RT_tmp = dat$RT_eff[dat$sub_i2 == n]
  Run_tmp = dat$runNum[dat$sub_i2 == n]
  N_trials_valid[n] = sum(dat$sub_i2 == n)
  RT_table[1:N_trials_valid[n],n] = RT_tmp
  Run_table[1:N_trials_valid[n],n] = Run_tmp  
}
Run_table2 = (Run_table-1)/19

stan_data = list(
  N_subjects = nlevels(factor(dat$subjectID2)),
  N_trials = N_trials,
  N_trials_valid = N_trials_per_sub,
  RT = RT_table,
  Run = Run_table2,
  theta_b0 = -3.1
)

write("// Stan model mixture gaus model - training with random theta
      data {
      int N_subjects;
      int N_trials;
      int N_trials_valid[N_subjects];
      real RT[N_trials, N_subjects];
      real Run[N_trials, N_subjects];
      real theta_b0;
      }
      
      parameters {
      real <upper=150> mu_fix_1; // fixed intercept - 2 distributions' means
      real <lower=200> mu_fix_2; // fixed intercept - 2 distributions' means
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
      ", 
      stan_file)

###  Stan model: run model ----
# Run the commented line to regenerate the stan model (will take a few good hours) or load trained Stan model
# Fit_RT_eff <- stan(file = stan_file, data = stan_data, iter = 2000, chains = 4,seed = 1)
# save(list = c("Fit_RT_eff", "stan_data"), file = paste0(models_path,"/Fit_RT_eff.Rda"))
load(file = paste0(models_path,"/Fit_RT_eff.Rda"))
print(Fit_RT_eff)
model_summary = summary(Fit_RT_eff)$summary %>% as.data.frame()
model_summary[1:6,] %>% select(mean, `2.5%`, `97.5%`) %>% 
  kbl(digits = 2) %>% kable_minimal()

###  Stan model: Trace plot ----
vars2pars_fix = pars = c("mu_fix_1","mu_fix_2", "sigma_e[1]", 
                         "sigma_e[2]", "theta_fix", "theta_sd")
vars2pars_rand = pars = c("theta_b1[1]", "theta_b1[2]", "theta_b1[3]", 
                          "theta_b1[4]", "theta_b1[5]", "theta_b1[6]")
vars2pars_name_fix = c("mu[1]", "mu[2]","sigma[epsilon[1]]", 
                       "sigma[epsilon[2]]", "theta[slope]","sigma[theta[slope]]")
vars2pars_name_rand = c("theta[slope[1]]", "theta[slope[2]]", "theta[slope[3]]", 
                        "theta[slope[4]]", "theta[slope[5]]", "theta[slope[6]]")
mcmc_fix = as.array(Fit_RT_eff)[,,vars2pars_fix]
mcmc_rand = as.array(Fit_RT_eff)[,,vars2pars_rand]
dimnames(mcmc_fix)[[3]] = vars2pars_name_fix
dimnames(mcmc_rand)[[3]] = vars2pars_name_rand
color_scheme_set("mix-brightblue-darkgray")

traceplot_fix = mcmc_trace(mcmc_fix, facet_args = list(labeller = ggplot2::label_parsed), size = 0.05,) + 
  theme(text = element_text(size=8), strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.spacing.y = unit(0, "lines"),
        strip.text = element_text(size=10, face = "bold"))+ 
  scale_x_continuous(breaks = seq(0, 1000, 500))+
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  ggtitle("")
traceplot_rand = mcmc_trace(mcmc_rand, facet_args = list(labeller = ggplot2::label_parsed), size = 0.05,) + 
  theme(text = element_text(size=8), strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.spacing.y = unit(0, "lines"),
        strip.text = element_text(size=10, face = "bold"))+ 
  scale_x_continuous(breaks = seq(0, 1000, 500))+
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  ggtitle("")
ggarrange(traceplot_fix, NULL, traceplot_rand, widths = c(1, 0.1, 1, 0.05), font.label = list(size = 10),
          labels = c("a.", "", "b.",""),
          common.legend = TRUE, ncol = 4, nrow = 1, legend = "bottom")
dev.copy2pdf(file=paste0(figures_path,"MetaAnalysis_StanTraceplot.pdf"), width = 7, height = 4)

###  Stan model: correlation with choice ----
model_summary = as.data.frame(summary(Fit_RT_eff)$summary)
theta_mean = model_summary$mean[grep("theta_b1",row.names(model_summary))]
probe_dat =  probe_data %>% filter(Experiment %in% exp2test , PairType<=2) %>%
  mutate(sub_i2 = as.numeric(factor(subjectID2)))
DM_dat = probe_dat %>% group_by(sub_i2) %>%
  summarise(probe_effect = mean(Outcome, na.rm=T),
            chose_go = sum(Outcome==1, na.rm=T),
            chose_nogo = sum(Outcome==0, na.rm=T),
            exp = as.factor(mean(Experiment))) %>%
  mutate(exp_num = as.numeric(exp),
         theta = theta_mean)

# lmer model (USED IN THE PAST, LESS APPROPRIATE THAN LOGISTIC REGRESSION BELOW)
# summary(lmer(scale(probe_effect) ~ 1 + scale(theta) + (1 + scale(theta)||exp), data = subset(DM_dat)))

# unscaled glmer model - to make predictions
glmer_model = glmer(cbind(chose_go,chose_nogo) ~ 1 + theta + (1 + theta|exp),data = subset(DM_dat), family = "binomial")
summary(glmer_model)
glmer_model_summary = summary(glmer_model)$coefficients %>% as.data.frame()

# Skip built-in CI function (takes a long time)
# PE = data.frame(Estimate = fixef(glmer_model))
# CI = confint(glmer_model, parm = row.names(PE))
# glmer_model_params = cbind(PE,CI)
# colnames(glmer_model_params) = c("Estimate", "CI_lower", "CI_upper")
# glmer_model_params
# exp(glmer_model_params)

# Confidence interval using normal distribution assumption
glmer_model_se = diag(vcov(glmer_model))^0.5
glmer_X = model.matrix(glmer_model)
glmer_model_beta + qnorm(0.025)*glmer_model_se # Alternative estimation with ~N
glmer_model_params = data.frame(Estimate =  fixef(glmer_model)) %>%
  mutate(CI_lower = Estimate + qnorm(0.025)*glmer_model_se,
         CI_upper = Estimate + qnorm(0.975)*glmer_model_se,
         OR = exp(Estimate),
         OR_CI_lower = exp(CI_lower),
         OR_CI_upper = exp(CI_upper),
         Z = glmer_model_summary$`z value`,
         p = glmer_model_summary$`Pr(>|z|)`)
print(glmer_model_params)

# Decision making data - merge probe and training
DM_dat = DM_dat %>% 
  mutate(pred_fixed = expit(glmer_X%*%glmer_model_params$Estimate),
         pred_fixed_lower = expit(glmer_X%*%(glmer_model_params$CI_lower)),
         pred_fixed_upper = expit(glmer_X%*%(glmer_model_params$CI_upper)),
         pred = predict(glmer_model,DM_dat, type = "response"),
         pred2 = predict(glmer_model,DM_dat, type = "response",re.form = NA)) # same as pred_fixed
SE_dat = data.frame(theta = c(sort(DM_dat$theta),sort(DM_dat$theta,decreasing =TRUE)), 
                    probe_effect = c(sort(DM_dat$pred_fixed_lower), sort(DM_dat$pred_fixed_upper,decreasing =TRUE)))

###  Stan model: Plot correlation with probe ------
# visualization - logistic regression model
ggplot(data = subset(DM_dat, theta < 200 ), aes(x = theta, y = probe_effect)) +
  theme_bw() + 
  labs(x = expression(theta[slope[i]]), y = 'proportion of trials Go stimuli were chosen') +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(alpha = 0.2, size = 2)  +
  geom_point(alpha = 0.8, size = 2, shape = 21)  + 
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  geom_polygon(data = SE_dat, alpha = 0.4, fill = "firebrick2") +
  geom_line(aes(x = theta, y = pred_fixed_lower), color = "black", size = .5, linetype =3) + 
  geom_line(aes(x = theta, y = pred_fixed_upper), color = "black", size = .5, linetype =3) +
  geom_line(aes(x = theta, y = pred_fixed), color = "black", size = 1) + 
  theme(text = element_text(size=12))
dev.copy2pdf(file=paste0(figures_path,"MetaAnalysis_CorrelationWithChoice.pdf"), width = 5, height = 4)

###  Stan model: simulate posterior data ----
params = extract_stan(Fit_RT_eff)
chains_length = length(params$lp__)
N_subjects = stan_data$N_subjects
N_trials = stan_data$N_trials
N_runs = stan_data$N_trials

simulations = data.frame(Sub = rep(1:N_subjects,rep(N_trials,N_subjects)),
                         Trial = rep(1:N_trials,N_subjects),
                         Run = as.vector(stan_data$Run),
                         RT_eff = as.vector(stan_data$RT)
)
simulations = simulations[simulations$Run>=0,]
n_sims = nrow(simulations)

set.seed(0)
simulations = simulations %>% 
  mutate(sim_i = sample(1:chains_length,n_sims,replace = TRUE),
         theta_b0 = stan_data$theta_b0,
         theta_b1 = NA)
for (sub_i in 1:N_subjects){
  theta_b1_varname = paste0("theta_b1[",sub_i,"]")
  simulations$theta_b1[simulations$Sub == sub_i] = params[simulations$sim_i[simulations$Sub == sub_i],theta_b1_varname]
  cat("\rprogress:",round(100*sub_i/N_subjects),"%")
  flush.console() 
}

set.seed(0)
simulations = simulations %>%
  mutate(theta =theta_b0 +theta_b1*Run,
         p = pnorm(theta),
         mu_1 = params$`mu_fix_1`[sim_i],
         mu_2 = params$`mu_fix_2`[sim_i],
         sigma_1 = params$`sigma_e[1]`[sim_i],
         sigma_2 = params$`sigma_e[2]`[sim_i],
         selection = rbinom(n_sims,1,p),
         RT_simulation_eff = rnorm(n_sims, mu_1 ,sigma_1) * selection +
           rnorm(n_sims, mu_2 ,sigma_2) * (1-selection),
         Run2 = as.factor(simulations$Run*19 + 1),
         selection_fac= as.factor(selection)) %>%
  group_by(Sub) %>% mutate(theta_b1_mean = mean(theta_b1)) %>%
  group_by(.) %>% 
  mutate(theta_b1_group = cut(theta_b1_mean, 
                              breaks = quantile(theta_b1_mean, probs = seq(0, 1, length.out = 5)),
                              include.lowest=TRUE))
levels(simulations$theta_b1_group) = paste0("Q", 1:nlevels(simulations$theta_b1_group),
                                            ": ", levels(simulations$theta_b1_group))

simulations_long = simulations %>%
  mutate(`Actual RT` = RT_eff, `Posterior Simulation` = RT_simulation_eff) %>%
  pivot_longer(cols = c("Actual RT","Posterior Simulation"), names_to = "RT measurement", values_to="RT")

###  Stan model: simulation plots ----
# Posterior simulations versus RT (by theta group and run)
ggplot(data = simulations_long, aes(x=RT, color = Run2)) +
  facet_grid(`RT measurement` ~ theta_b1_group) + 
  theme_bw() + 
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) + # trim space from x-axis
  geom_density(size = 0.3) + 
  geom_vline(xintercept = 0, linetype = 2) + 
  #annotate(geom = "text", x = 200, y = .0065,label = "Cue\nOnset") +  
  labs(color = "Training\nRun", x = "RT effective (from Cue onset)", y = "Density") + 
  theme(text = element_text(size=12), axis.text = element_text(size=6), legend.text = element_text(size=8),legend.key.size = unit(0.1,"inch"))  
dev.copy2pdf(file=paste0(figures_path,"MetaAnalysis_PosteriorSimulationVSActualData.pdf"), width = 7, height = 4)

# Posterior simulations versus RT (by theta group only)
simulations_long %>% filter(Run2 %in% c(1, seq(2,20,2))) %>%
  mutate(Run3 = sprintf("Run %02i",Run2)) %>%
  ggplot(aes(x=RT, color = `RT measurement`)) +
  facet_grid(Run3 ~ theta_b1_group) + 
  theme_bw() + 
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) + # trim space from x-axis
  geom_density(size = 0.3) + 
  geom_vline(xintercept = 0, linetype = 2) + 
  #annotate(geom = "text", x = 200, y = .0065,label = "Cue\nOnset") +  
  labs(x = "RT effective (from Cue onset)", y = "Density") + 
  theme(text = element_text(size=10),axis.text = element_text(size=6), legend.position = "top")  
dev.copy2pdf(file=paste0(figures_path,"MetaAnalysis_PosteriorSimulationVSActualData_matrix.pdf"), width = 7, height = 7)

# Stan model - individual stimuli ====
stan_file =  paste0(models_path,"training_mixture_individual_stimuli.stan")

# Some go stimuli were used in sanity trials. Exclude these from analysis
valid_go_stim = probe_data %>% filter(PairType <=2)  %>% 
  mutate(stim_go = ifelse(IsleftGo,as.character(ImageLeft),as.character(ImageRight)),
         sub_stim_i_name = sprintf("%02s_%02s",subjectID2,stim_go)) %>%
  pull(sub_stim_i_name) %>% unique()
training_data_valid_stim = training_data2 %>% 
  mutate(sub_stim_i_name = sprintf("%02s_%02s",subjectID2,itemName)) %>%
  filter(sub_stim_i_name %in% valid_go_stim)

# # Select sub-sample for speed testing
# n2select = 10 # how many participants to randomly sample from each experiment (set to low number for debugging)
# n_exp = count_unique(training_data_valid_stim$Experiment)
# set.seed(0)
# dat = c()
# for (exp_i in 1:n_exp){
#   dat_tmp = training_data_valid_stim %>% filter(Experiment == exp_i)
#   n_subs_i = count_unique(dat_tmp$sub_i)
#   selected_subs = sample(1:n_subs_i, size = n2select, replace = FALSE)
#   dat = rbind(dat, filter(dat_tmp, sub_i %in% selected_subs)) 
# }
dat = training_data_valid_stim
dat = dat %>% mutate(
  sub_i2 = as.numeric(as.factor(as.character(dat$subjectID2)))) %>%
  group_by(sub_i2) %>%
  mutate(StimNum = as.numeric(as.factor(as.character(itemName)))) %>% group_by() %>%
  mutate(sub_stim_i_name = sprintf("%02s_%02s",subjectID2,itemName),
         sub_stim_i = sprintf("%02i,%02i",sub_i2,StimNum))

N_trials_per_sub = dat %>% group_by(sub_i2) %>% summarise(n()) %>% pull
N_trials_max = max(N_trials_per_sub)
N_runs = max(dat$runNum)
N_subjects = count_unique(dat$sub_i2)
RT_table = matrix(data = 999000,nrow = N_trials_max,ncol = N_subjects)
Run_table = RT_table
Stim_table = RT_table
Cue_table = RT_table
N_trials_valid = c()
N_stim_valid = c()
prev_stim_count = 0
Stim_df = c()
for (n in 1:N_subjects){
  # find values
  RT_tmp = dat$RT_eff[dat$sub_i2 == n]
  Run_tmp = dat$runNum[dat$sub_i2 == n]
  Stim_tmp = dat$StimNum[dat$sub_i2 == n] 
  # Cue_tmp = dat$CueTime[dat$sub_i2 == n]
  
  # scale values
  Run_tmp = (Run_tmp-1)/19
  # organize together
  N_trials_valid[n] = sum(dat$sub_i2 == n)
  RT_table[1:N_trials_valid[n],n] = RT_tmp
  Run_table[1:N_trials_valid[n],n] = Run_tmp 
  Stim_table[1:N_trials_valid[n],n] = Stim_tmp + prev_stim_count
  # Cue_table[1:N_trials_valid[n],n] = Cue_tmp 
  N_stim_valid[n] = max(Stim_tmp)
  Stim_df_tmp = data.frame(sub_stim_i = sort(unique(Stim_tmp + prev_stim_count)),
                           sub_i2 = n,
                           stim_i = 1:max(Stim_tmp),
                           Experiment = unique(dat$Experiment[dat$sub_i2 == n]))
  Stim_df = rbind(Stim_df, Stim_df_tmp)
  prev_stim_count = prev_stim_count + N_stim_valid[n]
}
N_stim = max(Stim_df$sub_stim_i)

stan_data = list(
  Run = Run_table,
  theta_b0 = -3.1, # fix the proportion of early onset at pnorm(theta_b0) at the first run
  RT = RT_table,
  Stim = Stim_table,
  N_trials = N_trials_max,
  N_trials_valid = N_trials_valid,
  N_runs = N_runs,
  N_subjects = N_subjects,
  N_stim = N_stim
)

write("
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
	// Participant-level (random-effect) parameters: individualized learning parameters
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
", 
stan_file)

### Stan model - individual stimuli: run model ----

# fit_stan_individual_stimuli = stan(file = stan_file, data = stan_data, iter = 4000, chains = 4, seed = 0, control = list(adapt_delta =.85, max_treedepth=12)) 
# save(list = c("fit_stan_individual_stimuli","stan_data", "dat"), file = paste0(models_path,"/fit_stan_individual_stimuli.Rda"))
# print(fit_stan_individual_stimuli)

# For this part you have to have the full trained model (> 1G in size).
# Since GitHub has a limit of 100M, please contact the first author to get this file (or rerun the model for ~24h)
## PARTS yoU WILL NEED THE 1G FILE TO RECREATE :
# load(file = paste0(models_path,"/fit_stan_individual_stimuli.Rda"))
# model_summary_stimuli = summary(fit_stan_individual_stimuli)$summary %>% as.data.frame()
# vars2pars_fix = c("mu1_fix","mu2_fix", "sigma_e[1]",
#                          "sigma_e[2]", "theta_slope_fix", "theta_stim_sd")
# vars2pars_rand = c("theta_slope_stim[1]", "theta_slope_stim[2]", "theta_slope_stim[3]",
#                           "theta_slope_stim[4]", "theta_slope_stim[5]", "theta_slope_stim[6]")
# mcmc_fix = as.array(fit_stan_individual_stimuli)[,,vars2pars_fix]
# mcmc_rand = as.array(fit_stan_individual_stimuli)[,,vars2pars_rand]
# 
# ###  Stan model: simulate posterior data ----
# 
# # For this part you will need the 1G data file. Use previously made simulation instead
# params = extract_stan(fit_stan_individual_stimuli)
# chains_length = length(params$lp__)
# N_subjects = stan_data$N_subjects
# N_trials = stan_data$N_trials
# N_stim = stan_data$N_stim
# 
# simulations = data.frame(Sub = rep(1:N_subjects,rep(N_trials,N_subjects)),
#                          Run = as.vector(stan_data$Run),
#                          RT_eff = as.vector(stan_data$RT),
#                          Stim = as.vector(stan_data$Stim)
#                          
# )
# simulations = simulations[simulations$Run<99,]
# n_sims = nrow(simulations)
# set.seed(0)
# simulations = simulations %>% 
#   mutate(sim_i = sample(1:chains_length,n_sims,replace = TRUE),
#          theta_b0 = stan_data$theta_b0,
#          theta_b1 = NA)
# for (stim_i in 1:N_stim){
#   theta_b1_varname = paste0("theta_slope_stim[",stim_i,"]")
#   simulations$theta_b1[simulations$Stim == stim_i] = params[simulations$sim_i[simulations$Stim == stim_i],theta_b1_varname]
#   cat("\rprogress:",round(100*stim_i/N_stim),"%")
#   flush.console() 
# }
# 
# set.seed(0)
# simulations = simulations %>%
#   mutate(theta = theta_b0 +theta_b1*Run,
#          p = pnorm(theta),
#          mu_1 = params$mu1_fix[sim_i],
#          mu_2 = params$mu2_fix[sim_i],
#          sigma_1 = params$`sigma_e[1]`[sim_i],
#          sigma_2 = params$`sigma_e[2]`[sim_i],
#          selection = rbinom(n_sims,1,p),
#          RT_simulation_eff = rnorm(n_sims, mu_1 ,sigma_1) * selection +
#            rnorm(n_sims, mu_2 ,sigma_2) * (1-selection),
#          Run2 = as.factor(round(simulations$Run*19,0) + 1),
#          selection_fac= as.factor(selection)) %>%
#   group_by(Stim) %>% mutate(theta_b1_mean = mean(theta_b1)) %>%
#   group_by(.) %>% 
#   mutate(theta_b1_group = cut(theta_b1_mean, 
#                               breaks = quantile(theta_b1_mean, probs = seq(0, 1, length.out = 5)),
#                               include.lowest=TRUE))
# levels(simulations$theta_b1_group) = paste0("Q", 1:nlevels(simulations$theta_b1_group),
#                                             ": ", levels(simulations$theta_b1_group))

###  Stan model - individual stimuli: load pre processed large data ----
# For compatibility with GitHub 100M size limit, large data structures are pre processed and stored 
# save(list = c("dat","model_summary_stimuli","mcmc_fix", "mcmc_rand"), file = paste0(models_path,"/fit_stan_individual_stimuli_summary.Rda")) # Save the required data from the large Stan model result
load(paste0(models_path,"/fit_stan_individual_stimuli_summary.Rda"))

model_summary_stimuli[1:6,] %>% select(mean, `2.5%`, `97.5%`) %>% 
  kbl(digits = 2) %>% kable_minimal()

###  Stan model - individual stimuli: Trace plot ----
vars2pars_name_fix = c("mu[1]", "mu[2]","sigma[epsilon[1]]", 
                       "sigma[epsilon[2]]", "theta[slope]","sigma[theta[slope]]")
vars2pars_name_rand = c("theta[slope[list(1,1)]]", "theta[slope[list(1,2)]]", "theta[slope[list(1,3)]]", 
                        "theta[slope[list(1,4)]]", "theta[slope[list(1,5)]]", "theta[slope[list(1,6)]]")
dimnames(mcmc_fix)[[3]] = vars2pars_name_fix
dimnames(mcmc_rand)[[3]] = vars2pars_name_rand
color_scheme_set("mix-brightblue-darkgray")

traceplot_fix = mcmc_trace(mcmc_fix, facet_args = list(labeller = ggplot2::label_parsed), size = 0.05,) + 
  theme(text = element_text(size=8), strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.spacing.y = unit(0, "lines"),
        strip.text = element_text(size=10, face = "bold"))+ 
  scale_x_continuous(breaks = seq(0, 2000, 1000))+
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  ggtitle("")
traceplot_rand = mcmc_trace(mcmc_rand, facet_args = list(labeller = ggplot2::label_parsed), size = 0.05,) + 
  theme(text = element_text(size=8), strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.spacing.y = unit(0, "lines"),
        strip.text = element_text(size=10, face = "bold"))+ 
  scale_x_continuous(breaks = seq(0, 2000, 1000))+
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  ggtitle("")

ggarrange(traceplot_fix, NULL, traceplot_rand, widths = c(1, 0.1, 1, 0.05), font.label = list(size = 10),
          labels = c("a.", "", "b.",""),
          common.legend = TRUE, ncol = 4, nrow = 1, legend = "bottom")
dev.copy2pdf(file=paste0(figures_path,"MetaAnalysis_StanTraceplot_IndividualStimuli.pdf"), width = 7, height = 4)

### Stan model - individual stimuli: correlation with choice ----
theta_ind = grep(c("theta_slope_stim\\["),row.names(model_summary_stimuli))
thetas = Stim_df %>% 
  mutate(theta_sub_stim = model_summary_stimuli$mean[theta_ind],
         sub_i = sub_i2,
         sub_stim_i = sprintf("%02i,%02i",sub_i,stim_i)) %>% group_by(sub_i) %>%
  mutate(theta_sub_stim_mean = mean(theta_sub_stim)) %>% ungroup()

thetas_sub_model = DM_dat %>%  mutate(sub_i = sub_i2, theta_sub = theta) %>% # thetas from previous Stan model
  select(sub_i,theta_sub) 
stim_names = dat %>% select("sub_stim_i","sub_stim_i_name") %>%
  arrange(sub_stim_i) %>% unique()
thetas2 = thetas %>% left_join(thetas_sub_model, by = "sub_i") %>%
  left_join(stim_names, by = "sub_stim_i") %>%
  select(sub_stim_i_name, theta_sub, theta_sub_stim)
selected_subs_all = dat %>% pull(subjectID2) %>% unique()

DM_dat_long = probe_data %>% filter(subjectID2 %in% selected_subs_all) %>% filter(PairType <=2)  %>% 
  mutate(stim_go = ifelse(IsleftGo,as.character(ImageLeft),as.character(ImageRight)),
         sub_i2 = as.numeric(as.factor(as.character(subjectID2)))) %>% 
  group_by(sub_i2) %>%
  mutate(StimNum = as.numeric(as.factor(as.character(stim_go)))) %>% group_by() %>%
  mutate(sub_stim_i_name = sprintf("%02s_%02s",subjectID2,stim_go),
         sub_stim_i = sprintf("%02i,%02i",sub_i2,StimNum)) %>%
  left_join(thetas2, by = "sub_stim_i_name") %>% mutate(theta = theta_sub_stim)

# plot theta (subjects model) versus stimulus-thets
thetas2 %>% ggplot(aes(x = theta_sub, y = theta_sub_stim)) +
  geom_point(alpha = 0.4) + 
  geom_smooth(method = "lm", formula = "y ~ x ") + 
  theme_bw()

DM_dat_stimuli = DM_dat_long %>%
  group_by(Experiment,sub_stim_i, sub_i2, subjectID2) %>%
  summarise(probe_effect = mean(Outcome, na.rm=T),
            chose_go = sum(Outcome==1, na.rm=T),
            chose_nogo = sum(Outcome==0, na.rm=T),
            theta = mean(theta),
            theta_sub = mean(theta_sub)) %>%
  group_by(.) %>%
  mutate(theta_sub_scaled = scale(theta_sub),
         theta_scaled = scale(theta))

# Decision making models (might take a few minutes to run)
glmer_model = glmer(cbind(chose_go,chose_nogo) ~ 1 + theta + 
                      (1 + theta|Experiment/subjectID2),
                    data = DM_dat_stimuli, family = "binomial")
summary(glmer_model)
# For model comparison use nested model
glmer_model_theta_sub_and_stim = glmer(cbind(chose_go,chose_nogo) ~ 1 + theta_sub + theta + 
                                         (1 + theta |Experiment/subjectID2)+
                                         (0 + theta_sub|Experiment), 
                                       data = DM_dat_stimuli, family = "binomial")
glmer_model_theta_sub_and_stim %>% summary() %>% coefficients() %>%
  as.data.frame() %>%
  mutate(OR = exp(Estimate),
         OR_lower = exp(Estimate + `Std. Error`*qnorm(0.025)),
         OR_upper = exp(Estimate + `Std. Error`*qnorm(0.975)),
         p_one_sided = `Pr(>|z|)`/2
         )
glmer_model_theta_sub = glmer(cbind(chose_go,chose_nogo) ~ 1 + theta_sub + 
                                (1 + theta_sub|Experiment), 
                              data = DM_dat_stimuli, family = "binomial")
# Compare nested models
anova(glmer_model_theta_sub_and_stim, glmer_model_theta_sub)
glmer_model_X = model.matrix(glmer_model)
glmer_model_beta = fixef(glmer_model)
glmer_model_se = diag(glmer_model_X%*%vcov(glmer_model)%*%t(glmer_model_X))^0.5

## placeholder if not scaled:
sd_theta=1
m_theta=0

# Estimate the parameters - Confidence interval using normal distribution assumption
X_params = as.matrix(rbind(c(1,0), # intercept
                           c(0,1/sd_theta))) # theta 100% contrast
params_se = diag(X_params%*%vcov(glmer_model)%*%t(X_params))^0.5
params_CI = data.frame(#exper = rep(1,nrow(X_params)),
  Parameter = c("intercept", "theta 100%"),
  Estimate = X_params%*%glmer_model_beta,
  SE = params_se) %>%
  mutate(Z = Estimate/SE,
         p = pnorm(abs(Z),lower.tail = FALSE),
         p_two_sided = pnorm(abs(Z),lower.tail = FALSE)*2,
         OR = exp(Estimate),
         OR_CI_lower = exp(Estimate + qnorm(0.025)*SE),
         OR_CI_upper = exp(Estimate + qnorm(0.975)*SE))
params_CI
# Prediction line for all theta values in range
theta_100 = DM_dat_stimuli %>% pull(theta)
# extend the prediction line by a bit (will determine x-axis limits)
expander_100 = (max(theta_100) - min(theta_100))/50
theta_100_pred = seq(min(theta_100-expander_100),max(theta_100+expander_100), length.out =length(theta_100))
theta_100_scaled = (theta_100_pred-m_theta)/sd_theta
Pred_X = cbind(rep(1,nrow(glmer_model_X)), # intercept
               theta_100_scaled) %>% # theta
  as.matrix()
Pred_se = diag(Pred_X%*%vcov(glmer_model)%*%t(Pred_X))^0.5 # prediction variance = X*VCOV*X'
Pred = DM_dat_stimuli %>% 
  mutate(
    theta_scaled = Pred_X[,2],
    theta = theta_scaled*sd_theta + m_theta,
    pred = expit(Pred_X%*%glmer_model_beta),
    pred_lower = expit(Pred_X%*%glmer_model_beta + qnorm(0.025)*Pred_se),
    pred_upper = expit(Pred_X%*%glmer_model_beta + qnorm(0.975)*Pred_se)) %>%
  as.data.frame()

# a data frame which will be used to draw a CI polygon
Pred_SE = data.frame(theta =  c(Pred$theta,Pred$theta[nrow(Pred):1]),
                     probe_effect = c(Pred$pred_lower, Pred$pred_upper[nrow(Pred):1])) %>%
  as.data.frame()

### Stan model - individual stimuli: Correlation with choice plot ----
ggplot(data = subset(DM_dat_stimuli), aes(x = theta, y = probe_effect)) +
  ggtitle("") +
  theme_bw() + 
  labs(x = expression(theta[slope[list(i,s)]]), y = 'proportion of trials Go stimuli were chosen') +
  geom_hline(yintercept = 0.5, linetype = 2) +
  scale_y_continuous(limits = c(-0.05,1.05), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  # 95% CI - (min to max range)
  geom_polygon(data = Pred_SE, alpha = 0.3, size=0, fill = "red") +
  geom_line(aes(x = theta, y = pred_lower), data = Pred, color = "black", size = .5, linetype =3) + 
  geom_line(aes(x = theta, y = pred_upper), data = Pred, color = "black", size = .5, linetype =3) +
  geom_line(aes(x = theta, y = pred), data = Pred, color = "black", size = .5) +
  # Points
  set.seed(0) + geom_jitter(height = 0.01, alpha = .2, size = 2, stroke = .05)  + 
  set.seed(0) + geom_jitter(height = 0.01, alpha = 1, size = 2, shape = 21,  stroke = .05)  + 
  theme(text = element_text(size=8), legend.position = "top")

###  Stan model: simulation plots ----
simulations_long = simulations %>%
  mutate(`Actual RT` = RT_eff, `Posterior Simulation` = RT_simulation_eff) %>%
  pivot_longer(cols = c("Actual RT","Posterior Simulation"), names_to = "RT measurement", values_to="RT")

# Posterior simulations versus RT (by theta group and run)
ggplot(data = simulations_long, aes(x=RT, color = Run2)) +
  facet_grid(`RT measurement` ~ theta_b1_group) + 
  theme_bw() + 
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) + # trim space from x-axis
  geom_density(size = 0.3) + 
  geom_vline(xintercept = 0, linetype = 2) + 
  #annotate(geom = "text", x = 200, y = .0065,label = "Cue\nOnset") +  
  labs(color = "Training\nRun", x = "RT effective (from Cue onset)", y = "Density") + 
  theme(text = element_text(size=12), axis.text = element_text(size=6), legend.text = element_text(size=8),legend.key.size = unit(0.1,"inch"))  
dev.copy2pdf(file=paste0(figures_path,"MetaAnalysis_PosteriorSimulationVSActualData_IndividualStimuli.pdf"), width = 7, height = 4)
