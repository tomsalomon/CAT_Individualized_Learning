
library(lme4)
library(rstudioapi)    
library(rstan)
library(lmerTest)
library(dplyr)
library(tidyr)
library(kableExtra)
library(ggpubr)
library(bayesplot)
library(devEMF)

theme_set(theme_pubr())
rm(list=ls())

# Get current path
pwd = dirname(rstudioapi::getActiveDocumentContext()$path)
models_path = paste0(pwd,'/StanModels/')
figures_path = paste0(pwd,'/Figures/')
setwd(pwd)
# Assisting functions ----
# Save as PDF and EMF
save_plot = function(plot, file, width = 7, height = 4){
  emf(file = paste0(file,".emf"),width = width, height = height)
  print(plot)
  dev.off()
  pdf(file = paste0(file,".pdf"),width = width, height = height)
  print(plot)
  dev.off()
}
# Expit function (log-odds to p) 
expit = function(x){exp(x)/(1+exp(x))}
# p-value to asterisk
p2asterisk = function(p){ifelse(p<.001, "***", ifelse(p<.01, "**", ifelse(p<.05, "*", "")))}
# Confidence interval from logistic regression model / model summary function.
CI = function (in_model, var_num = '',one_sided = F) {
  # in_model = input model; [var_num] = variable in model to analyze (optional)
  
  # if the input is a model (not model summary), convert to summary
  if (class(in_model) == "glmerMod"){in_model = summary(in_model)}
  
  coef = as.data.frame(in_model$coefficients)
  coef$p = coef$`Pr(>|z|)`
  ci_output = c()
  ci_logs = as.data.frame(cbind(coef[,1],
                                coef[,1]+qnorm(0.025)*coef[,2],
                                coef[,1]+qnorm(0.975)*coef[,2]))
  ci_OR = exp(ci_logs)
  ci_prop = expit(ci_logs)
  colnames(ci_logs) = c("logOR","logOR_CI_lower","logOR_CI_upper")
  colnames(ci_OR) = c("OR","OR_CI_lower","OR_CI_upper")
  colnames(ci_prop) = c("prop","prop_CI_lower","prop_CI_upper")
  ci_output = cbind(coef,ci_OR,ci_logs,ci_prop)
  
  if (is.numeric(var_num)) {ci_output = ci_output[var_num,]}
  if (one_sided==TRUE){
    ci_output$p = ci_output$p/2
    ci_output$`Pr(>|z|)` = ci_output$`Pr(>|z|)`/2
  }
  ci_output$asterisk = p2asterisk(ci_output$p)
  
  return (ci_output)
}

# Load data ----
load("./Merged_Data.Rda")

CAT_data$Contingency2 = factor(CAT_data$Contingency, labels =  c("0% Contingency", "50% Contingency","100% Contingency"))
CAT_data$Run2 = as.factor(CAT_data$Run)
CAT_data$Response = 1 - is.na(CAT_data$RT)
CAT_data$RT_eff = CAT_data$RT - CAT_data$CueOnset
# Data frame of Go trials with response
CAT_data_GO = subset(CAT_data, Go==1 & Response == 1 & RT>=250)
Probe_data$PairType2 = as.factor(Probe_data$PairType)
Probe_data$Contingency2 = factor(Probe_data$Contingency, labels = c("50%", "100%"))

# Missing data
exclusion_total = CAT_data %>% filter(Go==1) %>%
  summarise(n = n(),
            missed = mean(is.na(RT)), 
            early_response = mean(RT<250, na.rm=TRUE), 
            excluded = mean(is.na(RT) | RT >1500 | RT<250, na.rm=TRUE),
            n_excluded = n*excluded) 
exclusion_total%>% kbl() %>% kable_minimal()

CAT_data %>% filter(Go==1) %>% group_by(experiment)%>%
  summarise(n = n(),
            missed = mean(is.na(RT)), 
            early_response = mean(RT<250, na.rm=TRUE), 
            excluded = mean(is.na(RT) | RT >1500 | RT<250, na.rm=TRUE),
            n_excluded = n*excluded) %>%
  rbind(cbind(data.frame(experiment = "Total"), exclusion_total)) %>%
  mutate(missed = paste0(round(missed*100,2),"%"),
         early_response = paste0(round(early_response*100,2),"%"),
         excluded = paste0(round(excluded*100,2),"%")) %>%
  kbl() %>% kable_minimal() %>% row_spec(3, bold = T)

RT_summary_by_run = CAT_data_GO %>% group_by(experiment, Contingency, Run) %>% 
  summarise(RT_eff_mean = mean(RT_eff),
            RT_eff_sd = sd(RT_eff), 
            RT_eff_1_per = quantile(RT_eff,0.01),
            RT_eff_2.5_per = quantile(RT_eff,0.025),
            RT_eff_50_per = quantile(RT_eff,0.5),
            RT_eff_97.5_per = quantile(RT_eff,0.975),
            RT_eff_99_per = quantile(RT_eff,0.99),
            prop_anticipatory = mean(RT_eff<145)) # use the 1% in Run 1 (145ms) as the threshold for anticipatory response
RT_summary_by_run %>% kbl %>% kable_minimal()

CAT_data_GO %>% filter(Run ==1 | Run==20) %>% 
  group_by(Run, experiment, Contingency) %>% 
  summarise(mean(RT_eff), sd(RT_eff)) %>% kbl(digits = 2) %>% kable_minimal()

RT_summary_by_run %>% filter(Run ==1 | Run==20) %>% 
  pivot_longer(cols = RT_eff_mean:prop_anticipatory) %>% 
  # gather(key = Run, value = value) %>%
  mutate(Run = paste0("Run_",Run)) %>%
  spread(key = Run, value = value) %>%
  kbl(digits = 2) %>% kable_minimal()

CAT_data_GO$anticipatory = CAT_data_GO$RT_eff<RT_summary_by_run$RT_eff_1_per[1]
CAT_data_GO %>% group_by(anticipatory, Run) %>% 
  summarise(n = n(),
            RT_eff_mean = mean(RT_eff),
            RT_eff_sd = sd(RT_eff), 
            RT_eff_1_per = quantile(RT_eff,0.01),
            RT_eff_2.5_per = quantile(RT_eff,0.025),
            RT_eff_50_per = quantile(RT_eff,0.5),
            RT_eff_97.5_per = quantile(RT_eff,0.975),
            RT_eff_99_per = quantile(RT_eff,0.99)) %>% 
  kbl(digits = 2) %>% kable_minimal()

# RT distribution plot ----
density_plot = CAT_data_GO %>% ggplot() +
  facet_grid(experiment ~ Contingency2 ) +
  geom_density(aes(x = RT_eff, y = ..density.., color = Run2),size = 0.3) +
  labs(color = "Training\nRun", x = "RT effective (from cue onset)") +
  theme_bw() +
  annotate(geom="text", x=-850+600, y=.0047, label="Cue onset\n(850ms)", size = 2.5) +
  annotate(geom = "text", x = -850+250, y = .0077,label = "Stim. onset",size = 2.5) + 
  annotate(geom = "text", x = -850+1250, y = .0077,label = "Stim. offset\n(1000ms)", size = 2.5) +  
  geom_vline(xintercept = -850+ 0, linetype =2) +
  geom_vline(xintercept = -850+850, linetype =2) +
  geom_vline(xintercept = -850+1000, linetype =2) +
  scale_y_continuous(limits =c(0,.0085), expand = c(0,0)) + 
  scale_x_continuous(limits =c(-850+(-50),-850+1500), expand = c(0,0),
                     breaks=seq(-800, 600, 200)) +   
  #theme(legend.key.size=unit(.4,"cm"),panel.spacing = unit(2, "lines")) 
  theme(text = element_text(size=10), legend.text =  element_text(size=8),legend.key.size = unit(0.1,"inch"))
density_plot
save_plot(density_plot,file=paste0(figures_path,"Study2_RT_density_locked_to_trial_onset"), width = 6, height = 4)

# Stan model - both contingency levels ----
# Load models
load(file = paste0(models_path,"fit_stan2_exp1.Rda"))
fit_stan_exp1 = fit_stan2
stan_data_exp1 = stan_data
load(file = paste0(models_path,"fit_stan2_exp2.Rda"))
fit_stan_exp2 = fit_stan2
stan_data_exp2 = stan_data
model_summary_exp1 = summary(fit_stan_exp1)$summary %>% as.data.frame() %>% mutate(Exp = 1)
model_summary_exp2 = summary(fit_stan_exp2)$summary %>% as.data.frame() %>% mutate(Exp = 2)
stan_data_both = list(stan_data_exp1, stan_data_exp2)
model_summary_both = list(model_summary_exp1,model_summary_exp2)

# Stan model - Individual stimuli: Load models
## Exp. 2 individual stimuli models exceed the 100M file limit of github. save the model summary only instead
# load(file = paste0(models_path,"/fit_stan3_exp1.Rda")) # 100% contingency
# fit_stan_exp1_stim100 = fit_stan3
# load(file = paste0(models_path,"/fit_stan4_exp1.Rda")) # 50% contingency
# fit_stan_exp1_stim50 = fit_stan4
# load(file = paste0(models_path,"/fit_stan3_exp2.Rda")) # 100% contingency
# fit_stan_exp2_stim100 = fit_stan3
# load(file = paste0(models_path,"/fit_stan4_exp2.Rda")) # 50% contingency
# fit_stan_exp2_stim50 = fit_stan4
# model_summary_exp1_50 = summary(fit_stan_exp1_stim50)$summary %>% as.data.frame() %>% mutate(Exp = 1)
# model_summary_exp2_50 = summary(fit_stan_exp2_stim50)$summary %>% as.data.frame() %>% mutate(Exp = 2)
# model_summary_exp1_100 = summary(fit_stan_exp1_stim100)$summary %>% as.data.frame() %>% mutate(Exp = 1)
# model_summary_exp2_100 = summary(fit_stan_exp2_stim100)$summary %>% as.data.frame() %>% mutate(Exp = 2)
# save(list = c("model_summary_exp1_50","model_summary_exp2_50","model_summary_exp1_100","model_summary_exp2_100"),
#      file = paste0(models_path,"/model_summary_stimuli_models.Rda"))
load(file = paste0(models_path,"/model_summary_stimuli_models.Rda"))
model_summary_individual_both_100 = list(model_summary_exp1_100,model_summary_exp2_100)
model_summary_individual_both_50 = list(model_summary_exp1_50,model_summary_exp2_50)

### Stan model - summary ----
vars2pars = pars = c("mu1_fix","mu2_fix", "sigma_e[1]", "sigma_e[2]", "mu1_sd[1]",
                     "theta_fix[1]", "theta_fix[2]", "theta_sd[1]", "theta_sd[2]")
model_summary = rbind(model_summary_exp1[vars2pars,c("Exp", "mean", "2.5%","97.5%")], 
                      model_summary_exp2[vars2pars,c("Exp", "mean", "2.5%","97.5%")])
model_summary %>% kbl(digits=2) %>% kable_minimal()
fit_stan_both = list(fit_stan_exp1, fit_stan_exp2)

### Stan model - Trace plot ----
vars2pars_name = c("mu[1]", "mu[2]","sigma[epsilon[1]]", "sigma[epsilon[2]]",
                   "sigma[mu[1]]", "theta[slope[100*'%\']]", "theta[slope[50*'%\']]",
                   "sigma[theta[slope100*'%\']]","sigma[theta[slope50*'%\']]")
traceplot_both = list()
for (exp_i in c(1:2)){
  mcmc_i = as.array(fit_stan_both[[exp_i]])[,,vars2pars]
  dimnames(mcmc_i)[[3]] = vars2pars_name
  color_scheme_set("mix-brightblue-darkgray")
  traceplot_i = mcmc_trace(mcmc_i, facet_args = list(labeller = ggplot2::label_parsed), size = 0.05,) + 
    theme(text = element_text(size=8), strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          panel.spacing.y = unit(0, "lines"),
          strip.text = element_text(size=10, face = "bold"))+ 
    scale_x_continuous(breaks = seq(0, 1000, 500))+
    guides(color = guide_legend(override.aes = list(size = 1.5))) +
    ggtitle("")
  traceplot_both[[exp_i]] = traceplot_i
}

trace_plot = ggarrange(traceplot_both[[1]], NULL, traceplot_both[[2]], widths = c(1, 0.1, 1, 0.05), font.label = list(size = 10),
          labels = c("a. Preliminary Exp.", "", "b. Replication Exp.",""),
          common.legend = TRUE, ncol = 4, nrow = 1, legend = "bottom")
trace_plot
save_plot(trace_plot, file=paste0(figures_path,"Study2_StanTraceplot_HyperParams.pdf"), width = 7, height = 4)

### Stan model - simulate posterior data ----
params_exp1 = params_exp2 = c()
for (chain in 1:4){
  params_exp1 = rbind(params_exp1,as.data.frame(rstan::extract(fit_stan_exp1, permuted=FALSE, inc_warmup = FALSE)[,chain,]))
  params_exp2 = rbind(params_exp2,as.data.frame(rstan::extract(fit_stan_exp2, permuted=FALSE, inc_warmup = FALSE)[,chain,]))
}
parms_both = list(params_exp1,params_exp2)
simulations_both = list(c(),c())
for (exp in 1:2){
  params = parms_both[[exp]] 
  stan_data = stan_data_both[[exp]]
  chains_length = length(params$lp__)
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
  set.seed(0)
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
  }
  
  set.seed(0)
  simulations = simulations %>%
    mutate(experiment_num = exp,
           Contingency2 = factor(Contingency, labels = c("50%", "100%")),
           theta = theta_b0 + (theta_b1*Contingency + (1-Contingency)*theta_b2)*Run,
           p = pnorm(theta),
           sigma_1 = params$`sigma_e[1]`[sim_i],
           sigma_2 = params$`sigma_e[2]`[sim_i],
           selection = rbinom(n_sims,1,p),
           RT_simulation = rnorm(n_sims, mu_1 ,sigma_1) * selection +
             rnorm(n_sims, mu_2 ,sigma_2) * (1-selection),
           RT_simulation_eff = RT_simulation - Cue,
           RT_eff = RT - Cue,
           Run2 = as.factor(simulations$Run*19 + 1),
           selection_fac= as.factor(selection)) %>%
    group_by(Sub) %>% mutate(theta_b1_mean = mean(theta_b1),
                             theta_b2_mean = mean(theta_b2)) %>%
    group_by(.) %>% 
    mutate(theta_b1_group = cut(theta_b1_mean, 
                                breaks = quantile(theta_b1_mean, probs = seq(0, 1, length.out = 5)),
                                include.lowest=TRUE),
           theta_b2_group = cut(theta_b2_mean, 
                                breaks = quantile(theta_b2_mean, probs = seq(0, 1, length.out = 5)),
                                include.lowest=TRUE))
  quantile_labs_lim_1 = quantile(round(simulations$theta_b1_mean,1), probs = seq(0, 1, length.out = 5))
  quantile_labs_lim_2 = quantile(round(simulations$theta_b2_mean,1), probs = seq(0, 1, length.out = 5))
  quantile_new_labs_1 = quantile_new_labs_2 = c()
  for (lab_i in c(1:4)) {
    quantile_new_labs_1[lab_i] = paste0("Q",lab_i," [",quantile_labs_lim_1[lab_i],", ",quantile_labs_lim_1[lab_i+1],"]")
    quantile_new_labs_2[lab_i] = paste0("Q",lab_i," [",quantile_labs_lim_2[lab_i],", ",quantile_labs_lim_2[lab_i+1],"]")
  }
  levels(simulations$theta_b1_group) = quantile_new_labs_1
  levels(simulations$theta_b2_group) = quantile_new_labs_2
  simulations = simulations %>% mutate(
    theta_b1_group = as.character(theta_b1_group),
    theta_b2_group = as.character(theta_b2_group),
    theta_group = ifelse(Contingency==1, theta_b1_group,theta_b2_group)
  )
  simulations_both[[exp]] = simulations
}

simulations = rbind(simulations_both[[1]],simulations_both[[2]])
simulations$experiment = factor(simulations$experiment_num, labels = levels(CAT_data$experiment))
simulations_long = simulations %>%
  mutate(`Actual RT` = RT, `Posterior Simulation` = RT_simulation) %>%
  pivot_longer(cols = c("Actual RT","Posterior Simulation"), names_to = "RT measurement", values_to="Reaction time")

### Stan model - simulation plots ----
plt_exp = c(1,1,2,2)
plt_contingency = c(0,1,0,1)
plot_simulation = list()
for (plt_i in 1:4){
  plot_simulation[[plt_i]] = 
    ggplot(data = subset(simulations_long, experiment_num==plt_exp[plt_i] & Contingency==plt_contingency[plt_i]), aes(x=`Reaction time`, color = Run2)) +
    facet_grid(`RT measurement` ~ theta_group) + 
    theme_bw() + 
    scale_y_continuous(expand = expansion(mult=c(0,0)), limits = c(0,0.01)) + # trim space from x-axis
    geom_density(size = 0.3) + 
    geom_vline(xintercept = 850, linetype = 2) + # Cue onset
    labs(color = "Training\nRun", x = "Reaction Time", y = "Density") + 
    ggtitle("\n")+
    theme(text = element_text(size=12), axis.text = element_text(size=6), 
          legend.text = element_text(size=8),legend.key.size = unit(0.1,"inch"),
          strip.text = element_text(size=8),
          strip.background = element_blank()) 
}

simulation_plot = ggarrange(plot_simulation[[1]] , plot_simulation[[2]],
          plot_simulation[[3]],  plot_simulation[[4]],
          font.label = list(size = 10),
          labels = c("a. Preliminary Exp.\n 50% Contingency","\n100% Contingency",
                     "b. Replication Exp.\n 50% Contingency","\n100% Contingency"),
          common.legend = TRUE, ncol = 2, nrow = 2, legend = "right",
          heights = c(1, 1), widths = c(1, 1))
simulation_plot
save_plot(simulation_plot, file=paste0(figures_path,"study2_PosteriorSimulationVSActualData"), width = 10, height = 8)

# Probe analysis ----
Probe_data_sub = Probe_data %>% group_by(experiment_num, subjectID,Contingency2) %>%
  summarise(value = mean(Outcome, na.rm=TRUE)) %>%   mutate(Contingency=Contingency2)
logistic_reg_groups_model = list()
logistic_reg_interaction_model = list()
logistic_reg_CI = list()
probe_plot = list()
probe_spaghetti_plot = list()
for (exp_i in 1:2){
  Probe_data_sub_exp_i = Probe_data_sub %>% filter(experiment_num==exp_i)
  diff = Probe_data_sub_exp_i$value[Probe_data_sub_exp_i$Contingency2=="100%"] -  Probe_data_sub_exp_i$value[Probe_data_sub_exp_i$Contingency2=="50%"]
  Probe_data_sub_exp_i$difference = rep(diff,each=2)
  # Examine each group is different from 50% chance level using no-intercept model
  probe_m1_exp_i = glmer(Outcome ~ 0 + Contingency2 + (1+Contingency2| subjectID), family = "binomial",
                         data = subset(Probe_data, experiment_num==exp_i))
  logistic_reg_groups_model[[exp_i]] = probe_m1_exp_i
  # Examine interaction (100% contingency is better than 50% contingency.
  # identical model overall, different interpretation for the slopes
  probe_m2_exp_i = glmer(Outcome ~ 1 + Contingency2 + (1+Contingency2| subjectID), family = "binomial",
                         data = subset(Probe_data, experiment_num==exp_i))
  logistic_reg_interaction_model[[exp_i]] = probe_m2_exp_i
  
  # X = rbind(c(1,0),
  #           c(1,1),
  #           c(0,1)) %>% as.matrix()
  # SE = diag(X %*% vcov(probe_m2_exp_i) %*% t(X))^0.5
  # expit(X %*% fixef(probe_m2_exp_i) + SE*qnorm(0.025))
  
  CI_exp_i = rbind(CI(probe_m1_exp_i,one_sided = TRUE),
                   CI(probe_m2_exp_i,one_sided = TRUE, var_num = 2)) 
  rownames(CI_exp_i) = c("50%", "100%", "interaction")
  CI_exp_i["interaction",c("prop","prop_CI_lower", "prop_CI_upper")] = NA
  logistic_reg_CI[[exp_i]] = CI_exp_i
  # Group summary statistics
  SummaryStats = data.frame(Contingency = factor(c(0.5,1), labels = c("50%", "100%")),
                            Mean = CI_exp_i$prop[1:2],
                            ymin = CI_exp_i$prop_CI_lower[1:2],
                            ymax = CI_exp_i$prop_CI_upper[1:2],
                            p = CI_exp_i$p[1:2],
                            asterisk = CI_exp_i$asterisk[1:2])
  h = max(SummaryStats$ymax) + 0.1
  InteractionDF = data.frame(x1 = c(1,1,2,2), y1 = c(h, h+0.02,h+0.02,h), # interaction line
                             x2 = 1.5, y2 = h+0.05, asterisk = CI_exp_i["interaction","asterisk"])
  
  ### Probe Bar Plot ----
  bar_plot = ggplot() +
    geom_bar(data = SummaryStats, aes(x = Contingency, y = Mean, fill = Contingency), stat = 'identity', color = "black", width=.5, size=.3) + 
    set.seed(1) + geom_jitter(data = Probe_data_sub_exp_i, shape=21 , aes(x = Contingency, y = value, fill= Contingency), 
                              size = 1, width = 0.15,  alpha= .4,  stroke = .2) +
    set.seed(1) + geom_jitter(data = Probe_data_sub_exp_i, shape=21 , aes(x = Contingency, y = value, group= Contingency), 
                              size = 1, width = 0.15,  alpha= 1, stroke = .2) +
    scale_y_continuous(limits = c(0,1), breaks=seq(0, 1, 0.1),expand = c(0,0)) +
    labs(x = "Go training contingency", y = "Proportion of trials Go stimuli were chosen") +
    theme_bw() +
    ggtitle("") + 
    theme(text = element_text(size=8)) + 
    geom_text(data = SummaryStats, aes(x=Contingency, y = ymax+0.05, label = asterisk), size = 5) +
    geom_errorbar(data = SummaryStats, aes(ymin=ymin, ymax=ymax,x = Contingency), width=.1,  size=.3) + 
    geom_abline(intercept = (0.5),slope=0,linetype =2) # 50% chance level reference line
  if (CI_exp_i["interaction","p"]<0.05) {
    bar_plot = bar_plot  +
      geom_path(data = InteractionDF, aes(x = x1, y = y1), size = 0.3) +
      geom_text(data = InteractionDF[1,], aes(x= x2, y = y2, label = asterisk), size = 5)
  }
  probe_plot[[exp_i]] = bar_plot
  
  ### Probe Spaghetti Plot ----
  spaghetti_plot = Probe_data_sub_exp_i %>% 
    arrange(difference, "decreasing") %>% 
    ggplot(aes(x = Contingency, y = value, color = difference, group=as.factor(difference))) +
    scale_x_discrete(expand = c(0.1,0.1)) + 
    scale_y_continuous(limits = c(0.25,1), breaks=seq(0, 1, 0.1),expand = c(0,0)) +
    labs(x = "Go training contingency", y = "Proportion of trials Go stimuli were chosen") +
    #geom_line(size = 1.5, color = "black") + 
    theme_bw() +
    theme(text = element_text(size=8)) + 
    ggtitle("") + 
    geom_hline(yintercept = 0.5, linetype = 2, size =0.5) + 
    geom_line(color = "black", size = .6, arrow = arrow(length = unit(0.1,"points"),ends = "both", type = "closed")) +
    geom_line(size = .5, arrow = arrow(length = unit(0.1,"points"),ends = "both", type = "closed")) +
    scale_color_gradient2(high = "dodgerblue2", mid = "gray90", low = "firebrick1",
                          limits = c(-.31, .31),oob = scales::squish, name = "Difference", 
                          labels = scales::label_percent(accuracy = 1))
  probe_spaghetti_plot[[exp_i]] = spaghetti_plot
}

probe_bar_plot = ggarrange(probe_plot[[1]], probe_plot[[2]], 
          font.label = list(size = 10),
          labels = c("a. Preliminary Exp.", "b. Replication Exp."),
          common.legend = TRUE, ncol = 2, nrow = 1, legend = "none")
probe_bar_plot
save_plot(probe_bar_plot, file=paste0(figures_path,"study2_Probe_barplot"), width = 5, height = 3)

probe_spaghetti_plot2save = ggarrange(probe_spaghetti_plot[[1]], probe_spaghetti_plot[[2]], 
          font.label = list(size = 8),
          labels = c("a. Preliminary Exp.", "b. Replication Exp."),
          common.legend = TRUE, ncol = 2, nrow = 1, legend = "bottom")
probe_spaghetti_plot2save
save_plot(probe_spaghetti_plot2save,file=paste0(figures_path,"study2_Probe_spaghettiplot"), width = 4, height = 3)

### Probe model summary table ----
cbind(experiment = 1, logistic_reg_CI[[1]]) %>% 
  rbind(cbind(experiment = 2, logistic_reg_CI[[1]])) %>%
  dplyr::select(experiment, prop, `z value`, p, OR, OR_CI_lower, OR_CI_upper, asterisk) %>%
  mutate(p = as.character(signif(p, 2)),
         prop = paste0(round(prop*100,2), "%")) %>%
  kbl(digits = 2) %>% kable_minimal()

### Examine value effect (make sure the effect is consistent above and beyond value difference) ----
Probe_data = Probe_data %>%
  mutate(value_diff = IsleftGo*(bidLeft - bidRight) + (1-IsleftGo)*(bidRight - bidLeft),
         value_ind_diff = IsleftGo*(bidIndexLeft - bidIndexRight) + 
           (1-IsleftGo)*(bidIndexRight - bidIndexLeft)) 
# glmer(Outcome ~ 0 + Contingency2 + value_diff + (1+value_diff+Contingency2| subjectID), family = "binomial",
#       data = subset(Probe_data, experiment_num==1)) %>% summary
# glmer(Outcome ~ 0 + Contingency2 + value_diff + (1+value_diff+Contingency2| subjectID), family = "binomial",
#       data = subset(Probe_data, experiment_num==2)) %>% summary
# glmer(Outcome ~ 1 + Contingency2 + value_diff + (1+value_diff+Contingency2| subjectID), family = "binomial",
#       data = subset(Probe_data, experiment_num==1)) %>% summary
# glmer(Outcome ~ 1 + Contingency2 + value_diff + (1+value_diff+Contingency2| subjectID), family = "binomial",
#       data = subset(Probe_data, experiment_num==2)) %>% summary

# Stan model - both contingency levels: correlation with choice ----
theta_probe_plot_both = list()
theta_probe_plot_same_axis_both = list()
theta_probe_model_both = list()
theta_probe_model_long_both = list()
params_CI_both = list()

for (exp_i in 1:2){
  model_summary = model_summary_both[[exp_i]] # stan model data
  theta_100 = model_summary$mean[grep("theta_b1_rand",row.names(model_summary))]
  theta_50 = model_summary$mean[grep("theta_b2_rand",row.names(model_summary))]
  thetas = data.frame(theta = c(theta_50,theta_100),
                      Contingency = rep(c(0.5,1),each = length(theta_50)),
                      sub_i = rep(1:length(theta_50),2)) %>%
    mutate(contingency_sub_i = paste0(Contingency*100,",",sub_i)) %>%
    select(contingency_sub_i, theta)
  
  DM_dat_long = Probe_data %>% filter(experiment_num == exp_i) %>% 
    mutate(sub_i = as.numeric(as.factor(subjectID)),
           Contingency2 = factor(Contingency, labels = c("50% Contingency","100% Contingency")),
           contingency_sub_i = paste0(Contingency*100,",",sub_i),
           stim_go = ifelse(IsleftGo,as.character(ImageLeft),as.character(ImageRight))) %>%
    left_join(thetas, by = "contingency_sub_i") %>%
    group_by(sub_i) %>%
    mutate(stim_go_i = as.numeric(as.factor(stim_go))) %>%
    group_by(.)
  
  DM_dat_med = DM_dat_long%>% # to fit the size of random slope per item
    group_by(Contingency, Contingency2, stim_go_i, sub_i, subjectID) %>%
    summarise(probe_effect = mean(Outcome, na.rm=T),
              chose_go = sum(Outcome==1, na.rm=T),
              chose_nogo = sum(Outcome==0, na.rm=T),
              theta = mean(theta)) %>%
    group_by(.) 
  
  DM_dat = DM_dat_long %>%
    group_by(Contingency, Contingency2, subjectID) %>%
    summarise(probe_effect = mean(Outcome, na.rm=T),
              chose_go = sum(Outcome==1, na.rm=T),
              chose_nogo = sum(Outcome==0, na.rm=T)) %>%
    group_by(.) %>% 
    mutate(theta = c(theta_50,theta_100),
           theta_scaled = c(scale(theta_50),scale(theta_100)))
  
  
  # GLMER model (same as modeling each condition separately)
  glmer_model = glmer(cbind(chose_go,chose_nogo) ~ 1 + theta*Contingency2 + 
                        (1 + Contingency2|subjectID), data = DM_dat, family = "binomial")
  glmer_model_no_contingency = glmer(cbind(chose_go,chose_nogo) ~ 1 + theta + 
                                       (1 + Contingency2|subjectID), data = DM_dat, family = "binomial")
  cat("\n\nExperiment", exp_i,"LRT (of contingency effect)",
      "\n=======================================\n")
  anova(glmer_model_no_contingency,glmer_model, test = "LRT") %>% print()
  cat("\nAIC difference (restricted - full):",
      AIC(glmer_model_no_contingency) - AIC(glmer_model),"\n\n",
      sep = "")
  # same model - dont group to compare with theta per individual stimuli model
  
  # glmer_model_long = glmer(cbind(chose_go,chose_nogo) ~ 1 + theta*Contingency2 + 
  #                            (1 + Contingency2|subjectID), data = DM_dat_med, family = "binomial")
  
  # glmer_model_100 = glmer(cbind(chose_go,chose_nogo) ~ 1 + theta + (1|subjectID), data = subset(DM_dat, Contingency == 1), family = "binomial")
  # glmer_model_50 = glmer(cbind(chose_go,chose_nogo) ~ 1 + theta + (1|subjectID), data = subset(DM_dat, Contingency == 0.5), family = "binomial")
  theta_probe_model_both[[exp_i]] = glmer_model
  # theta_probe_model_long_both[[exp_i]] = glmer_model_long
  glmer_model_X = model.matrix(glmer_model)
  glmer_model_beta = fixef(glmer_model)
  glmer_model_se = diag(glmer_model_X%*%vcov(glmer_model)%*%t(glmer_model_X))^0.5
  
  # Estimate the parameters - Confidence interval using normal distribution assumption
  X_params = as.matrix(rbind(c(0,1,0,0), # theta 50% contrast
                             c(0,1,0,1), # theta 100% contrast
                             c(0,0,0,1))) # interaction contrast
  params_se = diag(X_params%*%vcov(glmer_model)%*%t(X_params))^0.5
  params_CI = data.frame(exp = rep(exp_i,nrow(X_params)),
                         Parameter = c("theta 50%", "theta 100%", "interaction"),
                         Estimate = X_params%*%glmer_model_beta,
                         SE = params_se) %>%
    mutate(Z = Estimate/SE,
           p = pnorm(abs(Z),lower.tail = FALSE),
           p_two_sided = pnorm(abs(Z),lower.tail = FALSE)*2,
           OR = exp(Estimate),
           OR_CI_lower = exp(Estimate + qnorm(0.025)*SE),
           OR_CI_upper = exp(Estimate + qnorm(0.975)*SE))
  params_CI$p[params_CI$Parameter == "interaction"] = NA
  params_CI$p_two_sided[!(params_CI$Parameter == "interaction")] = NA
  
  params_CI_both[[exp_i]] = params_CI
  
  # Prediction line for all theta values in range
  expander_50 = (max(theta_50) - min(theta_50))/50
  expander_100 = (max(theta_100) - min(theta_100))/50
  theta_50_pred = seq(min(theta_50-expander_50),max(theta_50+expander_50), length.out =length(theta_50))
  theta_100_pred = seq(min(theta_100-expander_100),max(theta_100+expander_100), length.out =length(theta_100))
  
  Pred_X = cbind(rep(1,nrow(glmer_model_X)), # intercept
                 c(theta_50_pred,theta_100_pred), # theta
                 rep(c(0,1),each = length(theta_50)), # contingency
                 c(rep(0,length(theta_50)), theta_100_pred)) %>%
    as.matrix()
  Pred_se = diag(Pred_X%*%vcov(glmer_model)%*%t(Pred_X))^0.5
  
  Pred = DM_dat %>% 
    mutate(
      theta = Pred_X[,2],
      pred = expit(Pred_X%*%glmer_model_beta),
      pred_lower = expit(Pred_X%*%glmer_model_beta + qnorm(0.025)*Pred_se),
      pred_upper = expit(Pred_X%*%glmer_model_beta + qnorm(0.975)*Pred_se)) %>%
    as.data.frame()
  Pred_100 = filter(Pred, Contingency == 1)
  Pred_50 = filter(Pred, Contingency == .5)
  
  Pred_SE = data.frame(Contingency = rep(c(.5,1),each = nrow(DM_dat)),
                       theta =  c(sort(Pred_50$theta),sort(Pred_50$theta,decreasing =TRUE),
                                  sort(Pred_100$theta),sort(Pred_100$theta,decreasing =TRUE)),
                       probe_effect = c(sort(Pred_50$pred_lower), sort(Pred_50$pred_upper,decreasing =TRUE),
                                        sort(Pred_100$pred_lower), sort(Pred_100$pred_upper,decreasing =TRUE))) %>%
    mutate(Contingency2 = factor(Contingency, labels = c("50% Contingency","100% Contingency"))) %>%
    as.data.frame()
  
  DM_dat = DM_dat %>% 
    mutate(
      pred = expit(glmer_model_X%*%glmer_model_beta),
      pred_lower = expit(glmer_model_X%*%glmer_model_beta + qnorm(0.025)*glmer_model_se),
      pred_upper = expit(glmer_model_X%*%glmer_model_beta + qnorm(0.975)*glmer_model_se))
  
  DM_dat_100 = filter(DM_dat, Contingency == 1)
  DM_dat_50 = filter(DM_dat, Contingency == .5)
  SE_dat = data.frame(Contingency = rep(c(.5,1),each = nrow(DM_dat)),
                      theta =  c(sort(DM_dat_50$theta),sort(DM_dat_50$theta,decreasing =TRUE),
                                 sort(DM_dat_100$theta),sort(DM_dat_100$theta,decreasing =TRUE)),
                      probe_effect = c(sort(DM_dat_50$pred_lower), sort(DM_dat_50$pred_upper,decreasing =TRUE),
                                       sort(DM_dat_100$pred_lower), sort(DM_dat_100$pred_upper,decreasing =TRUE))) %>%
    mutate(Contingency2 = factor(Contingency, labels = c("50% Contingency","100% Contingency")))
  
  ### Correlation with choice plot ----
  corplot = ggplot(data = subset(DM_dat), aes(x = theta, y = probe_effect, color = Contingency2)) +
    facet_grid(. ~ Contingency2, scales = "free_x") +
    ggtitle("") +
    theme_bw() + 
    labs(x = expression(theta[slope[i]]), y = 'proportion of trials Go stimuli were chosen') +
    geom_hline(yintercept = 0.5, linetype = 2) +
    scale_y_continuous(limits = c(.2,1), expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    # # 95% CI (actual values)
    # geom_polygon(aes(fill = Contingency2), data = SE_dat, alpha = 0.3, size=0) +
    # geom_line(aes(x = theta, y = pred_lower), color = "black", linetype =3) + 
    # geom_line(aes(x = theta, y = pred_upper), color = "black", linetype =3) +
    # geom_line(aes(x = theta, y = pred), color = "black", size = 1) + 
    # 95% CI - (min to max range)
    geom_polygon(aes(fill = Contingency2), data = Pred_SE, alpha = 0.3, size=0) +
    geom_line(aes(x = theta, y = pred_lower, group = Contingency2), data = Pred, color = "black", size = .5, linetype =3) + 
    geom_line(aes(x = theta, y = pred_upper, group = Contingency2), data = Pred, color = "black", size = .5, linetype =3) +
    geom_line(aes(x = theta, y = pred, group = Contingency2), data = Pred, color = "black", size = .5) +
    # Points
    geom_point(aes(fill= Contingency2), color = "black", alpha = .4, size = 1, shape = 21,  stroke = .2)  + 
    geom_point(aes(group= Contingency2), color = "black", alpha = 1, size = 1, shape = 21,  stroke = .2)  + 
    theme(text = element_text(size=8), legend.position = "top")
  
  theta_probe_plot_both[[exp_i]] = corplot
  
  ### Correlation with choice plot (same axis) ----
  corplot_same_axis = ggplot(data = subset(DM_dat), aes(x = theta, y = probe_effect, color = Contingency2)) +
    facet_grid(. ~ Contingency2) +
    ggtitle("") +
    theme_bw() + 
    labs(x = expression(theta[slope[i]]), y = 'proportion of trials Go stimuli were chosen') +
    geom_hline(yintercept = 0.5, linetype = 2) +
    scale_y_continuous(limits = c(.2,1), expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    # # 95% CI (actual values)
    # geom_polygon(aes(fill = Contingency2), data = SE_dat, alpha = 0.3, size=0) +
    # geom_line(aes(x = theta, y = pred_lower), color = "black", linetype =3) + 
    # geom_line(aes(x = theta, y = pred_upper), color = "black", linetype =3) +
    # geom_line(aes(x = theta, y = pred), color = "black", size = 1) + 
    # 95% CI - (min to max range)
    geom_polygon(aes(fill = Contingency2), data = Pred_SE, alpha = 0.3, size=0) +
    geom_line(aes(x = theta, y = pred_lower, group = Contingency2), data = Pred, color = "black", size = .5, linetype =3) + 
    geom_line(aes(x = theta, y = pred_upper, group = Contingency2), data = Pred, color = "black", size = .5, linetype =3) +
    geom_line(aes(x = theta, y = pred, group = Contingency2), data = Pred, color = "black", size = .5) +
    # Points
    geom_point(aes(fill= Contingency2), color = "black", alpha = .4, size = 1, shape = 21,  stroke = .2)  + 
    geom_point(aes(group= Contingency2), color = "black", alpha = 1, size = 1, shape = 21,  stroke = .2)  + 
    theme(text = element_text(size=8), legend.position = "top")
  
  theta_probe_plot_same_axis_both[[exp_i]] = corplot_same_axis
}

theta_probe_corplot = ggarrange(theta_probe_plot_both[[1]], theta_probe_plot_both[[2]], 
          font.label = list(size = 10),
          labels = c("a. Preliminary Exp.", "b. Replication Exp."),
          common.legend = TRUE, ncol = 2, nrow = 1, legend = "none")
theta_probe_corplot
save_plot(theta_probe_corplot, file=paste0(figures_path,"study2_theta_probe_correlation"), width = 7, height = 2.5)

theta_probe_corplot_same_axis = ggarrange(theta_probe_plot_same_axis_both[[1]], theta_probe_plot_same_axis_both[[2]], 
          font.label = list(size = 10),
          labels = c("a. Preliminary Exp.", "b. Replication Exp."),
          common.legend = TRUE, ncol = 2, nrow = 1, legend = "none")
theta_probe_corplot_same_axis
save_plot(theta_probe_corplot_same_axis, file=paste0(figures_path,"study2_theta_probe_correlation_same_axis"), width = 7, height = 2.5)

### Probe correlation summary table ----
params_CI_both[[1]] %>% 
  rbind(params_CI_both[[2]]) %>%
  mutate(p  = as.character(signif(p, 2)),
         p_two_sided  = as.character(signif(p_two_sided, 2))) %>%
  kbl(digits = 2) %>% kable_minimal()

# Stan model - individual stimuli: correlation with choice ----
theta_cor_plots = list()
theta_probe_individual_plot_both = list()
theta_probe_individual_plot_both_same_axis = list()
theta_probe_model_individual_both = list()
theta_probe_model_sub_and_stim = list()
#theta_probe_model_long_individual_both = list()
params_individual_CI_both = list()
for (exp_i in 1:2){
  # previous model (slope per subject)
  model_summary = model_summary_both[[exp_i]] # stan model data
  theta_100 = model_summary$mean[grep("theta_b1_rand",row.names(model_summary))]
  theta_50 = model_summary$mean[grep("theta_b2_rand",row.names(model_summary))]
  thetas_sub = data.frame(theta_sub = c(theta_50,theta_100),
                          Contingency = rep(c(0.5,1),each = length(theta_50)),
                          sub_i = rep(1:length(theta_50),2)) %>%
    mutate(contingency_sub_i = paste0(Contingency*100,",",sub_i)) %>%
    select(contingency_sub_i, theta_sub)
  
  
  model_summary_50 = model_summary_individual_both_50[[exp_i]] # stan model data
  model_summary_100 = model_summary_individual_both_100[[exp_i]] # stan model data
  theta_ind = grep("theta_1_sub_stim",row.names(model_summary_50))
  thetas = data.frame(theta = c(model_summary_50$mean[theta_ind],model_summary_100$mean[theta_ind]),
                      ind_name = rep(row.names(model_summary_50)[theta_ind],2),
                      Contingency = rep(c(0.5,1),each = length(theta_ind))) %>%
    mutate(sub_stim_i = gsub(pattern = "theta_1_sub_stim|\\[|\\]", replacement = "", x=ind_name),
           contingency_sub_stim_i = paste0(Contingency*100, ",",sub_stim_i)) %>%
    select(contingency_sub_stim_i, theta)
  
  DM_dat_long = Probe_data %>% filter(experiment_num == exp_i) %>% 
    mutate(Contingency2 = factor(Contingency, labels = c("50% Contingency","100% Contingency")),
           stim_go = ifelse(IsleftGo,as.character(ImageLeft),as.character(ImageRight)),
           sub_i = as.numeric(as.factor(subjectID)),
           stim_go_i = NA) 
  # recode Go stimuli in each category alphabetically per participant (to match stan model indices)
  for (n in 1:max(DM_dat_long$sub_i)){
    Stim_tmp_50 = DM_dat_long %>%  filter(sub_i == n, Contingency==0.5) %>% 
      pull(stim_go) %>% as.factor() %>% as.numeric()
    Stim_tmp_100 = DM_dat_long %>%  filter(sub_i == n, Contingency==1) %>% 
      pull(stim_go) %>% as.factor() %>% as.numeric()
    DM_dat_long$stim_go_i[(DM_dat_long$sub_i == n) & (DM_dat_long$Contingency == 0.5)] = Stim_tmp_50
    DM_dat_long$stim_go_i[(DM_dat_long$sub_i == n) & (DM_dat_long$Contingency == 1)] = Stim_tmp_100
  }
  
  DM_dat_long = DM_dat_long %>% 
    mutate(contingency_sub_stim_i = paste0(Contingency*100, ",",sub_i,",",stim_go_i),
           contingency_sub_i = paste0(Contingency*100,",",sub_i)) %>%
    left_join(thetas, by = "contingency_sub_stim_i") %>% 
    mutate(theta_scaled = scale(theta)) %>% 
    left_join(thetas_sub, by = "contingency_sub_i")
  
  # # validate which participants' data are good for analysis
  # DM_dat_long %>%
  #   group_by(Contingency, Contingency2, stim_go_i, sub_i, subjectID) %>%
  #   summarise(probe_effect = mean(Outcome, na.rm=T),
  #             theta = mean(theta),
  #             theta_sub = mean(theta_sub)) %>%
  #   group_by(Contingency,sub_i) %>%
  #   mutate(sd_theta = sd(theta)) %>% 
  #   group_by(.) %>%
  #   mutate(theta_sub_scaled = scale(theta_sub),
  #          theta_scaled = scale(theta)) %>% 
  #   filter(theta_sub>0) %>%
  #   ggplot(aes(x = sd_theta, y = theta, color = subjectID)) + 
  #   facet_wrap(. ~ Contingency2) +
  #   theme(legend.position = "none") + 
  #   geom_point()
  
  # Examine Theta-sub and theta-sub-stim association
  DM_dat_full = DM_dat_long %>% group_by(Contingency, Contingency2, stim_go_i, sub_i, subjectID) %>%
    summarise(theta = mean(theta),
              theta_sub = mean(theta_sub))
  theta_cor = cor(DM_dat_full$theta,DM_dat_full$theta_sub)
  cat("\nExp", exp_i, "\n=======\nTheta-correlation: ", round(theta_cor,4), sep = "")
  
  ### Sub-stim and sub Thetas correlation plot (exploratory analysis) ----
  theta_cor_plots[[exp_i]] = DM_dat_full %>% 
    ggplot(aes(x= theta_sub, y = theta)) +
    facet_wrap(. ~ Contingency2, scales = "free") + 
    geom_vline(xintercept = 0.2, linetype = 2, size = 0.25) + 
    geom_smooth(method = "lm", formula = "y ~ x",  color = "black", se = FALSE,  show.legend = FALSE, size = 0.5) + 
    geom_point(aes(fill = Contingency2), color = "black", alpha = 1, size = .75, shape = 21,  stroke = .05) +
    theme_bw() + 
    labs(x = expression(theta[slope[i]]), y = expression(theta[slope[list(i,s)]])) +
    ggtitle("") + 
    theme(text = element_text(size=8), legend.position = "top")
  
  # Filter out participants with low variability in theta values (causes model convergence error)
  sub_2_exclude = DM_dat_long%>%
    group_by(Contingency,sub_i) %>%
    mutate(sd_theta = sd(theta)) %>% 
    filter(sd_theta <= 0.1 | theta_sub < 0) %>% 
    pull(sub_i) %>% unique()
  cat("\nexcluded n = ", length(sub_2_exclude),"\n\n", sep = "")
  
  DM_dat = DM_dat_long %>%
    filter(!(sub_i %in% sub_2_exclude)) %>%
    group_by(Contingency, Contingency2, stim_go_i, sub_i, subjectID) %>%
    summarise(probe_effect = mean(Outcome, na.rm=T),
              chose_go = sum(Outcome==1, na.rm=T),
              chose_nogo = sum(Outcome==0, na.rm=T),
              theta = mean(theta),
              theta_sub = mean(theta_sub)) %>%
    group_by(.) %>%
    mutate(theta_sub_scaled = scale(theta_sub),
           theta_scaled = scale(theta))
  
  ## USE THIS IF SCALE:
  # sd_theta = sd(DM_dat$theta)
  # m_theta = mean(DM_dat$theta)
  ## placeholder if not scaled:
  sd_theta=1
  m_theta=0
  
  # GLMER model (same as modeling each condition separately)
  # For prediction visualization exclude the theta_sub variable
  glmer_model = glmer(cbind(chose_go,chose_nogo) ~ 1 + theta*Contingency2 + 
                        (1 + theta*Contingency2|subjectID), data = DM_dat, family = "binomial")
  # For model comparison use nested model
  glmer_model_theta_sub_and_stim = glmer(cbind(chose_go,chose_nogo) ~ 1 + theta_sub*Contingency2 + theta*Contingency2 + 
                                           (1 + theta*Contingency2|subjectID), data = DM_dat, family = "binomial",
                                         control=glmerControl(optCtrl=list(maxfun=2e6)))
  glmer_model_theta_sub = glmer(cbind(chose_go,chose_nogo) ~ 1 + theta_sub*Contingency2 + 
                                  (1 + Contingency2|subjectID), data = DM_dat, family = "binomial")
  
  theta_probe_model_sub_and_stim[[exp_i]] = glmer_model_theta_sub_and_stim
  theta_probe_model_long_both[[exp_i]] = glmer_model_theta_sub
  theta_probe_model_individual_both[[exp_i]] = glmer_model
  glmer_model_X = model.matrix(glmer_model)
  glmer_model_beta = fixef(glmer_model)
  glmer_model_se = diag(glmer_model_X%*%vcov(glmer_model)%*%t(glmer_model_X))^0.5
  
  # Estimate the parameters - Confidence interval using normal distribution assumption
  X_params = as.matrix(rbind(c(0,1/sd_theta,0,0), # theta 50% contrast
                             c(0,1/sd_theta,0,1/sd_theta), # theta 100% contrast
                             c(0,0,0,1/sd_theta))) # interaction contrast
  params_se = diag(X_params%*%vcov(glmer_model)%*%t(X_params))^0.5
  params_CI = data.frame(exp = rep(exp_i,nrow(X_params)),
                         Parameter = c("theta 50%", "theta 100%", "interaction"),
                         Estimate = X_params%*%glmer_model_beta,
                         SE = params_se) %>%
    mutate(Z = Estimate/SE,
           p = pnorm(abs(Z),lower.tail = FALSE),
           p_two_sided = pnorm(abs(Z),lower.tail = FALSE)*2,
           OR = exp(Estimate),
           OR_CI_lower = exp(Estimate + qnorm(0.025)*SE),
           OR_CI_upper = exp(Estimate + qnorm(0.975)*SE))
  params_CI$p[params_CI$Parameter == "interaction"] = NA
  params_CI$p_two_sided[!(params_CI$Parameter == "interaction")] = NA
  params_individual_CI_both[[exp_i]] = params_CI
  
  # Prediction line for all theta values in range
  theta_50 = DM_dat %>% filter(Contingency==0.5) %>% pull(theta)
  theta_100 = DM_dat %>% filter(Contingency==1) %>% pull(theta)
  # extend the prediction line by a bit (will determine x-axis limits)
  expander_50 = (max(theta_50) - min(theta_50))/50
  expander_100 = (max(theta_100) - min(theta_100))/50
  theta_50_pred = seq(min(theta_50-expander_50),max(theta_50+expander_50), length.out =length(theta_50))
  theta_100_pred = seq(min(theta_100-expander_100),max(theta_100+expander_100), length.out =length(theta_100))
  theta_50_scaled = (theta_50_pred-m_theta)/sd_theta
  theta_100_scaled = (theta_100_pred-m_theta)/sd_theta
  Pred_X = cbind(rep(1,nrow(glmer_model_X)), # intercept
                 c(theta_50_scaled,theta_100_scaled), # theta
                 rep(c(0,1),each = length(theta_50)), # contingency
                 c(rep(0,length(theta_50)), theta_100_scaled)) %>% as.matrix()
  Pred_se = diag(Pred_X%*%vcov(glmer_model)%*%t(Pred_X))^0.5 # prediction variance = X*VCOV*X'
  Pred = DM_dat %>% 
    mutate(
      theta_scaled = Pred_X[,2],
      theta = theta_scaled*sd_theta + m_theta,
      pred = expit(Pred_X%*%glmer_model_beta),
      pred_lower = expit(Pred_X%*%glmer_model_beta + qnorm(0.025)*Pred_se),
      pred_upper = expit(Pred_X%*%glmer_model_beta + qnorm(0.975)*Pred_se)) %>%
    as.data.frame()
  Pred_100 = filter(Pred, Contingency == 1)
  Pred_50 = filter(Pred, Contingency == .5)
  # a data frame which will be used to draw a CI polygon
  Pred_SE = data.frame(Contingency = rep(c(.5,1),each = nrow(DM_dat)),
                       theta =  c(sort(Pred_50$theta),sort(Pred_50$theta,decreasing =TRUE),
                                  sort(Pred_100$theta),sort(Pred_100$theta,decreasing =TRUE)),
                       probe_effect = c(sort(Pred_50$pred_lower), sort(Pred_50$pred_upper,decreasing =TRUE),
                                        sort(Pred_100$pred_lower), sort(Pred_100$pred_upper,decreasing =TRUE))) %>%
    mutate(Contingency2 = factor(Contingency, labels = c("50% Contingency","100% Contingency"))) %>%
    as.data.frame()
  
  # Correlation with choice plot ----
  corplot = ggplot(data = subset(DM_dat), aes(x = theta, y = probe_effect, color = Contingency2)) +
    facet_grid(. ~ Contingency2, scales = "free_x") +
    ggtitle("") +
    theme_bw() + 
    labs(x = expression(theta[slope[list(i,s)]]), y = 'proportion of trials Go stimuli were chosen') +
    geom_hline(yintercept = 0.5, linetype = 2) +
    scale_y_continuous(limits = c(-0.05,1.05), expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    # 95% CI - (min to max range)
    geom_polygon(aes(fill = Contingency2), data = Pred_SE, alpha = 0.3, size=0) +
    geom_line(aes(x = theta, y = pred_lower, group = Contingency2), data = Pred, color = "black", size = .5, linetype =3) + 
    geom_line(aes(x = theta, y = pred_upper, group = Contingency2), data = Pred, color = "black", size = .5, linetype =3) +
    geom_line(aes(x = theta, y = pred, group = Contingency2), data = Pred, color = "black", size = .5) +
    # Points
    set.seed(0) + geom_jitter(aes(fill= Contingency2), height = 0.02, color = "black", alpha = .3, size = .5, shape = 21,  stroke = .05)  + 
    set.seed(0) + geom_jitter(aes(group= Contingency2), height = 0.02, color = "black", alpha = .5, size = .5, shape = 21,  stroke = .05)  + 
    theme(text = element_text(size=8), legend.position = "top")
  
  theta_probe_individual_plot_both[[exp_i]] = corplot
  
  # Correlation with choice plot (same axis) ----
  corplot_same_axis = ggplot(data = subset(DM_dat), aes(x = theta, y = probe_effect, color = Contingency2)) +
    facet_grid(. ~ Contingency2) +
    ggtitle("") +
    theme_bw() + 
    labs(x = expression(theta[slope[list(i,s)]]), y = 'proportion of trials Go stimuli were chosen') +
    geom_hline(yintercept = 0.5, linetype = 2) +
    scale_y_continuous(limits = c(-0.05,1.05), expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    # 95% CI - (min to max range)
    geom_polygon(aes(fill = Contingency2), data = Pred_SE, alpha = 0.3, size=0) +
    geom_line(aes(x = theta, y = pred_lower, group = Contingency2), data = Pred, color = "black", size = .5, linetype =3) + 
    geom_line(aes(x = theta, y = pred_upper, group = Contingency2), data = Pred, color = "black", size = .5, linetype =3) +
    geom_line(aes(x = theta, y = pred, group = Contingency2), data = Pred, color = "black", size = .5) +
    # Points
    set.seed(0) + geom_jitter(aes(fill= Contingency2), height = 0.02, color = "black", alpha = .3, size = .5, shape = 21,  stroke = .05)  + 
    set.seed(0) + geom_jitter(aes(group= Contingency2), height = 0.02, color = "black", alpha = .5, size = .5, shape = 21,  stroke = .05)  + 
    theme(text = element_text(size=8), legend.position = "top")
  
  theta_probe_individual_plot_both_same_axis[[exp_i]] = corplot_same_axis
}

thetas_cor_plot = ggarrange(theta_cor_plots[[1]], theta_cor_plots[[2]], 
          font.label = list(size = 10),
          labels = c("a. Preliminary Exp.", "b. Replication Exp."),
          common.legend = TRUE, ncol = 2, nrow = 1, legend = "none")
thetas_cor_plot
save_plot(thetas_cor_plot, file=paste0(figures_path,"study2_individual_and_sub_theta_correlation"), width = 7, height = 2)

probe_theta_stimuli_corplot = ggarrange(theta_probe_individual_plot_both[[1]], theta_probe_individual_plot_both[[2]], 
          font.label = list(size = 10),
          labels = c("a. Preliminary Exp.", "b. Replication Exp."),
          common.legend = TRUE, ncol = 2, nrow = 1, legend = "none")
probe_theta_stimuli_corplot
save_plot(probe_theta_corplot, file=paste0(figures_path,"study2_theta_probe_individual_stimuli_correlation"), width = 7, height = 2.5)

probe_theta_stimuli_corplot_same_axis = ggarrange(theta_probe_individual_plot_both_same_axis[[1]], theta_probe_individual_plot_both_same_axis[[2]], 
                                        font.label = list(size = 10),
                                        labels = c("a. Preliminary Exp.", "b. Replication Exp."),
                                        common.legend = TRUE, ncol = 2, nrow = 1, legend = "none")
probe_theta_stimuli_corplot_same_axis
save_plot(probe_theta_stimuli_corplot_same_axis, file=paste0(figures_path,"study2_theta_probe_individual_stimuli_correlation_same_axis"), width = 7, height = 2.5)

### Probe correlation summary table ----
params_individual_CI_both[[1]] %>% 
  rbind(params_individual_CI_both[[2]]) %>%
  mutate(p  = as.character(signif(p, 2)),
         p_two_sided  = as.character(signif(p_two_sided, 2))) %>%
  kbl(digits = 2) %>% kable_minimal()

### Compare model with individual theta to subject-level theta ----
for (exp_i in 1:2){
  cat("\nExperiment ",exp_i,"\n==========\nAIC difference:",
      AIC(theta_probe_model_long_both[[exp_i]]) - AIC(theta_probe_model_sub_and_stim[[exp_i]]),"\n\n",
      sep = "")
  anova(theta_probe_model_sub_and_stim[[exp_i]],theta_probe_model_long_both[[exp_i]], test = "LRT") %>% 
    print()
}
