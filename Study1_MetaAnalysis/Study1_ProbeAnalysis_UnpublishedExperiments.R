### Initialization ----

# Load rquired R packages
library(rstudioapi) # to define current working directory
library(ggplot2)
library(lme4)
library(lmerTest)
library(lm.beta)
library(dplyr)
#library(kableExtra)
# library(knitr)

# clear workspace
rm(list=ls())

# Define the current script location as the working directory
pwd = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(pwd)
figures_path = paste0(pwd,'/Figures/')

# Standard error function
se = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }
expit = function(u){exp(u)/(1+exp(u))}

# Load data from appropriate dir location (change if not saved in the same directory as the plotting script)
data_path = pwd # Change according to the path where the data is saved
load(paste0(data_path,"/probe_data.Rda"))
load(paste0(data_path,"/training_data.Rda"))

# count_unique function
count_unique = function(x) {length(unique(x))}

# training summary
training_summary = training_data %>% group_by(Experiment) %>% summarise(n_runs = max(runNum))
# count number of subjects per experiment
probe_data %>% group_by(Experiment, ExperimentName) %>% summarise(N = count_unique(subjectID)) %>% print(n=100)

# Summarize probe results ----
probe_data$unpublished = grepl(pattern = 'UnPub', x= probe_data$ExperimentName)
# Improve long experiment titles
probe_data$ExperimentName2 = gsub(pattern = 'UnPub',replacement = "\n", x = probe_data$ExperimentName) 
probe_data$ExperimentName2 = gsub(pattern = '_',replacement = " ", x = probe_data$ExperimentName2) 
probe_data$ExperimentName2 = gsub(pattern = 'p0',replacement = "p. ", x = probe_data$ExperimentName2) 
probe_data$ExperimentName2 = gsub(pattern = 'memory experiment',replacement = "mem.", x = probe_data$ExperimentName2) 
probe_data$ExperimentName2 = gsub(pattern = '12 training run Israel',replacement = "pilot I", x = probe_data$ExperimentName2) 
probe_data$ExperimentName2 = gsub(pattern = '20 training run Israel',replacement = "pilot II", x = probe_data$ExperimentName2) 

probe_summary_sub = probe_data %>% filter(PairType<=2) %>% group_by(Experiment, ExperimentName2,subjectID, PairType) %>% 
  summarise(propChoseGo_sub = mean(Outcome, na.rm=T), unpublished = mean(unpublished)) 

probe_summary =  probe_summary_sub %>% group_by(Experiment, ExperimentName2, PairType) %>% summarise(N = n(), unpublished = mean(unpublished), propChoseGo = mean(propChoseGo_sub), SEM = se(propChoseGo_sub)) %>%
  mutate(LogOdds = NA, LogOdds_e = NA, z = NA, p = NA) %>% as.data.frame()

for (i in 1:nrow(probe_summary)) {
  exp_i = probe_summary$Experiment[i]
  PairType_i = probe_summary$PairType[i]
  dat = probe_data %>% filter(Experiment==exp_i,PairType == PairType_i, Outcome <=1)
  model_i = glmer(Outcome ~ 1 + (1|subjectID2), data = dat, family = binomial)
  model_summary_i = summary(model_i)
  probe_summary[i,c("LogOdds","LogOdds_e","z","p")] = model_summary_i$coefficients 
  cat("\rProgress: ", 100*round(i/nrow(probe_summary),2),"%",
      sep="")
  flush.console() 
}
# summary statistics
probe_summary = left_join(probe_summary,training_summary,by = "Experiment")
probe_summary$p_one_sided = probe_summary$p/2
probe_summary$Odds = exp(probe_summary$LogOdds)
probe_summary$Odds_CI_lower = exp(probe_summary$LogOdds + probe_summary$LogOdds_e*qnorm(p=0.025))
probe_summary$Odds_CI_upper = exp(probe_summary$LogOdds + probe_summary$LogOdds_e*qnorm(p=0.975))
probe_summary$prop = probe_summary$Odds/(probe_summary$Odds+1)
probe_summary$prop_CI_lower = probe_summary$Odds_CI_lower/(probe_summary$Odds_CI_lower+1)
probe_summary$prop_CI_upper = probe_summary$Odds_CI_upper/(probe_summary$Odds_CI_upper+1)
probe_summary$asterisk <- cut(x = probe_summary$p_one_sided, breaks = c(0,.001,.01,.05,.1,1), labels = c("***","**","*","",""))
probe_summary$PairType2 = factor(probe_summary$PairType, labels = c("High Value", "Low Value"))

probe_summary_unpblished = probe_summary %>% filter(unpublished==1)
probe_unpblished_sub = probe_summary_sub %>% filter(unpublished==1)
probe_summary_unpblished$ExperimentName2 = as.factor(probe_summary_unpblished$ExperimentName2)
probe_unpblished_sub$ExperimentName2 = as.factor(probe_unpblished_sub$ExperimentName2)
probe_unpblished_sub$PairType2 = factor(probe_unpblished_sub$PairType, labels = c("High Value", "Low Value"))

# Plot probe results ----
ggplot(data = probe_summary_unpblished, aes(x = ExperimentName2, y = prop, fill = PairType2)) +
  geom_bar(stat = 'identity', position = 'dodge', width = 2/3, color = "black") + 
  geom_point(data = probe_unpblished_sub, aes(x = ExperimentName2,y = propChoseGo_sub, group = PairType2),
             position = position_jitterdodge(seed = 0), stroke = 0.2, size = 1, alpha  = .4, shape = 21, show.legend = FALSE, color = "black") +
  geom_hline(yintercept = 0.5, linetype = 2) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0,0,0.0), breaks = seq(0,1,0.1)) + 
  geom_errorbar(aes(ymin = prop_CI_lower, ymax = prop_CI_upper), width = 1/3, position = position_dodge(width = 2/3)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(y = "Proportion of trials Go stimuli were chosen" , x = "Experiment") + 
  theme_bw() +
  scale_fill_discrete(guide = guide_legend(reverse = FALSE) ) +
  #coord_flip() + 
  # scale_x_discrete(limits = rev(levels(probe_summary_unpblished$ExperimentName2))) + 
  theme(axis.title.x = element_blank(), legend.title = element_blank(), 
        legend.position = c(.01,.99), legend.background = element_rect(fill = NA), legend.direction = "horizontal",
        legend.justification = c("left","top")) +
  geom_text(aes(label = asterisk, y = prop_CI_upper + 0.05),position = position_dodge(width = 2/3), size = 6)
dev.copy2pdf(file=paste0(figures_path,"MetaAnalysis_UnpublishedData_ProbeResults.pdf"), width = 7, height = 3.5)

probe_table = probe_summary_unpblished %>% 
  mutate(OR_text = sprintf("%.2f [%.2f, %.2f]",Odds,Odds_CI_lower, Odds_CI_upper),
         propChoseGo_text =  sprintf("%s %.1f%%",PairType2, propChoseGo*100),
         p_text = sprintf("%.3f",p_one_sided)) 
probe_table[probe_table$p_one_sided<0.001,"p_text"] = format(probe_table[probe_table$p_one_sided<0.001,"p_one_sided"], digits = 2)
rows2modify = seq(from = 2,to = nrow(probe_table),by = 2)
probe_table[rows2modify, c("N","n_runs")] = ""
probe_table$ExperimentName = gsub(pattern = "\n",replacement = "", x = probe_table$ExperimentName2) 
probe_table$ExperimentNamePub = NA
probe_table$ExperimentNamePub[rows2modify-1] = substr(probe_table$ExperimentName[rows2modify-1],0,7)
probe_table$ExperimentNamePub[rows2modify] = substr(probe_table$ExperimentName[rows2modify],8, nchar(probe_table$ExperimentName[rows2modify]))


# summary table to save and for publication
probe_table_short = probe_table %>% select(ExperimentNamePub,  N, n_runs ,propChoseGo_text, p_text, OR_text) 
colnames(probe_table_short) = c("Experiment", "N", "Training runs", "Prop. Go were chosen", "p (one-sided)", "OR (95% CI)")

# Save table and print it to word using the "R_plots_to_word.Rmd" R markdown script
# save(probe_table_short,file = "probe_table_short.rda")
# rmarkdown::render("R_plots_to_word.Rmd","word_document")




