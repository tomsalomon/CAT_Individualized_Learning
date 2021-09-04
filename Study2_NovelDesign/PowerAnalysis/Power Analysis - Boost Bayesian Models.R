############################################################### -
# Power Analysis for Correlation coefficient using Bootstrap ----
############################################################### -
# By: Tom Salomon, January 2020
# This power analysis was used to calculate the minimal n required for the preregistered experiment at: https://osf.io/yswu7/
# The script will analyze the data from a pilot study with n = 20 valid participants. 
#
# Instruction: ----
# Install R and RStudio
# Install the required R packages(lines 13-18)
# Save the pilot data (unzip it first) and 'fit_stan2.Rda' (fitted Bayesian model results) and 'Power_analysis_Results.rda' (pre-run power analysis)
# on you local computer in the same directory as this script.
# In line 66, change the variable 'Skip_bootstrapping = TRUE' if you want to run the power analysis again (takes about 2-3 hours)

# Required packages
library(rstan)
library(lme4)
library(rstudioapi)    
library(ggplot2)    
library(reshape2)

rm(list=ls()) # clear workspace
start_time = Sys.time()
set.seed(1) # randomization seed for reproducibility

# Get current path
script_path = rstudioapi::getActiveDocumentContext()$path
pathSplit=strsplit(script_path, "/")
pathSplit=pathSplit[[1]]
current_path=paste0(pathSplit[1:(length(pathSplit)-1)],"/",collapse="")
setwd(current_path) # set working directory as the path where the script is saved
path_output = paste0(current_path,"PilotData/") # path of output data

# 1. Load probe and Bayesian model data ----

# Load Bayesian model data 
load(file = paste0(current_path,"fit_stan2.Rda")) # load the fitted Stan model
model_summary = as.data.frame(summary(fit_stan2)$summary)
# Load Theta slope parameter estimates for 100% and 50% contingency
theta1_mean = model_summary$mean[grep("theta_b1_rand",row.names(model_summary))] # 100% contingency theta slope
theta2_mean = model_summary$mean[grep("theta_b2_rand",row.names(model_summary))] # 50% contingency theta slope

# Load probe data 
subjects=c(101:110,112:121); # valid subjects
# Not really excluded:
# 111 - Eyetracker crashed the code during the experiment. 

filelist=c()
for (s in subjects){
  filelist=c(filelist,Sys.glob(paste(path_output, "BM_",s,"_Probe_*.txt",sep="")))
}
Probe_data=c()
for (f in filelist){
  Probe_data=rbind(Probe_data,read.table(f,header=T,na.strings=c(999,999000)))
}

Probe_data$sub_i2 = as.numeric(factor(Probe_data$subjectID))
probe_mean = tapply(Probe_data$Outcome, list(Probe_data$sub_i2, Probe_data$Contingency), mean, na.rm = T)
DM_dat = melt(probe_mean, varnames = c("sub","Contingency"), value.name = "probe_effect")
DM_dat$theta = c(theta2_mean,theta1_mean)
DM_dat$Contingency2 = as.factor(DM_dat$Contingency)
levels(DM_dat$Contingency2) = c("50%", "100%")
DM_dat50 = subset(DM_dat,Contingency==.5 ) 
DM_dat100 = subset(DM_dat,Contingency==1 ) 

# 2. Bootstrap -----
Skip_bootstrapping = TRUE # Change to 'TRUE' to skip the analysis, and load previous analysis results. Change to 'FALSE to run the full analysis (about 3 hours)
if (Skip_bootstrapping){
  # Load simulations workspace .rda file:
  load(file = "Power_analysis_Results.rda") # uncomment to load workspace
} else {
  iter = 1000  #num of simulations per sample size
  Sample_size_vector = c(seq(5,30,5),seq(32,50,2),seq(51,75,1)) #sample sizes: 5,10,15, ... , 30, 32, 34, ..., 50, 51, 52, ... , 74, 75.
  alpha = 0.05
  power = 0.95
  
  ## BOTTSTRAP ##
  Significant_Prop_Vector50 = length(Sample_size_vector)
  Significant_Prop_Vector100 = length(Sample_size_vector)
  Significant_Prop_Vector_log_reg50 = length(Sample_size_vector)
  Significant_Prop_Vector_log_reg100 = length(Sample_size_vector)
  PvalueVector_log_reg50 = matrix(0,nrow = iter,ncol = length(Sample_size_vector))
  PvalueVector_log_reg100 = matrix(0,nrow = iter,ncol = length(Sample_size_vector))
  PvalueVector_corr50 = matrix(0,nrow = iter,ncol = length(Sample_size_vector))
  PvalueVector_corr100 = matrix(0,nrow = iter,ncol = length(Sample_size_vector))
  
  for (K in 1:length(Sample_size_vector)){
    sample_size = Sample_size_vector[K]
    SampledTteta50 = numeric(sample_size)
    SampledTteta100 = numeric(sample_size)
    SmapledProbe50 = numeric(sample_size)
    SmapledProbe100 = numeric(sample_size)
    
    
    pb <- txtProgressBar()
    
    for (J in 1:iter){
      # Sample
      Sample = sample(seq(from=1,to=20,by=1), size=sample_size, replace = TRUE)
      probe_dat_itr = c()
      
      for (I in 1:sample_size){ 
        SampledTteta50[I]= DM_dat50$theta[Sample[I]]
        SmapledProbe50[I]= DM_dat50$probe_effect[Sample[I]]
        SampledTteta100[I]= DM_dat100$theta[Sample[I]]
        SmapledProbe100[I]= DM_dat100$probe_effect[Sample[I]]
        probe_dat_tmp = subset(Probe_data,sub_i2==Sample[I])
        probe_dat_tmp$sub_i2=I
        probe_dat_itr = rbind(probe_dat_itr,probe_dat_tmp)
      }
      
      Cor_Test50 =  cor.test(SampledTteta50,SmapledProbe50 ,method = "pearson")
      Cor_Test100 =  cor.test(SampledTteta100,SmapledProbe100 ,method = "pearson")
      PvalueVector_corr50[J,K] = Cor_Test50$p.value
      PvalueVector_corr100[J,K] = Cor_Test100$p.value
      
      probe_dat_itr$sub_i2 = as.factor(probe_dat_itr$sub_i2)
      log_reg_test50 = summary(suppressMessages(glmer(Outcome ~ 1 + (1|sub_i2), data  = subset(probe_dat_itr,Contingency==0.5),family = binomial)))
      log_reg_test100 = summary(suppressMessages(glmer(Outcome ~ 1 + (1|sub_i2), data  = subset(probe_dat_itr,Contingency==1),family = binomial)))
      PvalueVector_log_reg50[J,K]= log_reg_test50$coefficients[1,4]
      PvalueVector_log_reg100[J,K]= log_reg_test100$coefficients[1,4]
      
      setTxtProgressBar(pb, J/iter)
    }
    
    Significant_Prop_Vector50[K] = mean(PvalueVector_corr50[,K] < alpha)
    Significant_Prop_Vector100[K] = mean(PvalueVector_corr100[,K] < alpha)
    Significant_Prop_Vector_log_reg50[K] = mean(PvalueVector_log_reg50[,K] < alpha)
    Significant_Prop_Vector_log_reg100[K] =  mean(PvalueVector_log_reg100[,K] < alpha)
    
    if (min(c(Significant_Prop_Vector50[K],
              Significant_Prop_Vector100[K],
              Significant_Prop_Vector_log_reg50[K],
              Significant_Prop_Vector_log_reg100[K])) >= power){
      countdown2stop = countdown2stop-1
    } else {
      countdown2stop = 5 # stop sampling after 5 consequtive well powered N's
    }
    # Progress
    passed_time = difftime(Sys.time(),start_time, units = "mins")
    progress = K/length(Sample_size_vector)
    print("")
    print(paste0("Completed bootstrapping, sample size: n = ",sample_size ,", progress = ",round(100*progress,2), "%"))
    print(paste0("Estimated time until analysis is completed: ",round((1-progress)*passed_time/progress,2), " minutes"))
    if (countdown2stop==0){
      Significant_Prop_Vector_log_reg100 = Significant_Prop_Vector_log_reg100[1:K]
      Significant_Prop_Vector_log_reg50 = Significant_Prop_Vector_log_reg50[1:K]
      Significant_Prop_Vector100 = Significant_Prop_Vector100[1:K]
      Significant_Prop_Vector50 = Significant_Prop_Vector50[1:K]
      Sample_size_vector = Sample_size_vector[1:K]
      break
    }
  }
  # Save workspace
  save.image(file = "Power_analysis_Results.rda") # uncomment to save workspace
} # End of "Skip Bootstrapping"

# 3. Plot power analysis results ----
power_results = data.frame(N = rep(Sample_size_vector,4), 
                           power = c(Significant_Prop_Vector_log_reg100,
                                     Significant_Prop_Vector_log_reg50,
                                     Significant_Prop_Vector100,
                                     Significant_Prop_Vector50),
                           analysis_i = sort(rep(1:4,K)))
power_results$analysis = factor(power_results$analysis_i)
levels(power_results$analysis) = c("Probe logistic regression\n(100% Contingency)\n",
                                   "Probe logistic regression\n(50% Contingency)\n",
                                   "Theta-probe correlation\n(100% Contingency)\n",
                                   "Theta-probe correlation\n(50% Contingency)\n")
selected_N = Sample_size_vector[K-4]

p = ggplot(data = subset(power_results),aes(x=N, y = power, color = analysis)) + 
  theme_bw() +
  geom_hline(yintercept = .80, linetype=3) +
  geom_hline(yintercept = .95, linetype=3) +
  scale_x_continuous(breaks = seq(0,Sample_size_vector[K],5), expand = c(0.01,0.01),name= "n") + 
  scale_y_continuous(breaks = seq(0,1,.1), limits = c(0,1), name= expression(paste("power (1 - ",beta,")"))) + 
  geom_vline(xintercept = selected_N, linetype = 2) +
  annotate(geom="text", x=selected_N-4, y=.7, label=paste0("n = ",selected_N), size =5) +
  geom_point() +
  geom_line() +
  theme(text = element_text(size = 16))
print(p)

# Save plot
pdf(file = "Power_analysis_Results.pdf",width = 9, height = 4)
plot(p)
dev.off()
