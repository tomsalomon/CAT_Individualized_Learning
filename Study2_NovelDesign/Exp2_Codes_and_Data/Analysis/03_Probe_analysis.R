
library(lme4)
library("rstudioapi")    
library("ggplot2")    
library(reshape2)
library(dplyr)
rm(list=ls())

# Valid participants (n = 53):
subjects = c(101:105,107:114,117:118,121:127,131:154,156:162)

# Excluded participants:
# 115, 116, 130 - Transitivity
# 128, 155 - FA (max = 10.8%)
# 106, 119, 120, 155 - Miss (max = 28.12%)
# 129 - Technical issues

# Assisting functions ----
# Standard error function
se = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }

CI = function (in_model, var_num = '',one_sided = F) {
  # Confidence interval from logistic regression model / model summary function.
  # in_model = input model; [var_num] = variable in model to analyze (optional)
  
  # if the input is a model (not model summary), convert to summary
  if (class(in_model) == "glmerMod"){
    in_model = summary(in_model)
  }
  coef = as.data.frame(in_model$coefficients)
  coef$p = coef$`Pr(>|z|)`
  ci_output = c()
  ci_logs = as.data.frame(cbind(coef[,1],coef[,1]-1.96*coef[,2],coef[,1]+1.96*coef[,2]))
  ci_OR = exp(ci_logs)
  colnames(ci_logs) = c("logOR","logOR_CI_lower","logOR_CI_upper")
  colnames(ci_OR) = c("OR","OR_CI_lower","OR_CI_upper")
  ci_output = cbind(coef,ci_OR,ci_logs)
  
  if (is.numeric(var_num)) {
    ci_output = ci_output[var_num,]
  }
  if (one_sided==TRUE){
    ci_output$p = ci_output$p/2
    ci_output$`Pr(>|z|)` = ci_output$`Pr(>|z|)`/2
  }
  return (ci_output)
}

# Get current path
script_path = rstudioapi::getActiveDocumentContext()$path
pathSplit=strsplit(script_path, "/")
pathSplit=pathSplit[[1]]
main_path=paste0(pathSplit[1:(length(pathSplit)-2)],"/",collapse="")
path=paste0(main_path,"/Output/")

# load data
filelist=c()
for (s in subjects){
  filelist=c(filelist,Sys.glob(paste(path, "BM_",s,"_Probe_*.txt",sep="")))
}
Probe_data=c()
for (f in filelist){
  Probe_data=rbind(Probe_data,read.table(f,header=T,na.strings=c(999,999000)))
}

Probe_data$PairType2 = as.factor(Probe_data$PairType)
Probe_data$Contingency2 = factor(Probe_data$Contingency, labels = c("50%", "100%"))

# Descriptive Statistics ----

n = length(unique(Probe_data$subjectID))
writeLines(paste("n =",n)) # Number of Participants
writeLines ("Proportion of trials Go items were chosen over NoGo:")
tapply(Probe_data$Outcome,list(Probe_data$Contingency2),mean,na.rm=T)
writeLines ("Proportion of trials Go items were chosen over NoGo (By value groups):")
tapply(Probe_data$Outcome,list(Probe_data$PairType,Probe_data$Contingency2),mean,na.rm=T)

# Logistic Regression ----
Contingency50 = summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(Probe_data,(Probe_data$Contingency==0.5)),na.action=na.omit,family=binomial)) 
print(Contingency50)

Contingency100 = summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(Probe_data,(Probe_data$Contingency==1)),na.action=na.omit,family=binomial)) 
print(Contingency100)

ContingencDiff = summary(glmer(Outcome ~ 1 + Contingency2 + (1 + Contingency2 |subjectID),data=(Probe_data),na.action=na.omit,family=binomial)) 
print(ContingencDiff)

Contingency50_CI = CI(Contingency50, one_sided = T)
Contingency100_CI = CI(Contingency100, one_sided = T)
ContingencyDiff_CI = CI(ContingencDiff, var_num = 2)

# Summary Dataframe ----
# Mean proportions by subject
Probe_data_sub = as.data.frame(Probe_data %>% group_by(subjectID, Contingency=Contingency2) %>% summarise(value=(mean(Outcome,na.rm=T))))
Probe_data_sub$ContingencyNum = as.numeric(Probe_data_sub$Contingency)
Probe_data_sub$value_diff = rep(Probe_data_sub$value[Probe_data_sub$Contingency=="100%"]-Probe_data_sub$value[Probe_data_sub$Contingency=="50%"],each=2)

# Group summary statistics
SummaryStats = as.data.frame(Probe_data_sub %>% group_by(Contingency) %>% summarise(Mean = mean(value), SE = se(value)))
SummaryStats$ymin = SummaryStats$Mean - SummaryStats$SE
SummaryStats$ymax = SummaryStats$Mean + SummaryStats$SE
SummaryStats$p = c(Contingency50_CI$p,Contingency100_CI$p)
SummaryStats$asterisk=ifelse(SummaryStats$p<.001, "***", ifelse(SummaryStats$p<.01, "**", ifelse(SummaryStats$p<.05, "*", "")))

# Bar Plot ----
bar_plot = ggplot() +
  geom_bar(data = SummaryStats, aes(x = Contingency, y = Mean, fill = Contingency), stat = 'identity', color = "black", width=.5) + 
  set.seed(0) + geom_jitter(data = Probe_data_sub, aes(x = Contingency, y = value, group = subjectID), color = "black", size = 2, width = 0.1, alpha= 0.4) +
  set.seed(0) + geom_jitter(data = Probe_data_sub, aes(x = Contingency, y = value, group = subjectID, color= Contingency), size = 1, width = 0.1,  alpha= 0.4) +
  scale_y_continuous(limits = c(0,1), breaks=seq(0, 1, 0.1),expand = c(0,0),  name = "Proportion of trials Go stimuli were chosen") +
  xlab("Go training contingency") +
  theme_bw() +
  theme(legend.position="none") + 
  geom_text(data = SummaryStats, aes(x=Contingency, y = ymax+0.05, label = asterisk), size = 5) +
  geom_errorbar(data = SummaryStats, aes(ymin=ymin, ymax=ymax,x = Contingency), width=.1) + 
  geom_abline(intercept = (0.5),slope=0,linetype =2) # chace level 50% reference line
if (ContingencyDiff_CI$p<0.05) {
  h = max(SummaryStats$ymax) + 0.1
  asterisk = ifelse(ContingencyDiff_CI$p < 0.001, "***", ifelse(ContingencyDiff_CI$p < 0.01,"**","*"))
  bar_plot = bar_plot  + geom_path(aes(x = c(1,1,2,2), y = c(h, h+0.02,h+0.02,h))) + 
    geom_text(aes(x=1.5, y = h+0.05, label = asterisk), size = 5)
}
bar_plot

# Bar Plot2 ----
bar_plot = ggplot() +
  geom_bar(data = SummaryStats, aes(x = Contingency, y = Mean, fill = Contingency), stat = 'identity', color = "black", width=.5) + 
  set.seed(0) + geom_jitter(data = Probe_data_sub, shape=21 , aes(x = Contingency, y = value, fill= Contingency), size = 2, width = 0.1,  alpha= .4) +
  set.seed(0) + geom_jitter(data = Probe_data_sub, shape=21 , aes(x = Contingency, y = value, group= Contingency), size = 2, width = 0.1,  alpha= 1) +
  scale_y_continuous(limits = c(0,1), breaks=seq(0, 1, 0.1),expand = c(0,0),  name = "Proportion of trials Go stimuli were chosen") +
  xlab("Go training contingency") +
  theme_bw() +
  theme(legend.position="none") + 
  geom_text(data = SummaryStats, aes(x=Contingency, y = ymax+0.05, label = asterisk), size = 5) +
  geom_errorbar(data = SummaryStats, aes(ymin=ymin, ymax=ymax,x = Contingency), width=.1) + 
  geom_abline(intercept = (0.5),slope=0,linetype =2) # chace level 50% reference line
if (ContingencyDiff_CI$p<0.05) {
  h = max(SummaryStats$ymax) + 0.1
  asterisk = ifelse(ContingencyDiff_CI$p < 0.001, "***", ifelse(ContingencyDiff_CI$p < 0.01,"**","*"))
  bar_plot = bar_plot  + geom_path(aes(x = c(1,1,2,2), y = c(h, h+0.02,h+0.02,h))) + 
    geom_text(aes(x=1.5, y = h+0.05, label = asterisk), size = 5)
}
bar_plot

# Bar Plot3 ----
bar_plot = ggplot() +
  geom_bar(data = SummaryStats, aes(x = Contingency, y = Mean, fill = Contingency), stat = 'identity', color = "black", width=.5) + 
  set.seed(0) + geom_point(data = Probe_data_sub, shape=21 , aes(x = Contingency, y = value, group=Contingency, fill= Contingency), size = 2, alpha= 1, position = position_jitterdodge(jitter.width = .5)) +
  scale_y_continuous(limits = c(0,1), breaks=seq(0, 1, 0.1),expand = c(0,0),  name = "Proportion of trials Go stimuli were chosen") +
  xlab("Go training contingency") +
  theme_bw() +
  theme(legend.position="none") + 
  geom_text(data = SummaryStats, aes(x=Contingency, y = ymax+0.05, label = asterisk), size = 5) +
  geom_errorbar(data = SummaryStats, aes(ymin=ymin, ymax=ymax,x = Contingency), width=.1) + 
  geom_abline(intercept = (0.5),slope=0,linetype =2) # chace level 50% reference line
if (ContingencyDiff_CI$p<0.05) {
  h = max(SummaryStats$ymax) + 0.1
  asterisk = ifelse(ContingencyDiff_CI$p < 0.001, "***", ifelse(ContingencyDiff_CI$p < 0.01,"**","*"))
  bar_plot = bar_plot  + geom_path(aes(x = c(1,1,2,2), y = c(h, h+0.02,h+0.02,h))) + 
    geom_text(aes(x=1.5, y = h+0.05, label = asterisk), size = 5)
}
bar_plot
    
# Box Plot ----
ggplot(data = Probe_data_sub, aes(x = Contingency, y = value, fill = Contingency)) +
  geom_boxplot(data = Probe_data_sub, aes(x = Contingency, y = value, group = Contingency)) + 
  set.seed(0) + geom_jitter(size = 2, width = 0.1, alpha = .7, color = 'black') + 
  set.seed(0) + geom_jitter(size = 1, width = 0.1, alpha = .5, aes(color = Contingency)) + 
  theme_bw() + 
  geom_abline(intercept = (0.5),slope=0,linetype =2) + # chace level 50% reference line
  theme(legend.position="none") + 
  xlab("Go training contingency") +
  scale_y_continuous(limits = c(0,1), breaks=seq(0, 1, 0.1),expand = c(0,0),  name = "Proportion of trials Go stimuli were chosen") 

# Box Plot: colored by slope ----
ggplot(data = Probe_data_sub, aes(x = ContingencyNum+scale(value)/15, y = value, group = subjectID,color = value_diff)) +
  geom_boxplot(data = Probe_data_sub, aes(x = Contingency, y = value, group = Contingency), fill = "white") + 
  geom_line(size=.7, color="black",alpha = .9) +
  geom_line(size=.5, alpha = .9) +
  geom_point(size = 3, color="black",alpha = .9) + 
  geom_point(size = 2.5, alpha = .9) +
  theme_bw()+
  xlab("Go training contingency") +
  scale_color_gradient2(low = "firebrick2", mid = "white", high = "green4", midpoint = 0) +
  geom_abline(intercept = (0.5),slope=0,linetype =2) + # chace level 50% reference line
  scale_y_continuous(limits = c(0,1),  breaks=seq(0, 1, 0.1),expand = c(0,0), name = "Proportion of trials Go stimuli were chosen")


# Summary Dataframe: By value categories  ----
# Mean proportions by subject and value category
Probe_data_sub_value= dcast(data = Probe_data, formula =  subjectID + Contingency2 + PairType2 ~ ., fun.aggregate = mean, value.var = "Outcome", na.rm = T)
colnames(Probe_data_sub_value)[ncol(Probe_data_sub_value)] = "Mean"
Probe_data_sub_value$ContingencyNum = as.numeric(Probe_data_sub_value$Contingency2)
Probe_data_sub_value$PairType = as.numeric(Probe_data_sub_value$PairType2)

# Group summary statistics
SummaryStatsValue = dcast(data = Probe_data_sub_value, formula = Contingency2 + PairType2 ~ ., fun.aggregate = mean, value.var = "Mean", na.rm = T)
colnames(SummaryStatsValue)[ncol(SummaryStatsValue)] = "Mean"

SummaryStatsValueSE = dcast(data = Probe_data_sub_value, formula = Contingency2 + PairType2 ~ ., fun.aggregate = se, value.var = "Mean")
colnames(SummaryStatsValueSE)[ncol(SummaryStatsValueSE)] = "SE"
SummaryStatsValue = merge(SummaryStatsValue,SummaryStatsValueSE)

SummaryStatsValue$ymin = SummaryStatsValue$Mean - SummaryStatsValue$SE
SummaryStatsValue$ymax = SummaryStatsValue$Mean + SummaryStatsValue$SE

# Bar Plot: By value categories ----
ggplot() +
  geom_bar(data = SummaryStatsValue, aes(x = PairType2, y = Mean, fill = Contingency2), stat = 'identity', color = "black", width=.8, position = "dodge") + 
  scale_y_continuous(limits = c(0,1.01), breaks=seq(0, 1, 0.1),expand = c(0,0),  name = "Proportion of trials Go stimuli were chosen") +
  xlab("Value category") +
  theme_bw() +
  theme(legend.position="top") + 
  geom_errorbar(data = SummaryStatsValue, aes(ymin=ymin, ymax=ymax,x = PairType2, group = Contingency2), position = position_dodge(width=0.8), width = 0.4) + 
  geom_abline(intercept = (0.5),slope=0,linetype =2)  # chance level 50% reference line
  # geom_text(data = SummaryStatsValue, aes(x=PairType2, y = ymax+0.05, label = asterisk), size = 5) +
  # set.seed(0) + geom_jitter(data = subset(Probe_data_sub_value,ContingencyNum ==2), aes(x = PairType+0.2, y = Mean),color = "black", size = 1.5, width = 0.1,  alpha= 0.7) +
  # set.seed(0) + geom_jitter(data = subset(Probe_data_sub_value,ContingencyNum ==2), aes(x = PairType+0.2, y = Mean, color= Contingency2), size = 1, width = 0.1,  alpha= 0.4) + 
  # set.seed(0) + geom_jitter(data = subset(Probe_data_sub_value,ContingencyNum ==1), aes(x = PairType-0.2, y = Mean),color = "black", size = 1.5, width = 0.1,  alpha= 0.7) +  
  # set.seed(0) + geom_jitter(data = subset(Probe_data_sub_value,ContingencyNum ==1), aes(x = PairType-0.2, y = Mean, color= Contingency2), size = 1, width = 0.1,  alpha= 0.4) 
  
# Exploratory analyses: test value effect ----
m1 = glmer(Outcome ~ 1 + (1|subjectID), family = binomial, data = Probe_data)
m2 = glmer(Outcome ~ 1 + Contingency2 + (1|subjectID), family = binomial, data = Probe_data)
m3 = glmer(Outcome ~ 1 + Contingency2 + PairType2 + (1|subjectID), family = binomial, data = Probe_data)
m4 = glmer(Outcome ~ 1 + Contingency2 * PairType2 + (1|subjectID), family = binomial, data = Probe_data)

anova(m2,m3)
anova(m3,m4)

m1 = (glmer(Outcome ~ 1 + (1|subjectID), family = binomial, data = subset(Probe_data, Contingency2 == "100%")))
m2 = (glmer(Outcome ~ 1 + PairType + (1 + PairType|subjectID), family = binomial, data = subset(Probe_data, Contingency2 == "100%")))
summary(m2)
anova(m1,m2)


m1 = (glmer(Outcome ~ 1 + (1|subjectID), family = binomial, data = subset(Probe_data, Contingency2 == "50%")))
m2 = (glmer(Outcome ~ 1 + PairType + (1 + PairType|subjectID), family = binomial, data = subset(Probe_data, Contingency2 == "50%")))
summary(m2)
anova(m1,m2)





