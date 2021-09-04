
library(lme4)
library("rstudioapi")    
library("ggplot2")    
library(reshape2)
rm(list=ls())

## Sample
subjects=c(101:110,112:121); # 20 valid subjects
# Not really excluded:
# 111 - Eyetracker crashed code during the experiment. 

# Standard error function
se = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }

CI = function (in_model, var_num = '') {
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
  return (ci_output)
}

# Get current path
script_path = rstudioapi::getActiveDocumentContext()$path
pathSplit=strsplit(script_path, "/")
pathSplit=pathSplit[[1]]
main_path=paste0(pathSplit[1:(length(pathSplit)-2)],"/",collapse="")
path=paste0(main_path,"/Output/")

filelist=c()
for (s in subjects){
  filelist=c(filelist,Sys.glob(paste(path, "BM_",s,"_Probe_*.txt",sep="")))
}

Probe_data=c()
for (f in filelist){
  Probe_data=rbind(Probe_data,read.table(f,header=T,na.strings=c(999,999000)))
}

Probe_data$PairType2 = as.factor(Probe_data$PairType)
Probe_data$Contingency2 = as.factor(Probe_data$Contingency)
levels(Probe_data$Contingency2) = c("50%", "100%")
# Descriptive Statistics ----

n = length(unique(Probe_data$subjectID))
writeLines(paste("n =",n)) # Number of Participants
writeLines ("Proportion of trials Go items were chosen over NoGo:")
tapply(Probe_data$Outcome,list(Probe_data$Contingency2),mean,na.rm=T)
writeLines ("Proportion of trials Go items were chosen over NoGo (By subject):")
tapply(Probe_data$Outcome,list(Probe_data$subjectID,Probe_data$Contingency2),mean,na.rm=T)
writeLines ("Proportion of trials Go items were chosen over NoGo (By value groups):")
tapply(Probe_data$Outcome,list(Probe_data$PairType,Probe_data$Contingency2),mean,na.rm=T)

# Logistic Regression ----

Contingency50 = summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(Probe_data,(Probe_data$Contingency==0.5)),na.action=na.omit,family=binomial)) 
print(Contingency50)

Contingency100 = summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(Probe_data,(Probe_data$Contingency==1)),na.action=na.omit,family=binomial)) 
print(Contingency100)

ContingencDiff = summary(glmer(Outcome ~ 1 + Contingency2 + (1 |subjectID),data=(Probe_data),na.action=na.omit,family=binomial)) 
print(ContingencDiff)

Contingency50_CI = CI(Contingency50)
Contingency100_CI = CI(Contingency100)
ContingencyDiff_CI = CI(ContingencDiff,2)

# Summary Dataframe ----
# Mean proportions by subject
Probe_data_sub = as.data.frame(tapply(Probe_data$Outcome,list(Probe_data$subjectID,Probe_data$Contingency),mean,na.rm=T))
Probe_data_sub$SubjectID = rownames(Probe_data_sub)
Probe_data_sub2 = melt(Probe_data_sub, variable.name = 'Contingency', id.vars = 'SubjectID')
Probe_data_sub2$Contingency = factor(Probe_data_sub2$Contingency, labels = c("50%", "100%"))
Probe_data_sub2$ContingencyNum = as.numeric(Probe_data_sub2$Contingency)
Probe_data_sub2$value_diff = Probe_data_sub2$value - rep(Probe_data_sub2$value[Probe_data_sub2$Contingency=="50%"],2)
Probe_data_sub2$value_diff2 = rep(Probe_data_sub2$value[Probe_data_sub2$Contingency=="100%"]-Probe_data_sub2$value[Probe_data_sub2$Contingency=="50%"],2)

# Group summary statistics
SummaryStats = as.data.frame(tapply(Probe_data_sub2$value,Probe_data_sub2$Contingency,mean))
colnames(SummaryStats) = "Mean"
SummaryStats$Contingency = factor(row.names(SummaryStats), levels = c("50%","100%"))
SummaryStats$SE = (tapply(Probe_data_sub2$value,Probe_data_sub2$Contingency,se))
SummaryStats$ymin = SummaryStats$Mean - SummaryStats$SE
SummaryStats$ymax = SummaryStats$Mean + SummaryStats$SE
SummaryStats$p = c(Contingency50_CI$p,Contingency100_CI$p)
SummaryStats$asterisk = ""
SummaryStats$asterisk[SummaryStats$p <.05] = "*"
SummaryStats$asterisk[SummaryStats$p <.01] = "**"
SummaryStats$asterisk[SummaryStats$p <.001] = "***"


# Box Plot: colored by slope ----
ggplot(data = Probe_data_sub2, aes(x = Contingency, y = value, fill = Contingency)) +
  geom_boxplot(data = Probe_data_sub2, aes(x = Contingency, y = value, group = Contingency)) + 
  geom_jitter(size = 3, width = 0.1, alpha = 0.7) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,1),  name = "Proportion of trials Go stimuli were chosen")

# Box Plot: colored by slope ----
ggplot(data = Probe_data_sub2, aes(x = Contingency, y = value, group = SubjectID,color = value_diff2)) +
  geom_boxplot(data = Probe_data_sub2, aes(x = Contingency, y = value, group = Contingency), fill = "gray") + 
  geom_point(size = 3) + 
  geom_line() +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  scale_y_continuous(limits = c(0,1),  name = "Proportion of trials Go stimuli were chosen")

# Bar Plot ----
bar_plot = ggplot() +
  geom_bar(data = SummaryStats, aes(x = Contingency, y = Mean, fill = Contingency), stat = 'identity', color = "black", width=.5) + 
  set.seed(0) + geom_jitter(data = Probe_data_sub2, aes(x = Contingency, y = value, group = SubjectID, color= Contingency), size = 2, width = 0.1,  alpha= 0.4) +
  set.seed(0) + geom_jitter(data = Probe_data_sub2, aes(x = Contingency, y = value, group = SubjectID), color = "black", size = 1, width = 0.1, alpha= 0.4) +
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
dev.new()
bar_plot


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
  scale_y_continuous(limits = c(0,1), breaks=seq(0, 1, 0.1),expand = c(0,0),  name = "Proportion of trials Go stimuli were chosen") +
  xlab("Value category") +
  theme_bw() +
  theme(legend.position="top") + 
  geom_errorbar(data = SummaryStatsValue, aes(ymin=ymin, ymax=ymax,x = PairType2, group = Contingency2), position = position_dodge(width=0.8), width = 0.4) + 
  geom_abline(intercept = (0.5),slope=0,linetype =2) # chance level 50% reference line

#geom_text(data = SummaryStatsValue, aes(x=PairType2, y = ymax+0.05, label = asterisk), size = 5) +
#  set.seed(0) + geom_jitter(data = subset(Probe_data_sub_value,ContingencyNum ==2), aes(x = PairType+0.2, y = Mean, color= Contingency2), size = 2, width = 0.1,  alpha= 0.4) + 
#  set.seed(0) + geom_jitter(data = subset(Probe_data_sub_value,ContingencyNum ==2), aes(x = PairType+0.2, y = Mean, color= Contingency2),color = "black", size = 2, width = 0.1,  alpha= 0.4) +
#  set.seed(0) + geom_jitter(data = subset(Probe_data_sub_value,ContingencyNum ==1), aes(x = PairType-0.2, y = Mean, color= Contingency2), size = 2, width = 0.1,  alpha= 0.4) + 
#  set.seed(0) + geom_jitter(data = subset(Probe_data_sub_value,ContingencyNum ==1), aes(x = PairType-0.2, y = Mean, color= Contingency2),color = "black", size = 2, width = 0.1,  alpha= 0.4)  

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





