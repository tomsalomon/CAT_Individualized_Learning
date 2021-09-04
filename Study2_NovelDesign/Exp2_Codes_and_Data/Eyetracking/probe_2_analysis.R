# Load libraries
library(stringr)
library(lme4)
library(lmerTest)
library(rstudioapi)
library("reshape2")
library(ggplot2)

# Clear workspace
rm(list=ls())

# Define these variables:
sessionNum = 1 # 1- session 1 ; 2- follow-up
task_name = "Probe"
gaze_inclusion_thresh = 0.5 # only include scans with at least this proportion of valid eye-tracking 
valid_scans_inclusion_thresh = 2  # only include participants with at least this number of valid scans
session_names = c("First Session", "Follow-Up")
session_name = session_names[sessionNum]

# Standard error function
se = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }

# set paths
current_path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_path)
data_path = './../pre_processed_data/'
output_path = './analyzed_data/'

# Load Data, created by probe_1_organize_data.R script
output_filename_eyetrack_merged = paste0(output_path,'Eyetracking_Task_',task_name,'_Session',sessionNum,'.Rda')
output_filename_behave_merged = paste0(output_path,'Summary_Task_',task_name,'_Session',sessionNum,'.Rda')
load(output_filename_eyetrack_merged)
load(output_filename_behave_merged)

# Count how many scans have valid Eye-Tracking
prop_gaze_by_scan = data.frame(with(behave_merged,tapply(gaze_general , list(subjectID,scan),mean)))
colnames(prop_gaze_by_scan) = str_replace(colnames(prop_gaze_by_scan),"X","")
scans2exclude = which(is.na(prop_gaze_by_scan) | (prop_gaze_by_scan <= gaze_inclusion_thresh),arr.ind = T)
prop_gaze_by_scan$valid_scans = rowSums(prop_gaze_by_scan[,1:4] >= gaze_inclusion_thresh,na.rm = T)

# remove invalid subjects
subs2exclude = rownames(prop_gaze_by_scan[which (prop_gaze_by_scan$valid_scans <= valid_scans_inclusion_thresh),])
print("These participants are excluded due to insufficient number of valid scans:")
print(prop_gaze_by_scan[subs2exclude,])
data2analyze = behave_merged[!behave_merged$subjectID %in% subs2exclude,]
data2analyze$subjectID = droplevels(data2analyze$subjectID)
# remove invalid scans
for (i in 1:nrow(scans2exclude)){
  sub_i = row.names(scans2exclude)[i]
  scan_i = as.numeric(colnames(prop_gaze_by_scan)[scans2exclude[i,"col"]])
  ind2remove = (data2analyze$subjectID == sub_i) & (data2analyze$scan == scan_i)
  data2analyze = data2analyze[!ind2remove,]
}

# Only analyze HV and LV (without sannity trials)
data2analyze = data2analyze[data2analyze$PairType<=2,]
data2analyze = data2analyze[!is.na(data2analyze$Outcome),]
data2analyze$PairType2 = factor(data2analyze$PairType2, labels = c("High Value", "Low Value"))
data2analyze$gaze_go2 = data2analyze$gaze_go/(data2analyze$gaze_go + data2analyze$gaze_nogo) 
data2analyze$gaze_nogo2 = data2analyze$gaze_nogo/(data2analyze$gaze_go + data2analyze$gaze_nogo)

prop_gaze_by_scan2 = data.frame(with(data2analyze,tapply(gaze_general , list(subjectID,scan),mean, na.rm = TRUE)))
colnames(prop_gaze_by_scan2) = str_replace(colnames(prop_gaze_by_scan2),"X","")
prop_gaze_by_scan2$valid_scans = rowSums(prop_gaze_by_scan2[,1:4] >= gaze_inclusion_thresh,na.rm = T)
scans2exclude2 = which(is.na(prop_gaze_by_scan2) | (prop_gaze_by_scan2 <= gaze_inclusion_thresh),arr.ind = T)
n_valid = length(unique(data2analyze$subjectID))
n_recoreded = length(unique(behave_merged$subjectID))
n_excluded_subs = length(subs2exclude)
n_excluded_scans = nrow(scans2exclude2)
print(paste0("Recorded n = ",n_recoreded,"; Valid n = ",n_valid,"; Excluded n = ",n_excluded_subs,"; Excluded scans: ",n_excluded_scans))

data2analyze$gaze_chosen = NA
data2analyze$gaze_chosen[data2analyze$Outcome == 1] = data2analyze$gaze_go[data2analyze$Outcome == 1]
data2analyze$gaze_chosen[data2analyze$Outcome == 0] = data2analyze$gaze_nogo[data2analyze$Outcome == 0]
data2analyze$gaze_unchosen = NA
data2analyze$gaze_unchosen[data2analyze$Outcome == 1] = data2analyze$gaze_nogo[data2analyze$Outcome == 1]
data2analyze$gaze_unchosen[data2analyze$Outcome == 0] = data2analyze$gaze_go[data2analyze$Outcome == 0]

Data_by_sub_go = aggregate(gaze_go  ~ subjectID + PairType2 + Outcome, data = data2analyze, FUN =  mean, na.rm = TRUE)
Data_by_sub_nogo = aggregate(gaze_nogo ~ subjectID + PairType2 + Outcome, data = data2analyze, FUN =  mean, na.rm = TRUE)
Data_by_sub = cbind(Data_by_sub_go,gaze_nogo = Data_by_sub_nogo$gaze_nogo)
Data_by_sub$gaze_diff = Data_by_sub$gaze_go - Data_by_sub$gaze_nogo 

# Summary - with value category
summary_go_m = aggregate(gaze_go ~ PairType2 + Outcome, data = Data_by_sub, FUN =  mean, na.rm = TRUE)
summary_go_sd = aggregate(gaze_go ~ PairType2 + Outcome, data = Data_by_sub, FUN =  sd, na.rm = TRUE)
summary_go = cbind(summary_go_m, sd = summary_go_sd$gaze_go)
summary_go$viewed_stim = "Go"
colnames(summary_go)[[3]] = "gaze"

summary_nogo_m = aggregate(gaze_nogo ~ PairType2 + Outcome, data = Data_by_sub, FUN =  mean, na.rm = TRUE)
summary_nogo_sd = aggregate(gaze_nogo ~ PairType2 + Outcome, data = Data_by_sub, FUN =  sd, na.rm = TRUE)
summary_nogo = cbind(summary_nogo_m, sd = summary_nogo_sd$gaze_nogo)
summary_nogo$viewed_stim = "No-Go"
colnames(summary_nogo)[[3]] = "gaze"

summary_prop = rbind(summary_go,summary_nogo)
n = length(unique(data2analyze$subjectID))
summary_prop$se = summary_prop$sd / sqrt(n)
summary_prop$ymin = summary_prop$gaze - summary_prop$se
summary_prop$ymax = summary_prop$gaze + summary_prop$se
summary_prop$ChosenStimulus = factor(1-summary_prop$Outcome, labels = c("Go item\nWas Chosen", "No-Go item\nWas Chosen"))
summary_prop$PairType3 = factor(summary_prop$PairType2, labels = c("High-Value", "Low-Value"))
summary_prop$line_h = max(summary_prop$ymax) + 0.02 

summary_diff_m = aggregate(gaze_diff ~ PairType2 + Outcome, data = Data_by_sub, FUN =  mean, na.rm = TRUE)
summary_diff_sd = aggregate(gaze_diff ~ PairType2 + Outcome, data = Data_by_sub, FUN =  sd, na.rm = TRUE)
summary_diff = cbind(summary_diff_m, sd = summary_diff_sd$gaze_diff)
#summary_diff$viewed_stim = "No-Go"
colnames(summary_diff)[[3]] = "gaze"
summary_diff$se = summary_diff$sd / sqrt(n)
summary_diff$ChosenStimulus = factor(1-summary_diff$Outcome, labels = c("Go item\nWas Chosen", "No-Go item\nWas Chosen"))

HV_nogo_chosen = summary(lmer(gaze_delta ~ 1 + (1|subjectID), data = subset(data2analyze, PairType == 1 & Outcome == 0, na.action = na.omit)))
LV_nogo_chosen = summary(lmer(gaze_delta ~ 1 + (1|subjectID), data = subset(data2analyze, PairType == 2 & Outcome == 0, na.action = na.omit)))
HV_go_chosen = summary(lmer(gaze_delta ~ 1 + (1|subjectID), data = subset(data2analyze, PairType == 1 & Outcome == 1, na.action = na.omit)))
LV_go_chosen = summary(lmer(gaze_delta ~ 1 + (1|subjectID), data = subset(data2analyze, PairType == 2 & Outcome == 1, na.action = na.omit)))
HVLV_nogo_chosen = summary(lmer(gaze_delta ~ 1 + PairType2 + (1|subjectID), data = subset(data2analyze, PairType <= 2 & Outcome == 0, na.action = na.omit)))
HVLV_go_chosen = summary(lmer(gaze_delta ~ 1 + PairType2 + (1|subjectID), data = subset(data2analyze, PairType <= 2 & Outcome == 1, na.action = na.omit)))
summary_diff_stats = (rbind(HV_nogo_chosen$coefficients, LV_nogo_chosen$coefficients, HV_go_chosen$coefficients, LV_go_chosen$coefficients))
colnames(summary_diff_stats)[5] = "p"
row.names(summary_diff_stats) = 1:4
summary_diff = cbind(summary_diff,summary_diff_stats)
summary_diff$asterisk = ""
summary_diff$asterisk[summary_diff$p<.05] = "*"
summary_diff$asterisk[summary_diff$p<.01] = "**"
summary_diff$asterisk[summary_diff$p<.001] = "***"
summary_diff$viewed_stim = NA
summary_diff$line_h = summary_prop$line_h[1]
summary_diff$PairType3 = factor(summary_diff$PairType2, labels = c("High-Value", "Low-Value"))

plot_proportions = 
  ggplot(data = summary_prop, aes(x = PairType3, y = gaze, fill = viewed_stim)) +
  facet_grid(. ~ ChosenStimulus) +
  geom_bar(stat = 'identity', aes(fill = viewed_stim), position=position_dodge(.9), color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position=position_dodge(.9), width = 0.2) +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, .5)) + 
  scale_fill_manual(values=c("#585858","#D8D8D8")) + # color of bars
  theme(legend.position="top") +
  labs (x = "Value category", y = "Mean proportion of gaze time", fill = "Viewed stimulus", title =  session_name) + 
  theme(plot.title = element_text(hjust = 0.5)) + # center align title
  geom_errorbar(aes(ymin = line_h, ymax = line_h),width = 0.5) + 
  geom_text(data = summary_diff,aes(x = PairType3, y = line_h + .01, label = asterisk))

with(data2analyze,tapply(gaze_go,list(PairType2,Outcome),mean,na.rm=T))
with(data2analyze,tapply(gaze_nogo,list(PairType2,Outcome),mean,na.rm=T))

with(data2analyze,tapply(gaze_go,list(Outcome),mean,na.rm=T))
with(data2analyze,tapply(gaze_nogo,list(Outcome),mean,na.rm=T))

with(data2analyze,tapply(gaze_chosen,list(Outcome),mean,na.rm=T))
with(data2analyze,tapply(gaze_unchosen,list(Outcome),mean,na.rm=T))

dev.new(width=1, height=1)
plot_proportions
# Save plot as pdf
pdf(file=paste0('Probe_EyeTracking_plot_session',sessionNum,'.pdf'), width=7, height=5)
print(plot_proportions)
dev.off()

# Stats: prop. viewing Go when Go versus No-Go were chosen
print("prop. viewing Go when Go were chosen")
summary(lmer(gaze_delta ~ 1 + (1|subjectID), data = subset(data2analyze, Outcome ==1)))
print("prop. viewing Go when NoGo were chosen")
summary(lmer(gaze_delta ~ 1 + (1|subjectID), data = subset(data2analyze, Outcome ==0)))

# Models similar to Schonberg et al. (2014)
print("Difference in Go versus NoGo when the items were chosen")
summary(lmer(gaze_chosen ~ 1 + Outcome + (1 + Outcome|subjectID), data = data2analyze))
print("Difference in Go versus NoGo when the items were not chosen")
summary(lmer(gaze_unchosen ~ 1 + Outcome + (1 + Outcome|subjectID), data = data2analyze))

print("High versus Low value differential effect")
summary(lmer(gaze_delta ~ 1 + PairType2 + (1 + PairType2|subjectID), data = data2analyze))


# Plot - without separation to value category
summary_go_m2 = aggregate(gaze_go ~ Outcome, data = Data_by_sub, FUN =  mean, na.rm = TRUE)
summary_go_sd2 = aggregate(gaze_go ~ Outcome, data = Data_by_sub, FUN =  sd, na.rm = TRUE)
summary_go2 = cbind(summary_go_m2, sd = summary_go_sd2$gaze_go)
summary_go2$viewed_stim = "Go"
colnames(summary_go2)[[2]] = "gaze"

summary_nogo_m2 = aggregate(gaze_nogo ~ Outcome, data = Data_by_sub, FUN =  mean, na.rm = TRUE)
summary_nogo_sd2 = aggregate(gaze_nogo ~ Outcome, data = Data_by_sub, FUN =  sd, na.rm = TRUE)
summary_nogo2 = cbind(summary_nogo_m2, sd = summary_nogo_sd2$gaze_nogo)
summary_nogo2$viewed_stim = "No-Go"
colnames(summary_nogo2)[[2]] = "gaze"

summary_prop2 = rbind(summary_go2,summary_nogo2)
n = length(unique(data2analyze$subjectID))
summary_prop2$se = summary_prop2$sd / sqrt(n)
summary_prop2$ymin = summary_prop2$gaze - summary_prop2$se
summary_prop2$ymax = summary_prop2$gaze + summary_prop2$se
summary_prop2$ChosenStimulus = factor(1-summary_prop2$Outcome, labels = c("Go item\nWas Chosen", "No-Go item\nWas Chosen"))
summary_prop2$line_h = max(summary_prop2$ymax) + 0.02 


summary_diff_m2 = aggregate(gaze_diff ~ Outcome, data = Data_by_sub, FUN =  mean, na.rm = TRUE)
summary_diff_sd2 = aggregate(gaze_diff ~ Outcome, data = Data_by_sub, FUN =  sd, na.rm = TRUE)
summary_diff2 = cbind(summary_diff_m2, sd = summary_diff_sd2$gaze_diff)
#summary_diff$viewed_stim = "No-Go"
colnames(summary_diff2)[[2]] = "gaze"
summary_diff2$se = summary_diff2$sd / sqrt(n)
summary_diff2$ChosenStimulus = factor(1-summary_diff2$Outcome, labels = c("Go item\nWas Chosen", "No-Go item\nWas Chosen"))

go_chosen2 = summary(lmer(gaze_delta ~ 1 + (1|subjectID), data = subset(data2analyze, Outcome == 1, na.action = na.omit)))
nogo_chosen2 = summary(lmer(gaze_delta ~ 1 + (1|subjectID), data = subset(data2analyze, Outcome == 0, na.action = na.omit)))
summary_diff_stats2 = (rbind(go_chosen2$coefficients, nogo_chosen2$coefficients))
colnames(summary_diff_stats2)[5] = "p"
row.names(summary_diff_stats2) = 1:2
summary_diff = cbind(summary_diff2,summary_diff_stats2)
summary_diff2$asterisk = ""
summary_diff2$asterisk[summary_diff$p<.05] = "*"
summary_diff2$asterisk[summary_diff$p<.01] = "**"
summary_diff2$asterisk[summary_diff$p<.001] = "***"
summary_diff2$viewed_stim = NA
summary_diff2$line_h = summary_prop2$line_h[1]


plot_proportions2 = 
  ggplot(data = summary_prop2, aes(x = ChosenStimulus, y = gaze, fill = viewed_stim)) +
  geom_bar(stat = 'identity', aes(fill = viewed_stim), position=position_dodge(.9), color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position=position_dodge(.9), width = 0.2) +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, .5)) + 
  scale_fill_manual(values=c("#585858","#D8D8D8")) + # color of bars
  theme(legend.position="top") +
  labs (x = "Chosen Stimulus", y = "Mean proportion of gaze time", fill = "Viewed stimulus", title =  session_name) + 
  theme(plot.title = element_text(hjust = 0.5)) + # center align title
  geom_errorbar(aes(ymin = line_h, ymax = line_h),width = 0.5) + 
  geom_text(data = summary_diff2,aes(x = ChosenStimulus, y = line_h + .01, label = asterisk))

dev.new(width=1, height=1)
plot_proportions2
# Save plot as pdf
pdf(file=paste0('Probe_EyeTracking_plot_session',sessionNum,'_v2.pdf'), width=5, height=5)
print(plot_proportions2)
dev.off()




eyetrack_merged$unique_trial = paste0(eyetrack_merged$subjectID,"_",eyetrack_merged$scan,"_",eyetrack_merged$trial)
data2analyze$unique_trial = paste0(data2analyze$subjectID,"_",data2analyze$scan,"_",data2analyze$trial)
trials2analyze = unique(data2analyze$unique_trial)
eyetracking_data = eyetrack_merged[eyetrack_merged$unique_trial %in% data2analyze$unique_trial,]
eyetracking_data2 = merge(eyetracking_data,data2analyze,by = "unique_trial")

ggplot(data = subset(eyetracking_data2,Response == 'y' & trial.x >=30, na.rm = TRUE), aes(x = x, y =y)) +
  #annotation_raster(image, -Inf, Inf, -Inf, Inf, interpolate = TRUE) +
  stat_density2d(aes(x = x, y =y, fill = ..level.., alpha = ..level..), size= 100, bins= 50, geom='polygon') + 
  theme_bw() +scale_fill_gradient(low = "blue", high = "red") +
  scale_alpha_continuous(range=c(0.01,0.5) , guide = FALSE) +
  coord_cartesian(xlim= c(0,1920), ylim= c(0,1080))+
  scale_y_reverse() + 
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())


