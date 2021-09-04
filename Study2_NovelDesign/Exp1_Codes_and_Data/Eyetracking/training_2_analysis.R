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
task_name = "Training"
gaze_inclusion_thresh = 0.5 # only include scans with at least this proportion of valid eye-tracking 
valid_scans_inclusion_thresh = 0.5  # only include participants with at least this proportion of valid scans
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
prop_gaze_by_scan = data.frame(with(behave_merged,tapply(gaze_stim , list(subjectID,scan),mean)))
colnames(prop_gaze_by_scan) = str_replace(colnames(prop_gaze_by_scan),"X","")
scans2exclude = which(is.na(prop_gaze_by_scan) | (prop_gaze_by_scan <= gaze_inclusion_thresh),arr.ind = T)
n_scans = length(unique(behave_merged$scan))
valid_scans_inclusion_thresh_n = valid_scans_inclusion_thresh * n_scans
prop_gaze_by_scan$valid_scans = rowSums(prop_gaze_by_scan[,1:n_scans] >= gaze_inclusion_thresh,na.rm = T)

# remove invalid subjects
subs2exclude = rownames(prop_gaze_by_scan[which (prop_gaze_by_scan$valid_scans <= valid_scans_inclusion_thresh_n),])
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

# Only analyze probe items (without sannity trials and filler)
data2analyze = data2analyze[data2analyze$isprobeitem,]
data2analyze$PairType = 2-data2analyze$ishvitem
data2analyze$PairType2 = factor(data2analyze$PairType, labels = c("High Value", "Low Value"))
data2analyze$gaze_go = NA
data2analyze$gaze_nogo = NA
data2analyze$gaze_go[data2analyze$isgoitem] = data2analyze$gaze_stim[data2analyze$isgoitem]
data2analyze$gaze_nogo[!data2analyze$isgoitem] = data2analyze$gaze_stim[!data2analyze$isgoitem]

prop_gaze_by_scan2 = data.frame(with(data2analyze,tapply(gaze_stim , list(subjectID,scan),mean, na.rm = TRUE)))
colnames(prop_gaze_by_scan2) = str_replace(colnames(prop_gaze_by_scan2),"X","")
prop_gaze_by_scan2$valid_scans = rowSums(prop_gaze_by_scan2[,1:n_scans] >= gaze_inclusion_thresh,na.rm = T)
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

Data_by_sub_go = aggregate(gaze_go  ~ subjectID + PairType, data = data2analyze, FUN =  mean, na.rm = TRUE)
Data_by_sub_nogo = aggregate(gaze_nogo ~ subjectID + PairType, data = data2analyze, FUN =  mean, na.rm = TRUE)
Data_by_sub = cbind(Data_by_sub_go,gaze_nogo = Data_by_sub_nogo$gaze_nogo)
Data_by_sub$gaze_diff = Data_by_sub$gaze_go - Data_by_sub$gaze_nogo 
Data_by_sub$PairType2 = factor(Data_by_sub$PairType,labels = c("High-Value", "Low-Value"))
# Summary - with value category
summary_go_m = aggregate(gaze_go ~ PairType , data = Data_by_sub, FUN =  mean, na.rm = TRUE)
summary_go_sd = aggregate(gaze_go ~ PairType , data = Data_by_sub, FUN =  sd, na.rm = TRUE)
summary_go = cbind(summary_go_m, sd = summary_go_sd$gaze_go)
summary_go$viewed_stim = "Go"
colnames(summary_go)[[2]] = "gaze"

summary_nogo_m = aggregate(gaze_nogo ~ PairType, data = Data_by_sub, FUN =  mean, na.rm = TRUE)
summary_nogo_sd = aggregate(gaze_nogo ~ PairType, data = Data_by_sub, FUN =  sd, na.rm = TRUE)
summary_nogo = cbind(summary_nogo_m, sd = summary_nogo_sd$gaze_nogo)
summary_nogo$viewed_stim = "No-Go"
colnames(summary_nogo)[[2]] = "gaze"

summary_prop = rbind(summary_go,summary_nogo)
n = length(unique(data2analyze$subjectID))
summary_prop$se = summary_prop$sd / sqrt(n)
summary_prop$ymin = summary_prop$gaze - summary_prop$se
summary_prop$ymax = summary_prop$gaze + summary_prop$se
summary_prop$PairType2 = factor(summary_prop$PairType, labels = c("High-Value", "Low-Value"))
summary_prop$line_h = max(summary_prop$ymax) + 0.02 

summary_diff_m = aggregate(gaze_diff ~ PairType, data = Data_by_sub, FUN =  mean, na.rm = TRUE)
summary_diff_sd = aggregate(gaze_diff ~ PairType, data = Data_by_sub, FUN =  sd, na.rm = TRUE)
summary_diff = cbind(summary_diff_m, sd = summary_diff_sd$gaze_diff)
#summary_diff$viewed_stim = "No-Go"
colnames(summary_diff)[[2]] = "gaze"
summary_diff$se = summary_diff$sd / sqrt(n)

HV = summary(lm(gaze_diff ~ 1 , data = subset(Data_by_sub, PairType == 1 , na.action = na.omit)))
LV = summary(lm(gaze_diff ~ 1 , data = subset(Data_by_sub, PairType == 2 , na.action = na.omit)))
HVLV_diff = summary(lmer(gaze_diff ~ 1 + PairType2 + (1 |subjectID), data = Data_by_sub))
summary_diff_stats = (rbind(HV$coefficients, LV$coefficients))
colnames(summary_diff_stats)[4] = "p"
row.names(summary_diff_stats) = 1:nrow(summary_diff_stats)
summary_diff_stats = as.data.frame(summary_diff_stats)
summary_diff_stats$df = c(HV$df[2],LV$df[2])
summary_diff = cbind(summary_diff,summary_diff_stats)
summary_diff$asterisk = ""
summary_diff$asterisk[summary_diff$p<.05] = "*"
summary_diff$asterisk[summary_diff$p<.01] = "**"
summary_diff$asterisk[summary_diff$p<.001] = "***"
summary_diff$viewed_stim = NA
summary_diff$line_h = summary_prop$line_h[1]
summary_diff$PairType2 = factor(summary_diff$PairType, labels = c("High-Value", "Low-Value"))

plot_proportions = ggplot(data = summary_prop, aes(x = PairType2, y = gaze, fill = viewed_stim)) +
  geom_bar(stat = 'identity', aes(fill = viewed_stim), position=position_dodge(.9), color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position=position_dodge(.9), width = 0.2) +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1)) + 
  scale_fill_manual(values=c("#585858","#D8D8D8")) + # color of bars
  theme(legend.position="top") +
  labs (x = "Value category", y = "Mean proportion of gaze time", fill = "Viewed stimulus", title =  session_name) + 
  theme(plot.title = element_text(hjust = 0.5))  # center align title
#  geom_errorbar(aes(ymin = line_h, ymax = line_h),width = 0.5) + 
#  geom_text(data = summary_diff,aes(x = PairType2, y = line_h + .01, label = asterisk))

with(data2analyze,tapply(gaze_go,list(PairType2),mean,na.rm=T))
with(data2analyze,tapply(gaze_nogo,list(PairType2),mean,na.rm=T))

mean(data2analyze$gaze_go,na.rm=T)
mean(data2analyze$gaze_nogo,na.rm=T)


dev.new(width=1, height=1)
plot_proportions
# Save plot as pdf
pdf(file=paste0(task_name,'_EyeTracking_plot_session',sessionNum,'.pdf'), width=5, height=5)
print(plot_proportions)
dev.off()

# Stats: prop. viewing Go when Go versus No-Go were chosen
print("prop. viewing Go versus No-Go")
summary(lmer(gaze_diff ~ 1 + (1|subjectID), data = Data_by_sub))
print("Value category effect")
summary(lmer(gaze_diff ~ 1 + PairType2 + (1 |subjectID), data = Data_by_sub))

# Models similar to Schonberg et al. (2014)
print("Difference in Go versus NoGo when the items were chosen")
summary(lmer(gaze_chosen ~ 1 + Outcome + (1 + Outcome|subjectID), data = data2analyze))
print("Difference in Go versus NoGo when the items were not chosen")
summary(lmer(gaze_unchosen ~ 1 + Outcome + (1 + Outcome|subjectID), data = data2analyze))

Data_by_subscan_go = aggregate(gaze_go  ~ subjectID + PairType + scan, data = data2analyze, FUN =  mean, na.rm = TRUE)
Data_by_subscan_nogo = aggregate(gaze_nogo ~ subjectID + PairType + scan, data = data2analyze, FUN =  mean, na.rm = TRUE)
Data_by_subscan = cbind(Data_by_subscan_go,gaze_nogo = Data_by_subscan_nogo$gaze_nogo)
Data_by_subscan$gaze_diff = Data_by_subscan$gaze_go - Data_by_subscan$gaze_nogo 
Data_by_subscan$PairType2 = factor(Data_by_subscan$PairType,labels = c("High-Value", "Low-Value"))
Data_by_subscan$scan = as.numeric(Data_by_subscan$scan)


with(Data_by_subscan, tapply(gaze_go, list(scan,PairType2), mean))
with(Data_by_subscan, tapply(gaze_nogo, list(scan,PairType2), mean))
with(Data_by_subscan, tapply(gaze_diff, list(scan,PairType2), mean))

# Summary - with scan category
summary_go_m = aggregate(gaze_go ~ scan , data = Data_by_subscan, FUN =  mean, na.rm = TRUE)
summary_go_sd = aggregate(gaze_go ~ scan , data = Data_by_subscan, FUN =  sd, na.rm = TRUE)
summary_go = cbind(summary_go_m, sd = summary_go_sd$gaze_go)
summary_go$viewed_stim = "Go"
colnames(summary_go)[[2]] = "gaze"

summary_nogo_m = aggregate(gaze_nogo ~ scan, data = Data_by_subscan, FUN =  mean, na.rm = TRUE)
summary_nogo_sd = aggregate(gaze_nogo ~ scan, data = Data_by_subscan, FUN =  sd, na.rm = TRUE)
summary_nogo = cbind(summary_nogo_m, sd = summary_nogo_sd$gaze_nogo)
summary_nogo$viewed_stim = "No-Go"
colnames(summary_nogo)[[2]] = "gaze"

summary_prop = rbind(summary_go,summary_nogo)
n = length(unique(data2analyze$subjectID))
summary_prop$se = summary_prop$sd / sqrt(n)
summary_prop$ymin = summary_prop$gaze - summary_prop$se
summary_prop$ymax = summary_prop$gaze + summary_prop$se
summary_prop$line_h = max(summary_prop$ymax) + 0.02

plot_proportions2 = ggplot(data = summary_prop, aes(x = scan, y = gaze, fill = viewed_stim)) +
  geom_bar(stat = 'identity', aes(fill = viewed_stim), position=position_dodge(.9), color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position=position_dodge(.9), width = 0.2) +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1)) + 
  scale_fill_manual(values=c("#585858","#D8D8D8")) + # color of bars
  theme(legend.position="top") +
  labs (x = "Scan number", y = "Mean proportion of gaze time", fill = "Viewed stimulus", title =  session_name) + 
  theme(plot.title = element_text(hjust = 0.5))  # center align title


summary_diff_m = aggregate(gaze_diff ~ scan , data = Data_by_subscan, FUN =  mean, na.rm = TRUE)
summary_diff_sd = aggregate(gaze_diff ~ scan , data = Data_by_subscan, FUN =  sd, na.rm = TRUE)
summary_diff = cbind(summary_diff_m, sd = summary_diff_sd$gaze_diff)
summary_diff$viewed_stim = NA
colnames(summary_diff)[[2]] = "gaze"
summary_diff$se = summary_diff$sd / sqrt(n)
summary_diff$ymin = summary_diff$gaze - summary_diff$se
summary_diff$ymax = summary_diff$gaze + summary_diff$se
summary_diff$line_h = max(summary_diff$ymax) + 0.02
summary_diff$go = summary_go_m[,2]
summary_diff$nogo = summary_nogo_m[,2]

print(summary_diff)
print("Time effect")
summary(lmer(gaze_diff ~ 1 + scan + (1 + scan|subjectID), data = Data_by_subscan))

plot_proportions3 = ggplot(data = summary_diff, aes(x = scan, y = gaze)) +
  geom_bar(stat = 'identity', position=position_dodge(.9), color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position=position_dodge(.9), width = 0.2) +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent, limits = c(-.1, .1)) + 
  scale_fill_manual(values=c("#585858","#D8D8D8")) + # color of bars
  theme(legend.position="top") +
  labs (x = "Scan number", y = "Mean proportion of gaze time (Go minus No-Go)", fill = "Viewed stimulus", title =  session_name) + 
  theme(plot.title = element_text(hjust = 0.5))  # center align title

t_shift <- scales::trans_new("shift",
                             transform = function(x) {x-.8},
                             inverse = function(x) {x+.8})

plot_proportions4 = ggplot(data = summary_prop, aes(x = scan, y = gaze, fill = viewed_stim)) +
  geom_bar(stat = 'identity', aes(fill = viewed_stim), position=position_dodge(.9), color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position=position_dodge(.9), width = 0.2) +
  theme_bw() + 
  geom_line(data = summary_diff,  aes(x = scan, y = gaze + .85, color = "Delta Gaze"), size = 1.5) +
  #geom_errorbar(data = summary_diff,aes(ymin = ymin, ymax = ymax,color = "red"), width = 0.2) +
  #scale_y_continuous(labels = scales::percent, trans = t_shift, sec.axis = sec_axis(~.+(-.85), name = expression(paste(Delta," Gaze (Go minus No-Go)")), labels = scales::percent)) + 
  scale_y_continuous(labels = scales::percent, sec.axis = sec_axis(~.+(-.85), name = expression(paste(Delta," Gaze (Go minus No-Go)")), labels = scales::percent)) + 
  scale_fill_manual(values=c("#585858","#D8D8D8")) + # color of bars
  theme(legend.position="top") +
  labs (x = "Scan number", y = "Mean proportion of gaze time", fill = "Viewed stimulus", title =  session_name) + 
  theme(plot.title = element_text(hjust = 0.5)) +  # center align title
  
dev.new(width=1, height=1)
plot_proportions4
# Save plot as pdf
pdf(file=paste0(task_name,'_EyeTracking_plot_session',sessionNum,'_v4.pdf'), width=5, height=5)
print(plot_proportions4)
dev.off()
