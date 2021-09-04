
library(rstudioapi)  
library(tidyverse)


# Clear environment
rm(list=ls())

# Get current path
script_path = getActiveDocumentContext()$path
pwd = dirname(script_path)
pathSplit=strsplit(script_path, "/")
pathSplit=pathSplit[[1]]
main_path2=paste0(pathSplit[1:(length(pathSplit)-2)],"/",collapse="")
main_path1 = main_path2 %>% strtrim(nchar(main_path2)-1) %>% paste0("_pilot/")

## Sample
paths = c(paste0(main_path1,"Output/"), paste0(main_path2,"Output/"))
subjects_lists=list(c(101:110,112:121), #Experiment 2 (n = 20)
                    c(101:115,117:118,121:128,131:162,164:165)) #Experiment 2 (n = 59)
# Excluded participants:
# Exp 1: 111 (not really excluded. Code crashed)
# Exp 2:
# 116, 130 - Transitivity
# 163 - FA (max = 37.59%)
# 119, 120 - Miss (max = 38.75%)
# 129 - Technical issues


Ranking_data = c()
CAT_data = c()
Probe_data = c()
# Go over each task in each experiment and read the data
for (exp_i in c(1:2)){
  path = paths[exp_i]
  subjects = subjects_lists[exp_i]
  filelist=c()
  for (s in subjects){
    filelist=c(filelist,Sys.glob(paste(path, "BM_",s,"_ItemRanking*.txt",sep="")))
  }
  Ranking_data_tmp = c()
  for (f in filelist){
    Ranking_data_tmp=rbind(Ranking_data_tmp,read.table(f,header=T,na.strings=c(999,999000)))
  }
  # CAT
  filelist=c()
  for (s in subjects){
    filelist=c(filelist,Sys.glob(paste(path, "BM_",s,"_CAT_*.txt",sep="")))
  }
  CAT_data_tmp=c()
  for (f in filelist){
    CAT_data_tmp=rbind(CAT_data_tmp,read.table(f,header=T,na.strings=c(999,999000)))
  }
  # Probe
  filelist=c()
  for (s in subjects){
    filelist=c(filelist,Sys.glob(paste(path, "BM_",s,"_Probe_*.txt",sep="")))
  }
  Probe_data_tmp=c()
  for (f in filelist){
    Probe_data_tmp=rbind(Probe_data_tmp,read.table(f,header=T,na.strings=c(999,999000)))
  }
  
  Ranking_data_tmp = Ranking_data_tmp %>%
    mutate(experiment_num = exp_i,
           sub_i = as.numeric(as.factor(Subject)))
  CAT_data_tmp = CAT_data_tmp %>%
    mutate(experiment_num = exp_i,
           sub_i = as.numeric(as.factor(Subject)))
  Probe_data_tmp = Probe_data_tmp %>%
    mutate(experiment_num = exp_i,
           sub_i = as.numeric(as.factor(subjectID)))
  
  Ranking_data = rbind(Ranking_data,Ranking_data_tmp)
  CAT_data = rbind(CAT_data,CAT_data_tmp)
  Probe_data = rbind(Probe_data,Probe_data_tmp)
}

experiment_labels = c("Preliminary Exp.","Replication Exp.")
Ranking_data = Ranking_data %>% 
  mutate(experiment = factor(experiment_num, labels = experiment_labels),
         sub_i2 =experiment_num*100+sub_i)
CAT_data = CAT_data %>% 
  mutate(experiment = factor(experiment_num, labels = experiment_labels),
         sub_i2 =experiment_num*100+sub_i)
Probe_data = Probe_data %>% 
  mutate(experiment = factor(experiment_num, labels = experiment_labels),
         sub_i2 =experiment_num*100+sub_i)

save(Ranking_data,CAT_data,Probe_data, file = paste0(pwd,"/Merged_Data.Rda"))

rm(list = ls())
