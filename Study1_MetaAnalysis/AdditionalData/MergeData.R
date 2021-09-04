
# Load rquired R packages
library(rstudioapi) # to define current working directory
library(lme4)
library(plyr)
library(ggplot2)

# clear workspace
rm(list=ls())

# Define the current script location as the working directory
pwd = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(pwd)
data_path_main = './RawData/*'
experiment_paths = Sys.glob(paste0(data_path_main,'Exp*'))

# count_unique function
count_unique = function(x) {length(unique(x))}
count_na = function(x) {sum(is.na(x))}
count_na_prop = function(x) {sum(is.na(x)) / length(x)}


probe_data = c()
training_data = c()

old_headers = c("subjid","runnum","itemname","onsettime","bidindex","BDMtrialIndex","Audiotime","AudioTime","Ausdiotime","fixationtime","runtrial",
                "IsLefthigh","shuff_trialType","trialtype","bidValue","bidLeft","bidRight","bid")
new_headers = c("subjectID","runNum","itemName","onsetTime","valueIndex","valueIndex","CueTime","CueTime","CueTime","fixationTime","trial",
                "IsleftGo","trialType","trialType","value","valueLeft","valueRight","value")
header2remove = c("TypeLeft","TypeRight","X.ixationTime","scanner","test_comp","timeFixLeft","timeFixMid","timeFixRight","numLeftFix","numMidFix","numRightFix","firstFix" ,"firstFixTime")

for (exp_num in 1:length(experiment_paths)){
  exp_path = experiment_paths[exp_num]
  exp_name = substr(exp_path,nchar(data_path_main),nchar(exp_path))
  
  print("")
  print(paste0("Experiment ",exp_num,": ",exp_name))
  print("==============================")
  
  filelist_probe=Sys.glob(paste0(exp_path, "/*probe_block*.txt"))
  if (length(filelist_probe)==0){
    filelist_probe=Sys.glob(paste0(exp_path, "/*boostprobe*.txt"))
  }
  filelist_training = Sys.glob(paste0(exp_path, "/*training*run*.txt"))
  old_format = FALSE
  if (length(filelist_training)==0){
    old_format = TRUE
    filelist_training = Sys.glob(paste0(exp_path, "/*boosting*run*.txt"))
  }
  
  if (exp_num==26){ # Exp026_NN7_Snacks_LV had some participatns follow-up probe as blocks 5-6
    filelist_probe_sesII = c(Sys.glob(paste0(exp_path, "/*probe_block_5*.txt")),Sys.glob(paste0(exp_path, "/*probe_block_6*.txt")))
    filelist_probe = filelist_probe[!filelist_probe %in% filelist_probe_sesII]
  }
  # read data
  tmp_data_probe=c()
  for (f in filelist_probe){
    tmp_data_probe_i = read.table(f,header=T,na.strings=c(999,999000))
    # Exp026_NN7_Snacks_LV had for some participatns a binary variable called "bidIndexLeft" instead of "bidValueLeft"
    if (exp_num==26) {
      if("bidIndexLeft" %in% colnames(tmp_data_probe_i)){ 
        tmp_data_probe_i$bidIndexLeft = NA;
        tmp_data_probe_i$bidIndexRight = NA;
        colnames(tmp_data_probe_i)%in%"bidIndexLeft"
        colnames(tmp_data_probe_i) = gsub("bidIndex","bidValue",colnames(tmp_data_probe_i))
      }
    }
    if (exp_num==15) { #deal with subjects saved as numbers in CAT OFC older population (Canda)
    tmp_data_probe_i$subjectID = factor(tmp_data_probe_i$subjectID) 
    }
    #print(f)
    #print(count_na(tmp_data_probe_i$subjectID))
    #coln = rbind(coln,colnames(tmp_data_probe_i))
    tmp_data_probe=rbind(tmp_data_probe,tmp_data_probe_i)
  }

  tmp_data_training=c()
  for (f in filelist_training){
    tmp_data_training_i = read.table(f,header=T,na.strings=c(999,999000))
    # deal with demo run save in Adolescence experiment
    if (nrow(tmp_data_training_i)<=10){
      filelist_training=filelist_training[!filelist_training==f]
    }
  }
  
  for (f in filelist_training){
    tmp_data_training_i = read.table(f,header=T,na.strings=c(999,999000))
    run_num = tmp_data_training_i$runNum[1]
    if (is.null(run_num)) { 
      run_num = tmp_data_training_i$runnum[1]
    }
    n_training_trials = nrow(tmp_data_training_i)
    tmp_data_training_i$trial = (run_num - 1)*n_training_trials + 1:n_training_trials
    # deal with subjects saved as numbers in CAT OFC older population (Canda) and few other cases
    if (is.numeric(tmp_data_training_i$subjectID[1])) { 
      tmp_data_training_i$subjectID = factor(tmp_data_training_i$subjectID) 
    }
    tmp_data_training=rbind(tmp_data_training,tmp_data_training_i)
  }
  
  tmp_data_probe = tmp_data_probe[!colnames(tmp_data_probe) %in% header2remove]
  tmp_data_training = tmp_data_training[!colnames(tmp_data_training) %in% header2remove]
  for (str_i in 1:length(old_headers)){
    colnames(tmp_data_training) = gsub(old_headers[str_i],new_headers[str_i],colnames(tmp_data_training), ignore.case = TRUE)
    colnames(tmp_data_probe) = gsub(old_headers[str_i],new_headers[str_i],colnames(tmp_data_probe), ignore.case = TRUE)
    #col2replace_training = colnames(tmp_data_training) == old_headers[str_i]
    #col2replace_probe = colnames(tmp_data_probe) == old_headers[str_i]
    #if (sum(col2replace_training) >=1){
    #  colnames(tmp_data_training)[colnames(tmp_data_training) == old_headers[str_i]] = new_headers[str_i]
    #}
    #if (sum(col2replace_probe) >=1){
    #  colnames(tmp_data_probe)[colnames(tmp_data_probe) == old_headers[str_i]] = new_headers[str_i]
    #}
  }
  
  # count average number of files and number of trials per subject (should be nice and round number)
  n = count_unique(tmp_data_probe$subjectID)
  files_per_subject_probe = length(filelist_probe)/n
  files_per_subject_training = length(filelist_training)/n
  trials_per_subject_probe = nrow(tmp_data_probe)/n
  trials_per_subject_training = nrow(tmp_data_training)/n
  
  print(paste0("n = ",n))
  print(paste0("files per subject: probe = ",files_per_subject_probe," training = ",files_per_subject_training))
  print(paste0("trials per subject: probe = ",trials_per_subject_probe," training = ",trials_per_subject_training))
  
  tmp_data_probe$sub_i = as.numeric(as.factor(tmp_data_probe$subjectID))
  tmp_data_training$sub_i = as.numeric(as.factor(tmp_data_training$subjectID))
  tmp_data_probe$Experiment=exp_num
  tmp_data_probe$ExperimentName=exp_name
  tmp_data_training$Experiment=exp_num
  tmp_data_training$ExperimentName=exp_name
  
  if (exp_num ==27) {
    tmp_data_training$CueTime = tmp_data_training$CueTime*1000
    tmp_data_training$RT = tmp_data_training$RT*1000
  }
  probe_data=rbind.fill(probe_data,tmp_data_probe)
  training_data=rbind.fill(training_data,tmp_data_training)
}

probe_data$ExperimentNameFull=paste("Exp. ",formatC(probe_data$Experiment,width=2,flag = "0"),": ",probe_data$ExperimentName,sep="")
training_data$ExperimentNameFull=paste("Exp. ",formatC(training_data$Experiment,width=2,flag = "0"),": ",training_data$ExperimentName,sep="")

for (exp_num in 1:length(experiment_paths)){
  tmp_dat  = subset(probe_data, Experiment == exp_num )
  print("")
  print(tmp_dat$ExperimentNameFull[1])
  print(tapply(tmp_dat$Outcome, tmp_dat$PairType, mean, na.rm = T))
  #tmp = tmp_dat[1:10,c("IsleftGo","Response","Outcome")]
  #print(tmp)
}

probe_data$choseLeft = NA
probe_data$choseLeft[probe_data$Response %in% c("u","b")] = 1
probe_data$choseLeft[probe_data$Response %in% c("i","y")] = 0
probe_data$Outcome = NA
probe_data$Outcome[probe_data$choseLeft == probe_data$IsleftGo] = 1
probe_data$Outcome[!(probe_data$choseLeft == probe_data$IsleftGo)] = 0

# in experiment 27 each scan (with 2 run) was counted as 1 run. fix it by splitting each run into 2.
run_exp27 = training_data$runNum[training_data$Experiment==27]
n_subs_exp27 = max(training_data$sub_i[training_data$Experiment==27])
n_trials_exp27 = sum(run_exp27==1)/n_subs_exp27
run_exp27_correct = (run_exp27-1)*2 + c(rep(1,n_trials_exp27/2),rep(2,n_trials_exp27/2))
training_data$runNum[training_data$Experiment==27] = run_exp27_correct
# in experiment 27 only LV nogo was coded as [11,12,21,22]  instead of [11,12,22,24]
trialtype_exp27 = training_data$trialType[training_data$Experiment==27] 
trialtype_exp27_correct = trialtype_exp27
trialtype_exp27_correct[trialtype_exp27==22] = 24
trialtype_exp27_correct[trialtype_exp27==21] = 22
training_data$trialType[training_data$Experiment==27] = trialtype_exp27_correct

# In few cases, at a very high ladder the cue was so late that it did not appeared at all. for these set the planned time as the cue time
training_data$CueTime_original = training_data$CueTime
training_data$CueTime[is.na(training_data$CueTime_original) & training_data$trialType==11] = training_data$ladder1[is.na(training_data$CueTime_original) & training_data$trialType==11]
training_data$CueTime[is.na(training_data$CueTime_original) & training_data$trialType==22] = training_data$ladder2[is.na(training_data$CueTime_original) & training_data$trialType==22]


# Ignore RT > 2000 or RT < 0 (are all errors in code)
training_data$RT[training_data$RT>2000 | training_data$RT<=0] = NaN

training_data$RT_eff = training_data$RT - training_data$CueTime
training_data$RT_eff[abs(training_data$RT_eff)>1000] = NaN
training_data$WasCue = 0
training_data$WasCue[training_data$CueTime>=0] = 1
training_data$WasResponse = 0
training_data$WasResponse[training_data$RT>=40] = 1
training_data$WasResponse[training_data$RT<40] = NA # ignore very short RT
training_data$TP = training_data$WasCue & training_data$WasResponse
training_data$TN = !training_data$WasCue & !training_data$WasResponse
training_data$FP = !training_data$WasCue & training_data$WasResponse
training_data$FN = training_data$WasCue & !training_data$WasResponse
training_data$TP [training_data$WasCue]

probe_data$subjectID2 = as.factor(sprintf("E%02i_S%02i",probe_data$Experiment,probe_data$sub_i))
training_data$subjectID2 = as.factor(sprintf("E%02i_S%02i",training_data$Experiment,training_data$sub_i))
training_data$runNum2 = as.factor(training_data$runNum)

tapply(training_data$FP, training_data$Experiment, mean, na.rm = T) / tapply(!training_data$WasCue, training_data$Experiment, mean, na.rm = T)
tapply(training_data$FN, training_data$Experiment, mean, na.rm = T) / tapply(training_data$WasCue, training_data$Experiment, mean, na.rm = T)

hist(training_data$RT_eff,100)

save(probe_data,file=paste0(pwd,"/probe_data.Rda"))
save(training_data,file=paste0(pwd,"/training_data.Rda"))

# Count number of subjects
tapply(probe_data$subjectID,probe_data$Experiment,count_unique)
# Count invalid entries
tapply(probe_data$subjectID,probe_data$Experiment,count_na)
tapply(probe_data$valueIndexLeft,probe_data$Experiment,count_na)

tapply(probe_data$ExperimentName,probe_data$Experiment,unique)
tapply(probe_data$subjectID,probe_data$Experiment,length)/tapply(probe_data$subjectID,probe_data$Experiment,count_unique)

probe_data$PairType2[probe_data$PairType==1]="High Value"
probe_data$PairType2[probe_data$PairType==2]="Low Value"
probe_data$PairType2[probe_data$PairType==3]="Sanity Go"
probe_data$PairType2[probe_data$PairType==4]="Sanity NoGo"

probe_data$PairType3[probe_data$Experiment<=4 & probe_data$bidIndexLeft %in% c(7:14) & probe_data$bidIndexRight %in% c(7:14)]="Highest 4"
probe_data$PairType3[probe_data$Experiment<=4 & probe_data$bidIndexLeft %in% c(15:22) & probe_data$bidIndexRight %in% c(15:22)]="High middle 4"
probe_data$PairType3[probe_data$Experiment<=4 & probe_data$bidIndexLeft %in% c(39:46) & probe_data$bidIndexRight %in% c(39:46)]="Low middle 4"
probe_data$PairType3[probe_data$Experiment<=4 & probe_data$bidIndexLeft %in% c(47:54) & probe_data$bidIndexRight %in% c(47:54)]="Lowest 4"

training_data2 = training_data[!is.na(training_data$RT_eff),]
training_data2$runNum2 = as.factor(training_data2$runNum)

ggplot(data = training_data2, aes(x = CueTime, y = RT_eff, color = runNum)) +
  geom_point(alpha = 0.1)

ggplot(data = training_data2, aes(x = RT_eff)) + 
  geom_histogram(aes(x=RT_eff),bins = 200)

ggplot(data = training_data2, aes(x = RT_eff, color = runNum2)) + 
  geom_density(aes(x=RT_eff, y=..density.., color=runNum2, alpha = .01))

ggplot(data = training_data2, aes(x = RT, color = runNum2)) + 
  geom_density(aes(x=RT, y=..density.., color=runNum2, alpha = .01))

# plot 
n= data.frame(n = tapply(training_data$subjectID, training_data$Experiment, count_unique))
ggplot(data = n) + geom_histogram(aes(x = n),binwidth = 0.9, color="black", fill="white") + xlim(20,80)

ggplot(data = subset(training_data2,runNum<=4), aes(x = RT_eff, color = runNum2)) + 
  geom_density(aes(x=RT_eff, y=..density.., color=runNum2))
