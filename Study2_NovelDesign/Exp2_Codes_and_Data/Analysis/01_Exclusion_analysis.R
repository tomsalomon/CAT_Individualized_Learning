
library(lme4)
library(rstudioapi)  
library(tidyverse)

# Clear environment
rm(list=ls())

# Get current path
script_path = getActiveDocumentContext()$path
pathSplit=strsplit(script_path, "/")
pathSplit=pathSplit[[1]]
main_path=paste0(pathSplit[1:(length(pathSplit)-2)],"/",collapse="")

## Sample
path=paste0(main_path,"/Output/")
subjects=c(101:165); # all subjects

# technical issues
subjects_technical_issues = c(129)
subjects=subjects[!subjects%in%subjects_technical_issues]
# Define exclusion criteria (z scored)
exclusion_thresh = 3

# Load Data ----
# Binary ranking
filelist=c()
for (s in subjects){
  filelist=c(filelist,Sys.glob(paste(path, "BM_",s,"_ItemRanking*.txt",sep="")))
}
RankingData = c()
for (f in filelist){
  RankingData=rbind(RankingData,read.table(f,header=T,na.strings=c(999,999000)))
}
# CAT
filelist=c()
for (s in subjects){
  filelist=c(filelist,Sys.glob(paste(path, "BM_",s,"_CAT_*.txt",sep="")))
}
CAT_data=c()
for (f in filelist){
  CAT_data=rbind(CAT_data,read.table(f,header=T,na.strings=c(999,999000)))
}

# Transitivity Scores ----
Exclusion_Test = RankingData %>%
  group_by(Subject) %>%
  summarise(TransitivityScore = sd(Rank)) %>%
  mutate(subjectID_num = as.numeric(substr(Subject,4,6)))

# CAT performance ----
CAT_data_summary = CAT_data %>%
  mutate(Response = 1 - is.na(RT),
         Go = factor(Go, labels=c("FA", "TP"))) %>%
  group_by(Subject, Go) %>%
  summarise(prop_responded = mean(Response)) %>%
  spread(Go, prop_responded) %>%
  mutate(Miss = 1-TP)

# Define outliers ----
Exclusion_Test = left_join(Exclusion_Test, CAT_data_summary) %>%
  mutate(ZTransitivityScore = as.numeric(scale(TransitivityScore)),
         ZFA = as.numeric(scale(FA)),
         ZMiss = as.numeric(scale(Miss)),
         
         ExcludeTransitivityScore = ZTransitivityScore <= -exclusion_thresh,
         ExcludeFA = ZFA >= exclusion_thresh,
         ExcludeMiss = ZMiss >= exclusion_thresh,
         
         valid = (!(ExcludeTransitivityScore | ExcludeFA | ExcludeMiss))
         ) %>%
  as.data.frame()


# Exclude outliers ----
valid_codes = Exclusion_Test %>% filter(valid) %>%
  pull(subjectID_num)
n_valid = sum(Exclusion_Test$valid)
# add artificial edges to the subjects lists
valid_codes=c(0,valid_codes,10000)
# Go over valid subjects and rewrite as trimmed lost to be copied to further analysis scripts
valids_subjects_text = paste0('# Valid participants (n = ',n_valid,'):\nsubjects = c(');
for (i in 2:(length(valid_codes)-1)){
  starting_num = (valid_codes[i]-valid_codes[i-1])>1
  ending_num = (valid_codes[i+1]-valid_codes[i])>1
  if (starting_num & ending_num) {
    text_i=paste0(',',valid_codes[i])
  } else if (starting_num & (!ending_num)) {
    text_i=paste0(',',valid_codes[i],':')
  } else if ((!starting_num) & ending_num) {
    text_i=paste0('',valid_codes[i],'')
  } else if ((!starting_num) & (!ending_num)) {
    text_i=''
  } 
  if (i==2){ # remove ',' for the first subject in list
    text_i = substr(text_i,2,nchar(text_i))
  }
  valids_subjects_text = paste0(valids_subjects_text,text_i)
}
valids_subjects_text = paste0(valids_subjects_text,')')

invalids_subjects_text ='# Excluded participants:'
if (sum(Exclusion_Test$ExcludeTransitivityScore)>0){
  invalids_subjects_text = paste0(invalids_subjects_text, 
                                 '\n# ',paste0(Exclusion_Test$subjectID_num[Exclusion_Test$ExcludeTransitivityScore],collapse = ', '),' - Transitivity')}
if (sum(Exclusion_Test$ExcludeFA)>0){
  thresh = paste0(100*round(max(Exclusion_Test$FA[Exclusion_Test$valid]),4),'%')
  invalids_subjects_text = paste0(invalids_subjects_text, 
                                 '\n# ',paste0(Exclusion_Test$subjectID_num[Exclusion_Test$ExcludeFA],collapse = ', '),' - FA (max = ',thresh,')')}
if (sum(Exclusion_Test$ExcludeMiss)>0){
  thresh = paste0(100*round(max(Exclusion_Test$Miss[Exclusion_Test$valid]),4),'%')
  invalids_subjects_text = paste0(invalids_subjects_text, 
                                 '\n# ',paste0(Exclusion_Test$subjectID_num[Exclusion_Test$ExcludeMiss],collapse = ', '),' - Miss (max = ',thresh,')')}
if (length(subjects_technical_issues)>0){
  invalids_subjects_text = paste0(invalids_subjects_text, 
                                  '\n# ',paste0(subjects_technical_issues,collapse = ', '),' - Technical issues')}

# print outcomes valid participants and outliers ----

cat("\nNum of participant to exclude based on transitivity: ",sum(Exclusion_Test$ExcludeTransitivityScore),
    '\n',paste(Exclusion_Test$Subject[Exclusion_Test$ExcludeTransitivityScore]),': ',round(Exclusion_Test$TransitivityScore[Exclusion_Test$ExcludeTransitivityScore],4),
    '\n',"Num of participant to exclude based on FA: ",sum(Exclusion_Test$ExcludeFA),
    '\n',paste(Exclusion_Test$Subject[Exclusion_Test$ExcludeFA]),': ',round(Exclusion_Test$FA[Exclusion_Test$ExcludeFA],4),
    '\n',"Num of participant to exclude based on Miss: ",sum(Exclusion_Test$ExcludeMiss),
    '\n',paste(Exclusion_Test$Subject[Exclusion_Test$ExcludeMiss]),': ',round(Exclusion_Test$Miss[Exclusion_Test$ExcludeMiss],4),'\n',
    sep=' ')
cat(paste0(valids_subjects_text,'\n\n',invalids_subjects_text),sep='')

Exclusion_Valid = Exclusion_Test %>% filter(valid==1) 
options(digits = 4)
cat("Transitivity score: M = ", mean(Exclusion_Test$TransitivityScore), ", 3SD cutoff = ", 
    mean(Exclusion_Test$TransitivityScore) - 3*sd(Exclusion_Test$TransitivityScore), ", mean valid score = ",
    mean(Exclusion_Valid$TransitivityScore), ", min valid score = ",min(Exclusion_Valid$TransitivityScore), ";\n",
    "False alarm rate: M = ", 100*mean(Exclusion_Test$FA), "%, 3SD cutoff = ", 
    100*(mean(Exclusion_Test$FA) + 3*sd(Exclusion_Test$FA)), "%, M valid score = ",
    100*mean(Exclusion_Valid$FA), "%, max valid score = ", 100*max(Exclusion_Valid$FA), "%;\n",
    "Miss rate: M = ", 100*mean(Exclusion_Test$Miss), "%, 3SD cutoff = ", 
    100*(mean(Exclusion_Test$Miss) + 3*sd(Exclusion_Test$Miss)), "%, M valid score = ",
    100*mean(Exclusion_Valid$Miss), "%, max valid score = ", 100*max(Exclusion_Valid$Miss),"%",
    sep = ""
)

