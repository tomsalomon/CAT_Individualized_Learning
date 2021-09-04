
library(lme4)
library("rstudioapi")    

rm(list=ls())

# Get current path
script_path = rstudioapi::getActiveDocumentContext()$path
pathSplit=strsplit(script_path, "/")
pathSplit=pathSplit[[1]]
main_path=paste0(pathSplit[1:(length(pathSplit)-2)],"/",collapse="")

## Sample
path=paste0(main_path,"/Output/")
subjects=c(101:110,112:121); # 20 valid subjects

# Not really excluded:
# 111 - Eyetracker crashed code during the experiment. 

filelist=c()
for (s in subjects){
  filelist=c(filelist,Sys.glob(paste(path, "BM_",s,"_CAT_*.txt",sep="")))
}
CAT_data=c()
for (f in filelist){
  CAT_data=rbind(CAT_data,read.table(f,header=T,na.strings=c(999,999000)))
}

filelist=c()
for (s in subjects){
  filelist=c(filelist,Sys.glob(paste(path, "BM_",s,"_ItemRanking*.txt",sep="")))
}
RankingData = c()
for (f in filelist){
  RankingData=rbind(RankingData,read.table(f,header=T,na.strings=c(999,999000)))
}

# Transitivity Scores
Exclusion_Test = as.data.frame(tapply(RankingData$Rank, RankingData$Subject, sd))
colnames(Exclusion_Test) = "TransitivityScore"
Exclusion_Test$subjectID = rownames(Exclusion_Test)


CAT_data$Response = 1 - is.na(CAT_data$RT)
CAT_data_GO = subset(CAT_data, Go==1)
CAT_data_NoGO = subset(CAT_data, Go==0)

TP = with(data = subset(CAT_data, Go==1), tapply(Response, Subject, mean))
TN = 1 - with(data = subset(CAT_data, Go==0), tapply(Response, Subject, mean))
FA = with(data = subset(CAT_data, Go==0), tapply(Response, Subject, mean))
Miss = 1 - with(data = subset(CAT_data, Go==1), tapply(Response, Subject, mean))

Exclusion_Test$FA = FA
Exclusion_Test$Miss = Miss

Exclusion_Test$ZTransitivityScore = scale(Exclusion_Test$TransitivityScore)
Exclusion_Test$ZFA = scale(FA)
Exclusion_Test$ZMiss = scale(Miss)

Exclusion_Test$ExcludeTransitivityScore = Exclusion_Test$ZTransitivityScore <= -3
Exclusion_Test$ExcludeFA = Exclusion_Test$ZFA >= 3
Exclusion_Test$ExcludeMiss = Exclusion_Test$ZMiss >= 3

cat("Num of participant to exclude based on transitivity: ",sum(Exclusion_Test$ExcludeTransitivityScore),'\n',paste(Exclusion_Test$subjectID[Exclusion_Test$ExcludeTransitivityScore],':',Exclusion_Test$TransitivityScore[Exclusion_Test$ExcludeTransitivityScore],'\n'))
cat("Num of participant to exclude based on FA: ",sum(Exclusion_Test$ExcludeFA),'\n',paste(Exclusion_Test$subjectID[Exclusion_Test$ExcludeFA],':',Exclusion_Test$FA[Exclusion_Test$ExcludeFA],'\n'))
cat("Num of participant to exclude based on Miss: ",sum(Exclusion_Test$ExcludeMiss),'\n',paste(Exclusion_Test$subjectID[Exclusion_Test$ExcludeMiss],':',Exclusion_Test$Miss[Exclusion_Test$ExcludeMiss],'\n'))

options(digits = 4)
cat("Transitivity score: M = ", mean(Exclusion_Test$TransitivityScore), ", 3SD cutoff = ", 
    mean(Exclusion_Test$TransitivityScore) - 3*sd(Exclusion_Test$TransitivityScore), ", min valid score = ",
    min(Exclusion_Test$TransitivityScore), ";\n",
    "False alarm rate: M = ", 100*mean(Exclusion_Test$FA), "%, 3SD cutoff = ", 
    100*(mean(Exclusion_Test$FA) + 3*sd(Exclusion_Test$FA)), "%, max valid score = ",
    100*max(Exclusion_Test$FA), "%;\n",
    "Miss rate: M = ", 100*mean(Exclusion_Test$Miss), "%, 3SD cutoff = ", 
    100*(mean(Exclusion_Test$Miss) + 3*sd(Exclusion_Test$Miss)), "%, max valid score = ",
    100*max(Exclusion_Test$Miss), "%",
    sep = ""
    )
