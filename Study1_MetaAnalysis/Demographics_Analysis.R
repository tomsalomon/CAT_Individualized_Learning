
# Load rquired R packages
library(rstudioapi) # to define current working directory
library(lme4)
library(dplyr)
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


demographics = c()

for (exp_num in 1:length(experiment_paths)){
  exp_path = experiment_paths[exp_num]
  exp_name = substr(exp_path,nchar(data_path_main),nchar(exp_path))
  
  print("")
  print(paste0("Experiment ",exp_num,": ",exp_name))
  print("==============================")
  
  filelist=Sys.glob(paste0(exp_path, "/*personal*.txt"))
  
  if (length(filelist)>0){
    # read data
    tmp_data=c()
    for (f in filelist){
      tmp_data_i = read.delim(f,row.names = NULL, nrows=1, fill = TRUE, na.strings = "")
      if (colnames(tmp_data_i)[1]== "row.names"){
        colnames(tmp_data_i) = c(colnames(tmp_data_i)[2:(ncol(tmp_data_i)-1)],"NA")
      }
      subjectID = tmp_data_i[1,grepl(pattern = "sub", ignore.case = TRUE, x = colnames(tmp_data_i))]
      gender = tmp_data_i[1,grepl(pattern = "gender", ignore.case = TRUE, x = colnames(tmp_data_i))]
      age = tmp_data_i[1,grepl(pattern = "age", ignore.case = TRUE, x = colnames(tmp_data_i))]
      
      tmp_data=rbind(tmp_data,data.frame(subjectID = subjectID, gender = gender, age=age))
    }
    tmp_data
    tmp_data$Experiment = exp_num
    tmp_data$ExperimentName = exp_name
    demographics=rbind(demographics,tmp_data)
    
  }
}

demographics$unpublished = grepl(pattern = 'UnPub', x= demographics$ExperimentName)

demographics %>% filter(unpublished == TRUE) %>% group_by(Experiment, ExperimentName) %>% summarise(N = n(), 
                                                                    Gender = sum(gender==1, na.rm=T), GenderProp = mean(gender==1,na.rm=T),
                                                                    Age_mean = mean(age, na.rm = T), Age_SD = sd(age, na.rm = T), Age_min =min(age, na.rm = T),Age_max =max(age, na.rm = T) )

