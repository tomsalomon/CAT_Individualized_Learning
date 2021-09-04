# Load libraries
library(stringr)
library(lme4)
library(lmerTest)
library(rstudioapi)

# Clear workspace
rm(list=ls())

# Define these variables:
sessionNum = 1 # 1- session 1 ; 2- follow-up
task_name = "Probe"
buffer_around_blink = 150 # time in ms to remove around blink 
screen=c(1920,1080) # size in pixels
buffer_around_image_factor = 1.2 # in relative size
fixcross_size = c(300,300) * buffer_around_image_factor; # size in pixels
stim_size = c(400, 500) * buffer_around_image_factor; # size in pixels
dist_from_center = 300 # distance of the stimuli's edge from the screen center
is_binocular = FALSE

# location information
center=screen/2
center_stim_l = c(center[1] - (dist_from_center + stim_size[1]/2), center[2])
center_stim_r = c(center[1] + (dist_from_center + stim_size[1]/2), center[2])
location_stim_l=as.data.frame(matrix(unlist(c(center_stim_l - stim_size/2, center_stim_l + stim_size/2)),1,4,T))
location_stim_r=as.data.frame(matrix(unlist(c(center_stim_r - stim_size/2, center_stim_r + stim_size/2)),1,4,T))
location_fixcross=as.data.frame(matrix(unlist(c(center - fixcross_size/2, center + fixcross_size/2)),1,4,T))
colnames(location_stim_l) = c('min_x', 'min_y', 'max_x', 'max_y')
colnames(location_stim_r) = colnames(location_stim_l)
colnames(location_fixcross) = colnames(location_stim_l)

# CAT allocations
probe_items = c(7:18,43:54)
training_items = c(3:22,39:58)
go_items_order1 = c(7,10,12,13,15,18,44,45,47,50,52,53)
go_items_order2 = setdiff(probe_items,go_items_order1)

# define which subjects to analyze
  subjects=c(112)
  scans = c(1)
ignore_bad_eye_tracking = c()
subjects = subjects[!(subjects %in% ignore_bad_eye_tracking)]

# set paths
current_path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_path)
data_path = './../Output/'
output_path = './analyzed_data/'

# Conversion from eyelink text file to data frame
eyetracker_text2df = function(x,filter='',delimiter = "\t") {
  if (filter == '') {
    dataframe_headers = c("time","x","y","pupil","V1","phase")
  } else if (filter == 'SSACC') {
    dataframe_headers = c("event","filler","empty","time")
    delimiter = "\ "
  } else if (filter == 'ESACC') {
    dataframe_headers = c("event","time","duration","startx","starty","endx","endy","amplitude","velocity")
    if (is_binocular) {
      dataframe_headers = c(dataframe_headers,"startx2","starty2","endx2","endy2","amplitude2","velocity2")
    }
  } else if (filter == 'SFIX') {
    dataframe_headers = c("event","filler","e1","e2","time")
    delimiter = "\ "
  } else if (filter == 'EFIX') {
    dataframe_headers = c("event","time","duration","x","y","pupil")
    if (is_binocular) {
      dataframe_headers = c(dataframe_headers,"x2","y2","pupil2")
    }
  } else if (filter == 'SBLINK') {
    dataframe_headers = c("event","filler","time")
    delimiter = "\ "
  } else if (filter == 'EBLINK') {
    dataframe_headers = c("event","time","duration")
  } 
  
  x_filtered = x[grep(filter,x)]
  # deal with non-blinker subject
  if (length(x_filtered) ==0) { 
    x_filtered = paste0(rep('\t ',length(dataframe_headers)))
  }
  df=structure(data.frame(matrix(unlist(strsplit(x_filtered,delimiter)),length(x_filtered),length(dataframe_headers),T)),  names=dataframe_headers)
  for (col_i in 1:length(dataframe_headers)) {
    is_numeric_col = mean(!is.na(suppressWarnings(as.numeric(as.character(df[,col_i])))))>0.05
    if (is_numeric_col) {
      df[,col_i]=suppressWarnings(as.numeric(as.character(df[,col_i])))
    }
  }
  return(df)
}

filelist = data.frame(file_ind=integer(), subjectID=integer(), scan=integer(), behavioral=character(), eyetracking=character(),stringsAsFactors=F)
ind=0
for (s in subjects){
  for (scan in scans){
    ind=ind+1
    filelist[ind,'file_ind'] = ind
    filelist$subjectID[ind] = s
    filelist$scan[ind] = scan
    filelist$behavioral[ind]  = paste0("",Sys.glob(paste0(data_path,"*",s,"*",task_name,"*.txt")))
    filelist$eyetracking[ind] = paste0("",Sys.glob(paste0(data_path,"*",s,"*",task_name,"*.asc")))  }
}
filelist$valid_eyetracker = !filelist$eyetracking==""
valid_files = filelist$file_ind[filelist$valid_eyetracker]
eyetrack_merged = c()
behave_merged = c()

# Loop over all files  ----------------------------------------
for (file_i in valid_files){
  progress =round(which(valid_files==file_i)/length(valid_files),4)*100
  print(paste0('Progress: ',progress,'%, File index: ',file_i))
  # read behavioral data
  behavioral=read.table(filelist$behavioral[file_i],header=T,na.strings=c(999,999000,"NaN"),sep="\t")
  order = behavioral$order[1];
  behavioral$scan = filelist$scan[file_i]
   # read eyetracking data
  input = readLines(filelist$eyetracking[file_i])
  # get the frequency
  freq_line =unlist(strsplit(input[grep('RATE',input)[1]],"\t"))
  freq = as.numeric(freq_line[1+grep('RATE',freq_line)])
  
  ssacc = eyetracker_text2df(input,'SSACC')
  esacc = eyetracker_text2df(input,'ESACC')
  sfix = eyetracker_text2df(input,'SFIX')
  efix = eyetracker_text2df(input,'EFIX')
  sblink = eyetracker_text2df(input,'SBLINK')
  eblink = eyetracker_text2df(input,'EBLINK')
  
  # on rare occasions eyelink 2 SFIX at the beginning of the file
  if (nrow(sfix) - nrow(efix) == 1){
    if (sfix$time[2] < efix$time[1]){
      sfix=sfix[-2,]
    }
  }
  # skip files with missmatch between events start and end - these are usualy empty or corrupted files
  if (!(
    (nrow(ssacc) == nrow(esacc)) & (nrow(sfix) == nrow(efix)) & (nrow(sblink) == nrow(eblink)))){
    cat('\nWARNING! number of events do not match. skipping!\n',
        'you should look at this file:',filelist$eyetracking[file_i],'\n')
    next
  }
  
  #don't want header, figure out the first line that is not header
  for (line_i in 1:2000){
    if (!is.na(suppressWarnings(as.numeric(strsplit(input[line_i],"\t")[[1]][1])))){
      firstline=line_i
      break
    }
  }
  
  exclude=c(1:firstline,grep("MSG|SFIX|EFIX|SSACC|ESAC|SBLINK|EBLINK|SAMPLES|INPUT",input))
  eyetracking_data = eyetracker_text2df(input[-exclude])
  eyetracking_data$subjectID = behavioral$subjectID[1]
  eyetracking_data$session = sessionNum[1]
  eyetracking_data$scan = scan[1]
  eyetracking_data$trial = NA
  eyetracking_data$phase = NA
  eyetracking_data$event = 'fixation' #assume all non blink or saccade is fixation
  eyetracking_data$eventnumber = NA
  eyetracking_data$valid = (!is.na(eyetracking_data$x)) & (!is.na(eyetracking_data$y))
  
  # Analyze general eyetracking events: fixations, blinks and saccades =======================================
  eyetracking_events = data.frame(start = c(sfix$time,ssacc$time,sblink$time),end = c(efix$time,esacc$time,eblink$time))
  eyetracking_events$event = c(rep('fixation',nrow(sfix)), rep('saccade',nrow(ssacc)), rep('blink',nrow(sblink)))
  eyetracking_events$eventnumber = c(1:nrow(sfix),1:nrow(ssacc), 1:nrow(sblink))
  for (i in 1:length(eyetracking_events$event)){
    s1 = which.min(abs(eyetracking_data$time - eyetracking_events$start[i]))
    e1 = which.min(abs(eyetracking_data$time - eyetracking_events$end[i]))
    eyetracking_data$event[s1:e1]=eyetracking_events$event[i]
    eyetracking_data$eventnumber[s1:e1]=eyetracking_events$eventnumber[i]
    
    #Account for blinks, get rid of 150ms before and after blink
    if (eyetracking_events$event[i] == 'blink') {
      s2 = which.min(abs(eyetracking_data$time - (eyetracking_events$start[i] - buffer_around_blink)))
      e2 = which.min(abs(eyetracking_data$time - (eyetracking_events$end[i] + buffer_around_blink)))
      eyetracking_data$valid[s2:e2]=FALSE
    }
  }
  
  eyetracking_data$looking_at_stim_l = (eyetracking_data$x >= location_stim_l$min_x) & (eyetracking_data$y>=location_stim_l$min_y) & (eyetracking_data$x<=location_stim_l$max_x) & (eyetracking_data$y<=location_stim_l$max_y)
  eyetracking_data$looking_at_stim_l[(is.na(eyetracking_data$looking_at_stim_l) | !eyetracking_data$valid)] = FALSE
  eyetracking_data$looking_at_stim_r = (eyetracking_data$x >= location_stim_r$min_x) & (eyetracking_data$y>=location_stim_r$min_y) & (eyetracking_data$x<=location_stim_r$max_x) & (eyetracking_data$y<=location_stim_r$max_y)
  eyetracking_data$looking_at_stim_r[(is.na(eyetracking_data$looking_at_stim_r) | !eyetracking_data$valid)] = FALSE
  eyetracking_data$looking_at_fixation = (eyetracking_data$x >= location_fixcross$min_x) & (eyetracking_data$y>=location_fixcross$min_y) & (eyetracking_data$x<=location_fixcross$max_x) & (eyetracking_data$y<=location_fixcross$max_y)
  eyetracking_data$looking_at_fixation[(is.na(eyetracking_data$looking_at_fixation) | !eyetracking_data$valid)] = FALSE
  
  #eyetracking_data$looking_at_fixation = (eyetracking_data$x >= location_fixcross$min_x) & (eyetracking_data$y>=location_fixcross$min_y) & (eyetracking_data$x<=location_fixcross$max_x) & (eyetracking_data$y<=location_fixcross$max_y)
  
  # Analyze task behavioral data =======================================
  behavioral$freq = freq
  run_start_line=unlist(strsplit(input[grep("SYNCTIME",input)],"\t| "))
  run_start_time=as.numeric(run_start_line[2])
  behavioral$duration[is.na(behavioral$duration)] = 1.5
  trial_start_times=round(run_start_time+(behavioral$onset*1000))
  trial_end_times=round(trial_start_times+(behavioral$duration*1000))
  numtrials = length(behavioral$onset)
  # classify lines in eyetracking_data to choice/confirmation/fixation
  for (i in 1:numtrials){
    s1 = which.min(abs(eyetracking_data$time - trial_start_times[i]))
    e1 = which.min(abs(eyetracking_data$time - trial_end_times[i]))
    s2 = e1 + 1
    if (i<numtrials) {
      e2 = which.min(abs(eyetracking_data$time - trial_start_times[i+1])) - 1
    } else {
      e2 = length(eyetracking_data$time)
    }
    eyetracking_data$phase[s1:e1]="trial"
    eyetracking_data$phase[s2:e2]="ISI"
    eyetracking_data$trial[s1:e2]=i
    behavioral$gaze_stim_l [i] = mean(eyetracking_data$looking_at_stim_l[s1:e1])
    behavioral$gaze_stim_r [i] = mean(eyetracking_data$looking_at_stim_r[s1:e1])
    behavioral$fixation [i] = mean(eyetracking_data$looking_at_fixation[s1:e1])
    behavioral$gaze_ISI [i] = mean(eyetracking_data$looking_at_fixation[s2:e2])
    behavioral$time_duration [i] = eyetracking_data$time[e1] - eyetracking_data$time[s1]
  }

  eyetrack_tmp = eyetracking_data[eyetracking_data$phase %in% "trial",]
  eyetrack_merged = rbind(eyetrack_merged,eyetrack_tmp)
  behave_merged = rbind(behave_merged,behavioral)
  # output_filename = paste0(output_path,tail(unlist(strsplit( filelist$behavioral[file_i],'/')),n=1))
  # write.table(eyetracking_data,file=output_filename,sep = '\t',row.names = F)
}
behave_merged$scan = as.factor(behave_merged$scan)
behave_merged$PairType2 = as.factor(behave_merged$PairType)
behave_merged$PairType2 = NA
behave_merged$PairType2[behave_merged$PairType %in% c(1,2)] = behave_merged$PairType[behave_merged$PairType %in% c(1,2)]
behave_merged$PairType2 = factor(behave_merged$PairType2, labels = c('High-Value', 'Low-Value'))
behave_merged$gaze_general = behave_merged$gaze_stim_l + behave_merged$gaze_stim_r + behave_merged$fixation
behave_merged$gaze_go = NA
behave_merged$gaze_nogo = NA
behave_merged$gaze_go[behave_merged$IsleftGo == 1] = behave_merged$gaze_stim_l [behave_merged$IsleftGo == 1]
behave_merged$gaze_nogo[behave_merged$IsleftGo == 1] = behave_merged$gaze_stim_r [behave_merged$IsleftGo == 1]
behave_merged$gaze_go[behave_merged$IsleftGo == 0] = behave_merged$gaze_stim_r [behave_merged$IsleftGo == 0]
behave_merged$gaze_nogo[behave_merged$IsleftGo == 0] = behave_merged$gaze_stim_l [behave_merged$IsleftGo == 0]
behave_merged$gaze_delta = behave_merged$gaze_go - behave_merged$gaze_nogo
print ('Completed! Saving all data into R data frame')
output_filename = paste0(output_path,'Eyetracking_Task_',task_name,'_Session',sessionNum,'.Rda')
save(eyetrack_merged,file=output_filename)
output_filename = paste0(output_path,'Summary_Task_',task_name,'_Session',sessionNum,'.Rda')
save(behave_merged,file=output_filename)
print (paste0('Done. You can find the data at: ',output_filename))
