function [] = organize_CAT(subjectID,outputPath)
% Created by: Tom Salomon, July 2019.
%
% This function sorts the stimuli according to the subjective preference
%% define variables 
if ~exist('outputPath','var')
    outputPath = './Output';
end
subjectID_num = str2double(subjectID(end-2:end));
order = rem(subjectID_num,4);
if order ==0
    order = 4;
end

n_runs = 20; % number of training repetitions
n_runs_full_association = 2; % number of training repetitions where all stimuli are associated
n_groups = 8; % number of value group
n_Go_per_group = 4; % number of Go and NoGo stimuli per group
n_sanity = 2;
n_probe_blocks = 2; % to how many block split the probe
n_probe_runs = 1; % how many repetitions of each probe choice - CURRENTLY NOT APPLIED CHANGE CODE TO HAVE MORE THAN 1 REPETITION

%% define Go/Nogo and contingency design
% Go allocation - which items will be Go and which NoGo
% Hard coded for 4 Go versus NoGo stimuli - Change if use different numbers
Go_allocation_opt = [
    [1;0;0;1;1;0;0;1] ,...
    [1;0;0;1;0;1;1;0] ,...    
    [0;1;1;0;0;1;1;0] ,...
    [0;1;1;0;1;0;0;1] ];
% Value categories contingency allocation - which value category will be
% associated with which contingency
% Hard coded for 8 value categories with 100%/50% contingency  - Change if use different numbers
contingency_opt = [
    [1 ;.5;.5;1 ;1 ;.5;.5;1] ,...
    [1 ;.5;.5;1 ;.5;1 ;1 ;.5] ,...    
    [.5;1 ;1 ;.5;.5;1 ;1 ;.5] ,...
    [.5;1 ;1 ;.5;1 ;.5;.5;1 ] ];

%% read data file
subj_rank_file=dir(fullfile(outputPath,[subjectID,'_ItemRankingResults*']));
subj_rank_data = readtable([outputPath,'/',subj_rank_file.name],'delimiter','\t');
% sort descending by value
subj_rank_data = sortrows(subj_rank_data,'Rank','descend');
n_stim = length(subj_rank_data.Rank); % 80
n_non_filler = n_groups*n_Go_per_group*2 ; % 8*4*2 + 4 - 64 : 32 Go 32 NoGo
n_filler = n_stim - n_non_filler; % 12

fillers = [1:(n_filler/2), (1 + n_stim - n_filler/2):n_stim]; % edges of value distribution used as fillers. 
value_group_tmp = repmat(1:n_groups,[n_Go_per_group*2,1]);
value_group_tmp = value_group_tmp(:);
value_group = zeros(n_stim,1);
value_group((n_filler/2+1) : (n_stim - n_filler/2)) = value_group_tmp;

%% Create table with Go/NoGo allocation

Go_contingency = nan(n_stim,1);
% Value-Contingency association is determined according to participant's code
contingency_per_category = contingency_opt(:,order); 
for group_i = 1:n_groups
    rnd_selection = randi(size(Go_allocation_opt,2)); % randomly select Go allocation order
    loc = value_group == group_i;
    Go_contingency(loc,1) = Go_allocation_opt(:,rnd_selection)*contingency_per_category(group_i);
end
Go_contingency(fillers,1) = 0;

training_table = subj_rank_data(:,1); % Subject ID
training_table.StimRank = (1:n_stim)'; % ordered-rank
training_table.StimValue = subj_rank_data.Rank; % value (Bid/Colley rank)
training_table.StimName = subj_rank_data.StimName; % stim name 
training_table.StimNum = subj_rank_data.StimNum; % stim code (ordered alphabetically) 
training_table.ValueGroup = value_group; % value group
training_table.Contingency = Go_contingency; % Go Contingency

% build a matrix with whether each stimulus will be Go in each training run
% last (2) runs associated for all stimuli with p>0, to be scanned in fMRI
training_asso = zeros(n_stim,n_runs);
for stim_i = 1:n_stim
    p = training_table.Contingency(stim_i);
    if p>0
        n_assoc = p*n_runs - n_runs_full_association;
        
        training_asso(stim_i,1+end-n_runs_full_association:end) = 1;
        training_asso(stim_i,1:end-n_runs_full_association) = Shuffle([ones(1,n_assoc),zeros(1,n_runs - n_assoc - n_runs_full_association)]);
    end
end

training_table_full = training_table([],:);
training_table_full.Go = zeros(0,1);
training_table_full.Run = zeros(0,1);
training_table_full.TrialInRun = zeros(0,1);

for run = 1:n_runs
   training_table_tmp = training_table;
   training_table_tmp.Go = training_asso(:,run);
   training_table_tmp.Run = repmat(run,[n_stim,1]);
   training_table_tmp.TrialInRun = randperm(n_stim)';
   training_table_tmp = sortrows(training_table_tmp,'TrialInRun','ascend');
   training_table_full = [training_table_full;training_table_tmp];
end
training_table_full.Trial = (1:size(training_table_full,1))';

%% Probe design
stim_Go_ind = repmat(1:n_Go_per_group,[n_Go_per_group,1]);
stim_Go_ind = stim_Go_ind(:);
stim_NoGo_ind = repmat(1:n_Go_per_group,[1,n_Go_per_group]);
stim_NoGo_ind = stim_NoGo_ind(:);

probe_Go_stim = training_table([],:);
probe_NoGo_stim = training_table([],:);
probe_Go_stim.block = ones(0,1);
for group_i = 1:n_groups
    % extract Go / NoGo stim info
    Go_stim_tmp = training_table(training_table.ValueGroup == group_i & training_table.Contingency > 0,:);
    NoGo_stim_tmp = training_table(training_table.ValueGroup == group_i & training_table.Contingency == 0,:);
    % duplicate info to create all possible combinations
    Go_stim_tmp2 = Go_stim_tmp(stim_Go_ind,:);
    NoGo_stim_tmp2 = NoGo_stim_tmp(stim_NoGo_ind,:);
    
    Go_stim_tmp2.block = Shuffle(repmat((1:n_probe_blocks)',[size(Go_stim_tmp2,1)/n_probe_blocks,1]));
    probe_Go_stim = [probe_Go_stim;Go_stim_tmp2];
    probe_NoGo_stim = [probe_NoGo_stim;NoGo_stim_tmp2];
end
n_probe_trials = size(probe_Go_stim,1);

% Build probe table
probe_table_headers = {'subjectID','scanner','order','block','run','trial','onsettime','ImageLeft','ImageRight','bidIndexLeft','bidIndexRight','IsleftGo','Response','PairType','Outcome','RT','bidLeft','bidRight','FixationTime','Contingency'};
probe_table_full = array2table(nan(0,length(probe_table_headers)),'VariableNames',probe_table_headers);

for run = 1:n_probe_runs
    probe_table = array2table(nan(n_probe_trials,length(probe_table_headers)),'VariableNames',probe_table_headers);
    % Shuffle left and right
    is_left_go = Shuffle([ones(n_probe_trials/2,1);zeros(n_probe_trials/2,1)]);
    probe_left_stim = (probe_Go_stim.StimRank).*is_left_go + (probe_NoGo_stim.StimRank).*(1-is_left_go);
    probe_right_stim = (probe_NoGo_stim.StimRank).*is_left_go + (probe_Go_stim.StimRank).*(1-is_left_go);
    
    probe_table.subjectID = repmat(subjectID,[n_probe_trials,1]);
    probe_table.order = repmat(order,[n_probe_trials,1]);
    probe_table.block = probe_Go_stim.block;
    probe_table.run = repmat(run,[n_probe_trials,1]); % probe run currently not applied. modify code to have more than 1 repetition
    probe_table.ImageLeft = training_table.StimName(probe_left_stim);
    probe_table.ImageRight = training_table.StimName(probe_right_stim);
    probe_table.bidIndexLeft =  training_table.StimRank(probe_left_stim);
    probe_table.bidIndexRight =  training_table.StimRank(probe_right_stim);
    probe_table.IsleftGo = is_left_go;
    probe_table.PairType = probe_Go_stim.ValueGroup;
    probe_table.bidLeft = training_table.StimValue(probe_left_stim);
    probe_table.bidRight = training_table.StimValue(probe_right_stim);
    probe_table.FixationTime = nan(n_probe_trials,1);
    probe_table.Contingency = training_table.Contingency(probe_left_stim) + training_table.Contingency(probe_right_stim);
    
    n_probe_trials_per_block = n_probe_trials/n_probe_blocks;
    probe_table = sortrows(probe_table,'block','ascend');
        
    probe_table.trial = (probe_table.run-1)*n_probe_trials +  (probe_table.block-1)*n_probe_trials_per_block + repmat(randperm(n_probe_trials_per_block)',[n_probe_blocks,1]);
   
    probe_table = sortrows(probe_table,'trial','ascend');
    probe_table_full = [probe_table_full;probe_table];
end

%% save files
c = clock;
timestamp = sprintf('%04i%02i%02i_%02dh%02dm',c(1),c(2),c(3),c(4),c(5));
training_file_name = sprintf('%s/%s_design_CAT_%s',outputPath,subjectID,timestamp);
stimlu_file_name = sprintf('%s/%s_design_stimuli_%s',outputPath,subjectID,timestamp);
probe_file_name = sprintf('%s/%s_design_probe_%s',outputPath,subjectID,timestamp);

writetable(training_table_full,training_file_name,'Delimiter','\t');
writetable(training_table,stimlu_file_name,'Delimiter','\t');
writetable(probe_table_full,probe_file_name,'Delimiter','\t');
