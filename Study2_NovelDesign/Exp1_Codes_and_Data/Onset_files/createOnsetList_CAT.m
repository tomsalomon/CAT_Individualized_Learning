%function [averageDiff, maximumDiff, minimumDiff, diff] = createOnsetList()
% -------------------------------------
% Edited by Tom Salomon, July 2019. Based on code from Rotem Bitvinik-Nezer
% -------------------------------------

% function [averageDiff, maximumDiff, minimumDiff, diff] = createOnsetList(mean_t,min_t,max_t,interval,requiredLength,numTrials,numOnsetlists)
% Creates and saves numOnsetlists onsetlists for the probe of the cue approach task.
% This function makes sure that the mean of diff (variable averageDiff) is close enough (+- 2%) to
% the requested mean value and that the length of the onsetlist (that is,
% the last onset time) is the same in all the lists (as requested in the
% requiredLength variable).

clear all
close all

%% Define Variables
task_name = 'CAT';
trial_duration = 1; % duration of non ISI trial: training = 1; probe = 2;
mean_t = 2; % the mean requested ISI jitter (not including the trial duration)
min_t = 1; % min_t: the minimum requested jitter (not including the 2 seconds of stimuli presentation)
max_t = 10; % max_t: the maximum requested jitter (not including the 2 seconds of stimuli presentation)
interval = 0.1; % interval: the interval to ceil to (1 for integers)
numTrials = 80; % number trials: training = 40/60/80
numOnsetlists = 4; % how many onsetlists should the function create

% ----------------------------------
%% Creating onsetlists
% ----------------------------------
requiredLength = (mean_t + trial_duration)*(numTrials-1);

for list = 1:numOnsetlists
differences = 0;
%while (sum(diff) > requiredLength) || (abs(sum(diff)-requiredLength)> max_t/2)
while sum(differences) ~= requiredLength % can only work if using large enough interval
    tmp=exprnd(mean_t*0.7,[numTrials*10,1]); % random exponential
    tmp=tmp-mod(tmp,interval); % round to nearest interval
    tmp = tmp(tmp<=max_t & tmp>=min_t); % truncate the exp to desired limits
    tmp = Shuffle(tmp);
    differences = trial_duration + tmp(1:(numTrials-1));
end

onsets = [0;cumsum(differences)];
save(sprintf('%s_onsets_%i',task_name,list),'onsets');
fprintf('complete list %i out of %i\n',list,numOnsetlists);
end

