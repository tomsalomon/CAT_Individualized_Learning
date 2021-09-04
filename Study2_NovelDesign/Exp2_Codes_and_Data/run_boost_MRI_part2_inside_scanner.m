%% MRI Script
clear all;

%% Do the screen resolution magic thingy

%% Response to stimuli - Demo
responseToStimuliDemo;

%% Response to stimuli - 2 runs
responseToStimuli_short;

%% Training - Demo
trainingDemo_Israel;

%% Training - 5 runs
training_Israel;

%% MPRAGE + FLAIR

%% Response to stimuli - 2nd scan - 2 runs
responseToStimuli_short;

%% Anatomical scans time

%% organize probe
subjectID = input('Subject code: ','s');
[subjectID_num,okID]=str2num(subjectID(end-2:end));
while okID==0
    disp('ERROR: Subject code must contain 3 characters numeric ending, e.g "BMI_bf_101". Please try again.');
    subjectID = input('Subject code:','s');
    [subjectID_num,okID]=str2num(subjectID(end-2:end));
end

% Assign order
% --------------------------
% give order value of '1' or '2' for subjects with odd or even ID, respectively
if mod(subjectID_num,2) == 1 % subject code is odd
    order = 1;
else % subject code is even
    order = 2;
end

sessionNum=1;

numBlocks = 2;
numRunsPerBlock = 2;
for block = 1:numBlocks
    organizeProbe_Israel(subjectID, order, pwd, block, numRunsPerBlock);
end
trialsPerRun=38;

%% Probe - Demo
probeDemo_Israel;

%% Probe - 4 runs (2 blocks, 2 runs each)
probe_Israel;

%% Face Localizer - 2 runs
cd './FaceObjectLocalizer';
dynamicLocalizer;
cd './../';

%% End of fMRI scans



