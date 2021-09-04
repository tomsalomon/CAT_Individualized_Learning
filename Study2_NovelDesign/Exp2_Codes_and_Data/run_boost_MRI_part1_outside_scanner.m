% function run_boost_Israel()

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ======================= By Tom Salomon, July 2019 =======================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% Runs the CAT: Part1 - tasks outside the scanner:
% Subjective preference evaluation (Binary ranking)
% Training (18 runs)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FUNCTIONS REQUIRED TO RUN PROPERLY: ----------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   --- Cue-approach codes: ---
% % %   'binary_ranking'
% % %   'sort_binary_ranking'
% % %   'CAT'
% % %   'organizeProbe_Israel'
% % %   'probe'

% % %   --- Other codes: ---
% % %  'CenterText'

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FOLDERS REQUIRED TO RUN PROPERLY: ------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   'Misc': a folder with the audio file.
% % %   'Onset_files': a folder with the onset files for the training and
% % %    for the probe.
% % %   'Output': the folder for the output files- results.
% % %   'stim': with the image files of all the stimuli for the cue-approach
% % %    task (old stimuli).
% % %   'stim/recognitionNew': with the image files of the new stimuli
% % %   (stimuli that are not included in the cue-approach tasks, only in the
% % %   recognitionNewOld task, as new stimuli)

% =========================================================================
%% Get input args and check if input is ok
% =========================================================================
close all;
clear;

% variables to define
exp_code = 'BM';
sessionNum=1;
test_comp=0; % 1 MRI, 0 if testooom
mainPath = '.'; % Change if you don't run from the experimnet folder - not recomended.
outputPath = [mainPath '/Output'];
subjectID = GetSubjectID(exp_code);
order = 1;

% ask user whether to use eye tracker
ask_eyetracker = questdlg('Do you want to use the eye tracker?','Use eye tracker','Yes','No','Yes');
use_eyetracker = strcmp(ask_eyetracker,'Yes');

% timestamp
c = clock;
hr = sprintf('%02d', c(4));
min = sprintf('%02d', c(5));
timestamp=[date,'_',hr,'h',min,'m'];

% essential for randomization
rng('shuffle');

%% Personal Details
% =========================================================================
personal_details(subjectID, order, outputPath, sessionNum)

%% Task 1 - Binary Ranking (including demo)
% demo
demo = true;
while demo
        binary_ranking(subjectID,demo,use_eyetracker);
        repeat_demo = questdlg('Do you want to run the demo again?','Repeat Demo','Yes','No','No');
        demo = strcmp(repeat_demo,'Yes');
end

% full task
binary_ranking(subjectID,demo,use_eyetracker);

%% Sort stimuli according to the binary ranking
% -------------------------------------------------
sort_binary_ranking(subjectID,order,outputPath);

%% CAT - outside scanner
CAT(subjectID);

