% =========================================================================
%% Get input args and check if input is ok
% =========================================================================
close all;
clear;

% variables to define
exp_code = 'BM';
sessionNum=1;
is_MRI=0; % 1 MRI, 0 if testooom
mainPath = '.'; % Change if you don't run from the experimnet folder - not recomended.
outputPath = [mainPath '/Output'];
subjectID = GetSubjectID(exp_code);

% Increase the size of the UI fonts
old_UI_font_size=get(0,'defaultUicontrolFontSize');
set(0,'defaultUicontrolFontSize', 14);

% ask user whether to use eye tracker
ask_eyetracker = questdlg('Do you want to use the eye tracker?','Use eye tracker','Yes','No','Yes');
use_eyetracker = strcmp(ask_eyetracker,'Yes');

% essential for randomization
rng('shuffle');

% Assign order: 1-4
subjectID_num = str2double(subjectID(end-2:end));
order = rem(subjectID_num,4);
if order ==0
    order = 4;
end

% =========================================================================
%% Select Starting Point
% =========================================================================
StartingPoint=1; % start at the very first part
% Check for pre-existing data
ExistingSubjectData=dir([outputPath,'/',subjectID,'*.txt']);
if ~isempty(ExistingSubjectData)
    DifferentStartingPoint = questdlg(...
        sprintf(['I noticed this subject code (%s) have some pre-existing data.\n'...
        'Do you want to start the experiment from the beggining?'],subjectID)...
        ,'WARNING - subject code already exist','Yes','No','No');
    if strcmp(DifferentStartingPoint,'No')
        task_parts = {
            '1: Beginning';
            '2: Binary choices I';
            '3: CAT';
            '4: Binary choices II';
            '';
            'Oops. Stop the Experiment'};
        StartingPoint = listdlg('PromptString','Select the part you want to start from:',...
            'SelectionMode','single','ListString',task_parts,'ListSize',[300,300]);
    else
        StartingPoint=1;
    end
end

% =========================================================================
%% Task 1 - Personal Details
% =========================================================================
if StartingPoint<=1
    personal_details(subjectID, order, outputPath, sessionNum)
end

% =========================================================================
%% Task 2 - Binary Ranking (including demo)
% =========================================================================
if StartingPoint<=2
    % demo
    demo = true;
    while demo
        binary_ranking(subjectID,demo,use_eyetracker);
        repeat_demo = questdlg('Do you want to run the demo again?','Repeat Demo','Yes','No','No');
        demo = strcmp(repeat_demo,'Yes');
    end
    % full task
    binary_ranking(subjectID,demo,use_eyetracker);
    
    %Sort stimuli according to the binary ranking
    % -------------------------------------------------
    organize_CAT(subjectID,outputPath);
    questdlg('Great! Please call the experimenter to continue to the next part','Next part','Continue','Continue');
end

% =========================================================================
%% Task 3 - CAT
% =========================================================================
if StartingPoint<=3
    % demo
    demo = true;
    while demo
        CAT(subjectID,demo,use_eyetracker,is_MRI);
        repeat_demo = questdlg('Do you want to run the demo again?','Repeat Demo','Yes','No','No');
        demo = strcmp(repeat_demo,'Yes');
    end
    % full task
    CAT(subjectID,demo,use_eyetracker,is_MRI);
    questdlg('Great! Please call the experimenter to continue to the next part','Next part','Continue','Continue');
end

% =========================================================================
%% Task 4 - Probe
% =========================================================================
if StartingPoint<=4
    % demo
    demo = true;
    while demo
        probe(subjectID,demo,use_eyetracker,is_MRI);
        repeat_demo = questdlg('Do you want to run the demo again?','Repeat Demo','Yes','No','No');
        demo = strcmp(repeat_demo,'Yes');
    end
    % full task
    probe(subjectID,demo,use_eyetracker,is_MRI);
end

% =========================================================================
%% End of experiment
% =========================================================================
% Copy the data to cloud location
sync_with_cloud(subjectID,1)

questdlg('Thank you! The experiment is done.','End of Experiment','Continue','Continue');

% Return to the defualt UI font size
set(0,'defaultUicontrolFontSize', old_UI_font_size)
