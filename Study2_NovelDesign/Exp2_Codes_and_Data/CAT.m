

function CAT(subjectID,demo,use_eyetracker,is_MRI,MRI_scan)
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ==================== by Tom Salomon, July 2019 =====================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% This function runs the boost (cue-approach) training session,
% in which the items are shown on the screen while some of them (GO items) are
% paired with a cue. The subject should press a predefined button as fast
% as possible after hearing the beep.
% This session is composed of total_num_runs_training number of runs.
% After two runs there is a short break.
% =========================================================================
% Get input args and check if input is ok
% =========================================================================

%==============================================
%% INPUTS
%==============================================
% demo
if (~exist('demo','var'))|| isempty(demo)
    demo_ans = questdlg('Do you want to run to full task or demo?','Ask Demo',...
        'Full Task','Demo','Full Task');
    demo = strcmp(demo_ans,'Demo'); % set to 1/0 to run demo/full task
end

% SubjectID
exp_code = 'BM';
if (~exist('subjectID','var')) || isempty(subjectID)
    if (~demo)
        suggested_subID = whats_the_next_sub('read'); % get suggestion for the next participant code
        okID = false;
        while ~okID % make sure code is 3 digit numeric
            subjectID_num = str2double(inputdlg('Please insert subject code (3 digits only). 999 for debugging',...
                'subject code',1,{num2str(suggested_subID)}));
            if isnumeric(subjectID_num) && (~isempty(subjectID_num))% input checkers
                okID = subjectID_num<=999;
            else
                okID = false;
            end
        end
    elseif demo
        subjectID_num = 999;
    end
    subjectID = sprintf('%s_%03i',exp_code, subjectID_num);
end
% Ask if want eye tracker
if ((~exist('use_eyetracker','var')) || isempty(use_eyetracker) )&& (~demo)
    ask_eyetracker = questdlg('Do you want to use the eye tracker?','Use eye tracker','Yes','No','Yes');
    use_eyetracker = strcmp(ask_eyetracker,'Yes');
elseif demo
    use_eyetracker=0;
end
% Ask if this is an MRI environment
if (~exist('is_MRI','var')) || isempty(is_MRI)
    ask_MRI = questdlg('Are you running in the MRI?','MRI environment','Yes','No','Yes');
    is_MRI = strcmp(ask_MRI,'Yes'); % 1 MRI, 0 if testooom
end
% Ask which MRI scan is it (and define defualr runs used for non MRI
% environment and demo)

if is_MRI && ~demo
    if ~exist('MRI_part','var')
        MRI_scan = str2double(inputdlg('Which MRI CAT scan are you running? (1 or 2)','Get run',1,{'1'}));
    end
    if MRI_scan == 1 % first scan
        run_start = 19;
        run_end = 19;
    elseif MRI_scan == 2 % first scan
        run_start = 20;
        run_end = 20;
    end
elseif demo
    run_start = 1;
    run_end = 3;
else
    run_start = 1;
    run_end = 20;
end

%==============================================
%% GLOBAL VARIABLES
%==============================================
task_name = 'CAT';
%paths
mainPath=pwd;
outputPath = [mainPath '/Output'];
stim_dir = [mainPath '/Stim/'];
instructions_dir = [mainPath '/Instructions/'];
onsetsPath = [mainPath '/Onset_files/'];
if demo
    stim_dir = [stim_dir 'demo/'];
    instructions_dir = [instructions_dir 'demo/'];
end

% essential for randomization
rng('shuffle');

% about timing
c = clock;
timestamp = sprintf('%04i%02i%02i_%02dh%02dm',c(1),c(2),c(3),c(4),c(5));

% timing
trial_duration = 1; % time to present each image
late_resp_time = 0.45; % time to get late responses
baseline_fixation_duration = 2; % time before the first image
postexp_fixation_duration = 6 - late_resp_time; % time after last image
Initial_Ladder = 0.850; % start with 850 ms Ladder
Ladder_update_insuccessful = -0.0; % reduce ladder if unseccessful response (in seconds. e.g 0.05)
Ladder_update_successful = abs(Ladder_update_insuccessful)/3; % increase ladder if unseccessful response

% give break every X runs
runs_untill_break = 4;
break_counter = runs_untill_break; % ADD BREAK WITH THIS COUNTER
if is_MRI
    total_num_runs = 20; % used to define end of script message
else
    total_num_runs = 18; % used to define end of script message
end
% Assign order
subjectID_num = str2double(subjectID(end-2:end));
order = rem(subjectID_num,4);
if order ==0
    order = 4;
end

%% Load dataframe with info from previous segments
if demo
    stimuli=dir([stim_dir '*.jpg' ]);
    all_stimname={stimuli.name};
    all_stimname = all_stimname(1:6); % use only 6 images for demo
else
    if run_start == 1 % if this is the first run - read the design.
        output_filename = sprintf('%s/%s_%s_%s',outputPath,subjectID,task_name,timestamp);
        info_data_file = dir(sprintf('%s/%s_design_%s_*.txt',outputPath,subjectID,task_name));
        info_data = readtable([info_data_file.folder,'/',info_data_file.name],'Delimiter','\t');
        % add additional info
        info_data.TrialOnset(:)=nan;
        info_data.TrialOffset(:)=nan;
        info_data.Ladder(:)=Initial_Ladder;
        info_data.CueOnset(:)=nan;
        info_data.RT(:)=nan;
        %info_data.KeyPressed(:) = nan;
        %info_data.KeyPressed = char(info_data.KeyPressed);
        info_data.KeyPressed(:) = {'x'};
        info_data.onsetPlanned(:)=nan;
        info_data.onset(:)=nan;
        info_data.duration(:)=nan;
    else % if this is not the first run - read the latest CAT file. This file will be updated
        info_data_file = dir(sprintf('%s/%s_%s_*.txt',outputPath,subjectID,task_name));
        if length(info_data_file)>1
            selected_file = listdlg('PromptString','Found multiple files. Select a file:',...
                'SelectionMode','single','ListSize',[350,250],...
                'ListString',{info_data_file.name});
        else
            selected_file = 1;
        end
        output_filename = [outputPath,'/',info_data_file(selected_file).name];
        info_data = readtable(output_filename,'Delimiter','\t');
        output_filename = output_filename(1:end-4); %remove '.txt' ending
    end
    all_stimname = unique(info_data.StimName);
end

n_trials = length(all_stimname); % number of stimuli per run
%n_total_trials = n_trials*length(run_start:run_end);

%% Load image arrays
all_Images=cell(n_trials,1);
for i=1:n_trials
    all_Images{i}=imread([stim_dir,all_stimname{i}]);
end

% Load Hebrew instructions image files
if is_MRI
    Instructions_Image=imread(sprintf('%s%s_MRI.jpg',instructions_dir,task_name));
else
    Instructions_Image=imread(sprintf('%s%s.jpg',instructions_dir,task_name));
end
%==============================================
%% 'INITIALIZE Screen variables'
%==============================================
Screen('Preference', 'VisualDebuglevel', 3); %No PTB intro screen
screennum = max(Screen('Screens'));
pixelSize=32;
try
    % [w] = Screen('OpenWindow',screennum,[],[0 0 640 480],pixelSize);% %debugging screensize
    [w] = Screen('OpenWindow',screennum,[],[],pixelSize);
catch
    % Skip sync test and suppress visual warning
    Screen('Preference','VisualDebugLevel', 0);
    Screen('Preference', 'SkipSyncTests', 1);
    [w] = Screen('OpenWindow',screennum,[],[],pixelSize);
    % Remove skip sync test
    Screen('Preference', 'SkipSyncTests', 0);
end
clc; %clear psychotoolbox text from command window

% Define Colors
% - - - - - - - - - - - - - - -
black=BlackIndex(w); % Should equal 0.
white=WhiteIndex(w); % Should equal 255.
green=[0 255 0];
Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);

% text stuffs
% - - - - - - - - - - - - - - -
Screen('TextFont',w,'Arial');
Screen('TextSize',w, 60);

% Frame and stack properties
% - - - - - - - - - - - - - - -
[wWidth, wHeight]=Screen('WindowSize', w);
%HideCursor;

% Refresh rate of the screen
% - - - - - - - - - - - - - - -
RefRate = Screen('GetFlipInterval', w) *(0.5);

% Maximum priority level
topPriorityLevel = MaxPriority(w);
Priority(topPriorityLevel);

%% Cue Options
CueDuration=0.100;
% Set up alpha-blending
Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% Load cue image
[CueImage, ~, CueImageAlpha] = imread('./Misc/gabor_cue.png');
% Place the transparency layer of the foreground image into the 4th (alpha)
% image plane. This is the bit which is needed, else the tranparent
% background will be non-transparent
transparencyFactor=0.7;  %change this factor to manipulate cue's transparency
CueImage(:, :, 4) = CueImageAlpha*transparencyFactor ;
imageTextureCue = Screen('MakeTexture', w, CueImage);

% Get the size of the Images
images=dir('./stim/demo/*.jpg');
theImage=imread(['./stim/demo/',images(1).name]);
[image_size,~,~]=size(theImage);
clear images theImage

% Scale the size of the cue as a function of the size of the back image
scaleFactor = 0.4    ; %change this factor to manipulate cue's size
dstRect = CenterRectOnPointd([0 0 image_size image_size].* scaleFactor ,...
    wWidth / 2, wHeight / 2);

% create image textures
Images_texture_all = cell(n_trials,1);
for i = 1:n_trials
    Images_texture_all{i}=Screen('MakeTexture', w, all_Images{i});
end

%%---------------------------------------------------------------
%%  'FEEDBACK VARIABLES'
%%---------------------------------------------------------------
KbName('UnifyKeyNames');
resp_buttons = {'b','y','g','r'}; % optional response keys

%-----------------------------------------------------------------
%% Initializing eye tracking system %
%-----------------------------------------------------------------
if use_eyetracker
    % STEP 2
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    el=EyelinkInitDefaults(w);
    el.backgroundcolour = black;
    el.backgroundcolour = black;
    el.foregroundcolour = white;
    el.msgfontcolour    = white;
    el.imgtitlecolour   = white;
    el.calibrationtargetcolour = el.foregroundcolour;
    EyelinkUpdateDefaults(el);
    
    % STEP 3
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit(0, 1)
        fprintf('Eyelink Init aborted.\n');
    end
    
    % make sure we're still connected.
    if Eyelink('IsConnected')~=1
        fprintf('not connected, clean up\n');
        Eyelink( 'Shutdown');
        Screen('CloseAll');
        return;
    end
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,HREF,AREA');
    % open file to record data to
    edfFile = ['CAT_', subjectID(end-2:end), '.edf']; % 8 characters limit
    i = Eyelink('Openfile', edfFile);
    open_edf_attempt = 0;
    while i~=0 && open_edf_attempt<=5
        i = Eyelink('Openfile', edfFile);
        open_edf_attempt = open_edf_attempt + 1;
        fprintf('\nFailed to open EDF. Trying again. Attempt: %i\n',open_edf_attempt);
    end
    if i~=0
        fprintf('\nFailed to open EDF too many times. Aborting script.');
        return
    end
    
    % SET UP TRACKER CONFIGURATION
    % Setting the proper recording resolution, proper calibration type,
    % as well as the data file content;
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, wWidth-1, wHeight-1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, wWidth-1, wHeight-1);
    % set calibration type.
    Eyelink('command', 'calibration_type = HV9');
    % set parser (conservative saccade thresholds)
    
    % STEP 4
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    % do a final check of calibration using driftcorrection
    EyelinkDoDriftCorrection(el);
end

%% onsets
%---------------------------
aftertrialfixation = .5;
TotalTrialDuration = trial_duration + aftertrialfixation;
if is_MRI && ~demo
    % load 1 out of 4 random onsets lists
    %load(sprintf('%s/CAT_onsets_%i',onsetsPath, randperm(4,1)),'onsets');
    load(sprintf('%s/%s_onsets_%i',onsetsPath,task_name, order),'onsets');
else
    % To be used in demo and behavioral
    onsets=0:TotalTrialDuration:TotalTrialDuration*(n_trials-1);
end
run_onsets = onsets;
planned_onsets = run_onsets + baseline_fixation_duration;

%---------------------------------------------------------------
%% 'Display Main Instructions'
%---------------------------------------------------------------
Screen('PutImage',w,Instructions_Image);
Screen(w,'Flip');
WaitSecs(0.01);

KbQueueCreate;
noresp = 1;
while noresp
    [keyIsDown] = KbCheck(-1); % deviceNumber=keyboard
    if keyIsDown && noresp
        noresp = 0;
    end
end

% Wait for 't' in the MRI
if is_MRI == 1 && (~demo)
    CenterText(w,'GET READY! Waiting for trigger', white, 0, 0);
    Screen(w,'Flip');
    escapeKey = KbName('t');
    while 1
        [keyIsDown,~,keyCode] = KbCheck(-1);
        max(keyCode);
        if keyIsDown && keyCode(escapeKey)
            break;
        end
    end
    DisableKeysForKbCheck(KbName('t')); % So trigger is no longer detected
end % end if test_comp == 1

WaitSecs(0.01);
CenterText(w,'+', white,0,0);
TaskStartTime = Screen(w,'Flip');

KbQueueCreate;
if use_eyetracker
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.05);
    Eyelink('StartRecording');
    WaitSecs(.01);
end
%for run_number = run_start:min(run_start+2,20) % for debugging 2 runs starting with runInd
for run_number = run_start:run_end %this for loop allows all runs in block to be completed
    
    CenterText(w,'+');
    Screen('Flip',w);
    if use_eyetracker
        Eyelink('Message', Eventflag(GenFlags.RunStart.str,task_name,run_number,0,TaskStartTime)); % mark start time in file
    end
    
    %% get info for for the current run
    %---------------------------
    if ~demo
        run_info_ind = info_data.Run == run_number;
        run_stim_names = info_data.StimName(run_info_ind);
        isGoTrial = info_data.Go(run_info_ind);
        run_ladders = info_data.Ladder(run_info_ind);
    elseif demo
        run_stim_names = all_stimname;
        % in the demo - have last 2 trials go in the first run and 1 in the
        % second run
        n_go_stim = 4-run_number;
        isGoTrial = [zeros(n_trials-n_go_stim,1);ones(n_go_stim,1)];
        run_ladders = repmat(Initial_Ladder,[n_trials,1]);
    end
    % images index
    [~,img_location] = ismember(run_stim_names,all_stimname);
    
    % preallocation for recorded variable
    [actual_onset_time,actual_offset_time,actual_cue_time, RT] = deal(nan(n_trials,1));
    key_pressed = repmat('x',[n_trials,1]);
    
    %% Trials
    % for trialNum = 1:6 % shorter version for debugging (max 6 for demo)
    for trialNum = 1:n_trials   % To cover all the items in one run.
        CenterText(w,'+');
        Screen('Flip',w);
        
        %% Display image
        Screen('DrawTextures',w,Images_texture_all{img_location(trialNum)});
        image_start_time = Screen('Flip',w, TaskStartTime + planned_onsets(trialNum) - RefRate); % display images according to Onset times
        actual_onset_time(trialNum) = image_start_time - TaskStartTime ;
        
        if use_eyetracker
            %   Eyelink MSG: event, task, run, trial, runStart
            Eyelink('Message', Eventflag(GenFlags.TrialStart.str,task_name,run_number,trialNum,TaskStartTime)); % mark start time in file
        end
        
        %% 'EVALUATE RESPONSE & ADJUST LADDER ACCORDINGLY'
        %---------------------------------------------------
        noresp = 1;
        KbQueueFlush;
        KbQueueStart;
        
        t1_trial_onset = image_start_time; % part 1 - show image
        t2_cue_onset = image_start_time + run_ladders(trialNum) - RefRate; %  part 2 - show cue
        t3_cue_offset = image_start_time + run_ladders(trialNum) + CueDuration ; %  part 3 - show image
        t4_fixation_onset = image_start_time + trial_duration ; %  part 4 - show fixation
        completed_stage = 0;
        
        while GetSecs < image_start_time + trial_duration + late_resp_time
            %% look for response while listening
            if noresp
                [pressed, firstPress, ~, ~, ~] = KbQueueCheck;
                if pressed
                    firstKeyPressed = find(firstPress>0);
                    if length(firstKeyPressed)>=2
                        firstKeyPressed=firstKeyPressed(1);
                    end
                    tmp = KbName(firstKeyPressed);
                    if ~ischar(tmp) % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed{runnum} to be a char, so this converts it and takes the first key pressed
                        tmp = char(tmp);
                    end
                    key_pressed(trialNum) = tmp(1);
                    
                    if ismember(key_pressed(trialNum),resp_buttons)
                        noresp = 0;
                        RT(trialNum) = firstPress(firstKeyPressed)-image_start_time;
                        if use_eyetracker
                            Eyelink('Message', Eventflag(GenFlags.Response.str,task_name,run_number,0,TaskStartTime)); % mark start time in file
                        end
                    end
                end
            end
            
            if GetSecs > t4_fixation_onset && completed_stage < 4
                %% Stage 4: Show fixation
                CenterText(w,'+');
                TrialOffset = Screen('Flip',w);
                actual_offset_time(trialNum) = TrialOffset - TaskStartTime;
                if use_eyetracker
                    Eyelink('Message', Eventflag(GenFlags.Fixation.str,task_name,run_number,trialNum,TaskStartTime)); % mark start time in file
                end
            elseif GetSecs > t3_cue_offset && completed_stage < 3
                %% Stage 3: Display image and prepare next stage: Stim offset
                Screen('DrawTextures',w,Images_texture_all{img_location(trialNum)});
                Screen('Flip',w);
                CenterText(w,'+', white,0,0);
            elseif GetSecs > t2_cue_onset && completed_stage < 2  && isGoTrial(trialNum) == 1
                %% Stage 2 (Go trials only): Display image and prepare next stage: Cue offset
                cue_time = Screen('Flip',w);
                actual_cue_time(trialNum) = cue_time - image_start_time;
                if use_eyetracker
                    Eyelink('Message', Eventflag(GenFlags.CueStart.str,task_name,run_number,trialNum,TaskStartTime)); % mark start time in file
                end
                completed_stage = 2;
                % next step preparation
                Screen('DrawTextures',w,Images_texture_all{img_location(trialNum)});
            elseif GetSecs > t1_trial_onset && completed_stage < 1  && isGoTrial(trialNum) == 1
                %% Stage 1 (Go trials only): Display image and prepare next stage: Cue onset
                completed_stage = 1;
                % next step preparation
                Screen('DrawTextures',w,Images_texture_all{img_location(trialNum)});
                Screen('DrawTextures', w, imageTextureCue, [], dstRect, 0);
            end
        end %%% End big while waiting for response within 1450 msec
        
        % Update ladder after successful/insuccessful Go trials
        if isGoTrial(trialNum) == 1 && (RT(trialNum) >= trial_duration || isnan(RT(trialNum)))
            ladder_update = Ladder_update_insuccessful;
        elseif isGoTrial(trialNum) == 1 && (RT(trialNum) < trial_duration)
            ladder_update = Ladder_update_successful;
        else
            ladder_update = 0;
        end
        try % if there are any more trials left in the run
            run_ladders(trialNum+1:end) = run_ladders(trialNum) + ladder_update;
        catch
        end
        if ~ demo
            % update data
            trial_location = (info_data.Run == run_number) & (info_data.TrialInRun == trialNum);
            info_data.onsetPlanned(trial_location) = planned_onsets(trialNum);
            info_data.TrialOnset(trial_location) = actual_onset_time(trialNum);
            info_data.TrialOffset(trial_location) = actual_offset_time(trialNum);
            try % if there are any more trials left in the task
                info_data.Ladder(find(trial_location)+1:end) = info_data.Ladder(trial_location) + ladder_update;
            catch
            end
            info_data.CueOnset(trial_location) = actual_cue_time(trialNum)*1000;
            info_data.RT(trial_location) = RT(trialNum)*1000;
            info_data.KeyPressed{trial_location} = key_pressed(trialNum);
        end
    end %	End the big trialNum loop showing all the images in one run.
    
    %% Break every few run
    break_counter = break_counter - 1;
    if break_counter == 0 && run_number~=run_end
        CenterText(w,'This is a short break',white, 0,-270);
        CenterText(w,'press when you are ready to continue', green, 0, -180);
        Screen('Flip',w);
        
        KbQueueFlush;
        KbQueueCreate;
        noresp = 1;
        while noresp
            [keyIsDown] = KbCheck(-1); % deviceNumber=keyboard
            if keyIsDown && noresp
                noresp = 0;
            end
        end
        KbQueueFlush;
        % reset the break counter
        break_counter = runs_untill_break;
        CenterText(w,'+', white,0,0);
        Screen(w,'Flip');
    end
    RunEndTime = planned_onsets(trialNum) + TotalTrialDuration;
    planned_onsets = run_onsets + RunEndTime ;
end % End the run loop to go over all the runs

postexperiment = GetSecs;
CenterText(w,'+', white,0,0);
Screen(w,'Flip');
while GetSecs < postexperiment+postexp_fixation_duration
end

%---------------------------------------------------------------
%%   save data to a .mat file & close out
%---------------------------------------------------------------
vars_in_workspace = who;
vars2clear = vars_in_workspace(contains(vars_in_workspace,'Image'));
for i = 1:length(vars2clear)
    clear(vars2clear{i});
end

if ~demo
    % compatible with eye tracker: add onset and duration
    info_data.onset = info_data.TrialOnset;
    info_data.duration = info_data.TrialOffset - info_data.TrialOnset;
    % save data and workspace
    writetable(info_data,[output_filename,'.txt'],'Delimiter','\t');
    if is_MRI
        % save data only of MRI scan
        output_filename = sprintf('%s/%s_CATMRI_scan_%i_%s',outputPath,subjectID,MRI_scan,timestamp);
        info_data_MRI_scan = info_data(info_data.Run == run_number,:);
        writetable(info_data_MRI_scan,[output_filename,'.txt'],'Delimiter','\t');
    end
    save(output_filename);
end

if use_eyetracker
    % finish up: stop recording eye-movements,
    % close graphics window, close data file
    Eyelink('StopRecording');
    WaitSecs(.1);
    Eyelink('CloseFile');
    
    % download data file
    if ~demo
        try
            fprintf('Receiving data file ''%s''\n', edfFile );
            status=Eyelink('ReceiveFile');
            if status > 0
                fprintf('ReceiveFile status %d\n', status);
            end
            if 2==exist(edfFile, 'file')
                fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
            end
        catch
            fprintf('Problem receiving data file ''%s''\n', edfFile );
        end
        edf_filename = [output_filename,'.edf'];
        movefile(edfFile,edf_filename);
        [~,tmp] = system([mainPath '/edf2asc ',edf_filename]);
        converted_ok = contains(tmp,'successfully');
        if ~converted_ok
            disp('Could not download EDF file!\n');
        end
    end
    Eyelink('ShutDown');
end
%   outgoing msg & closing
% ------------------------------
if demo
    CenterText(w,sprintf('End of demo. Any questions?') ,white,0,-170);
elseif is_MRI && (run_number < total_num_runs) % if this is not the last run
    CenterText(w,sprintf('Another run will begin soon'), white, 0,-300);
else % if this is the last run
    CenterText(w,'Great Job. Thank you!',green, 0,-270);
    CenterText(w,'Now we will continue to the next part', white, 0, -180);
end
DisableKeysForKbCheck([]); % Revert disable key
Screen('Flip',w);
WaitSecs(4);
KbQueueFlush;
Screen('CloseAll');
ShowCursor;

