
function binary_ranking(subjectID,demo,use_eyetracker)

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ========================== by Tom Salomon 2019 ==========================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% Ranking using a set of binary choices.
% This function was written based on the probe code. It is used to rank 60
% different items, based on subject's selection. In each trial two images
% will appear on screen, and the subject is asked to choose one, using the
% 'u' and 'i' keys. Each stimulus will be presented exactly 10 times, each
% time in a different comparison, for a total number of 300 unique
% comparisons.
% Finally, Colley ranking code is run on the subject's choices to create a
% ranking list of all items.

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FUNCTIONS REQUIRED TO RUN PROPERLY: ----------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   'random_stimlist_generator'
% % %   'colley'
% % %   'edf2asc' - if using eyetracker

%=========================================================================
%% dummy info for testing purposes
%=========================================================================

% subjectID='BM2_000';
% test_comp=0;
% path=pwd

%==============================================
%% GLOBAL VARIABLES
%==============================================
% default input
% - - - - - - - - - - - - - - - - -
mainPath='.';
outPath = [mainPath,'/Output'];
task_name = 'binary_ranking';
run_number = 1;
MRI_comp=0;
if nargin<3
    ask_eyetracker = questdlg('Do you want to use the eye tracker?','Use eye tracker','Yes','No','Yes');
    use_eyetracker = strcmp(ask_eyetracker,'Yes');
    if nargin<2
    demo_ans = questdlg('Do you want to run to full task or demo?','Ask Demo',...
        'Full Task','Demo','Full Task');
        demo = strcmp(demo_ans,'Demo'); % set to 1/0 to run demo/full task
    end
end
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

% timestamp
c = clock;
timestamp=sprintf('%04i%02i%02i_%02dh%02dm', c(1),c(2),c(3),c(4),c(5));

% define duration variables
maxtime=2;      % 2 second limit on each selection
baseline_fixation_dur=4; % Need to modify based on if first few volumes are saved or not
aftertrialoutcome_dur = 0.5; % time to show 'responde faster' / confirmation green rect
aftertrialfixation_dur = 0.5;
aftertaskfixation = 2;
TotalTrialDuration = maxtime +  aftertrialfixation_dur + aftertrialoutcome_dur;
tic

% essential for randomization
rng('shuffle');

%==============================================
%% Read in data
%==============================================

% load image arrays
% - - - - - - - - - - - - - - -
% Read all image files in the 'stim' folder. make sure the ending is suitable
% to your stimuli, and the folder contains only the experiment's 60 images.
stim_dir = [mainPath '/Stim/'];
instructions_dir = [mainPath '/Instructions/'];

if demo
    stim_dir = [stim_dir 'demo/'];
    instructions_dir = [instructions_dir 'demo/'];
end
stimuli=dir([stim_dir '*.jpg' ]);

stimname={stimuli.name};
Images=cell(length(stimname),1);
for i=1:length(stimname)
    Images{i}=imread([stim_dir,stimname{i}]);
end

% Load Hebrew instructions image files
Instructions_Image=imread(sprintf('%s%s.jpg',instructions_dir,task_name));

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

% Define image scale - Change according to your stimuli
% - - - - - - - - - - - - - - -
stackH = size(Images{1},1);
stackW = size(Images{1},2);

% Frame and stack properties
% - - - - - - - - - - - - - - -
[wWidth, wHeight]=Screen('WindowSize', w);
xcenter=wWidth/2;
ycenter=wHeight/2;
xDist = 300; % distance from center in the x axis. can be changed to adjust to smaller screens
leftRect=[xcenter-stackW-xDist ycenter-stackH/2 xcenter-xDist ycenter+stackH/2]; % left stack location
rightRect=[xcenter+xDist ycenter-stackH/2 xcenter+stackW+xDist ycenter+stackH/2]; % right stack location
penWidth=10; % frame width
HideCursor;

% Refresh rate of the screen
% - - - - - - - - - - - - - - -
% used to flip one frame before you want the stim to appeare
RefRate = Screen('GetFlipInterval', w) *(0.5);

% Maximum priority level
topPriorityLevel = MaxPriority(w);
Priority(topPriorityLevel);

%==============================================
%% 'ASSIGN response keys'
%==============================================
KbName('UnifyKeyNames');
switch MRI_comp
    case 0 % Experiment room
        leftstack='u';
        rightstack= 'i';
        badresp='x';
    case 1 % MRI response box
        leftstack='b';
        rightstack= 'y';
        badresp='x';
end
%==============================================
%%   'PRE-TRIAL DATA ORGANIZATION'
%==============================================
% Stimuli lists
% - - - - - - - - - - - - - - -
number_of_stimuli=length(stimname); % number of stimuli
if demo
    number_of_trials=8; % desired number of trials (number of unique comparisons)
else
    number_of_trials=400; % desired number of trials (number of unique comparisons)
    %number_of_trials=10; % for debugging
end

% Define onsets
% - - - - - - - - - - - - - - -
onsets=0:TotalTrialDuration:TotalTrialDuration*(number_of_trials-1);
planned_onsets = onsets + baseline_fixation_dur;

% IMPORTANT NOTE: every stimulus will be presented exactly 2*number_of_trials/n times. Therefore, number_of_trials should be a multiple of n/2. e.g if n=60, number_of_trials can be 30,60,90,120....
% ==============                                                                                  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[shuffle_stimlist1,shuffle_stimlist2]=random_stimlist_generator(number_of_stimuli,number_of_trials);

% Colley's ranking analysis variables
% - - - - - - - - - - - - - - -
N=zeros(number_of_stimuli,3); %create the N matrix of (win,loses,total) for each stimulus (each row represent a specific stimulus). To be used in Colley's ranking
T=zeros(number_of_stimuli,number_of_stimuli); %create the T symetric matrix of "competitions"  (each row represent a specific stimulus). T(a,b)=T(b,a)=-1 means a "competition" had taken place between stimuli a&b. To be used in Colley's ranking
for stimulus=1:number_of_stimuli
    T(stimulus,shuffle_stimlist1(shuffle_stimlist2==stimulus))=-1;
    T(shuffle_stimlist1(shuffle_stimlist2==stimulus),stimulus)=-1;
    T(stimulus,shuffle_stimlist2(shuffle_stimlist1==stimulus))=-1;
    T(shuffle_stimlist2(shuffle_stimlist1==stimulus),stimulus)=-1;
    N(stimulus,3)=(sum(shuffle_stimlist1==stimulus)+sum(shuffle_stimlist2==stimulus));
end

%-----------------------------------------------------------------
%% 'Write output file header'
%-----------------------------------------------------------------
if ~demo
    output_filename=sprintf('%s/%s_binary_ranking_%s',outPath,subjectID,timestamp);
    fid1=fopen([output_filename,'.txt'], 'a');
    fprintf(fid1,'subjectID\truntrial\tonsettime\tImageLeft\tImageRight\tStimNumLeft\tStimNumRight\tResponse\tOutcome\tRT\n'); %write the header line
end
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
    % Disable key output to Matlab window:
    
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
        return;
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
    edfFile = ['ip_', subjectID(end-2:end), '.edf']; % 8 characters limit
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
%==============================================
%% 'Display Main Instructions'
%==============================================
KbQueueCreate;
Screen('PutImage',w,Instructions_Image);
Screen(w,'Flip');

noresp=1;
while noresp
    [keyIsDown] = KbCheck(-1); % deviceNumber=keyboard
    if keyIsDown && noresp
        noresp=0;
    end
end

%==============================================
%% 'Run Trials'
%==============================================
CenterText(w,'+', white,0,0);
TaskStartTime = Screen(w,'Flip');

if use_eyetracker
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.05);
    Eyelink('StartRecording');
    WaitSecs(0.05);
    % Event,Task,Run,Trial,runStart
    Eyelink('Message', Eventflag(GenFlags.RunStart.str,task_name,run_number,0,TaskStartTime)); % mark start time in file
end

for trial=1:number_of_trials
    
    chose_rand=rand; % randomly chose left/right location od stimuli
    if chose_rand<=0.5
        leftname=stimname(shuffle_stimlist1(trial));
        rightname=stimname(shuffle_stimlist2(trial));
    else
        leftname=stimname(shuffle_stimlist2(trial));
        rightname=stimname(shuffle_stimlist1(trial));
    end
    
    out=999000;
    
    %-----------------------------------------------------------------
    % display images
    %-----------------------------------------------------------------
    if chose_rand<=0.5
        Screen('PutImage',w,Images{shuffle_stimlist1(trial)}, leftRect);
        Screen('PutImage',w,Images{shuffle_stimlist2(trial)}, rightRect);
    else
        Screen('PutImage',w,Images{shuffle_stimlist2(trial)}, leftRect);
        Screen('PutImage',w,Images{shuffle_stimlist1(trial)}, rightRect);
    end
    
    
    
    CenterText(w,'+', white,0,0);
    StimOnset=Screen(w,'Flip', TaskStartTime + planned_onsets(trial) - RefRate);
    if use_eyetracker
        %   Eyelink MSG: event, task, run, trial, runStart
        Eyelink('Message', Eventflag(GenFlags.TrialStart.str,task_name,run_number,trial,TaskStartTime)); % mark start time in file
    end
    
    %-----------------------------------------------------------------
    % get response
    %-----------------------------------------------------------------
    KbQueueFlush;
    KbQueueStart;
    noresp=1;
    goodresp=0;
    while noresp
        
        % check for response
        [keyIsDown, firstPress] = KbQueueCheck;
        
        if keyIsDown && noresp
            keyPressed=KbName(firstPress);
            if ischar(keyPressed)==0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                keyPressed=char(keyPressed);
                keyPressed=keyPressed(1);
            end
            % valid response
            if strcmp(keyPressed,leftstack) || strcmp(keyPressed,rightstack)
                if use_eyetracker
                    %   Eyelink MSG: event, task, run, trial, runStart
                    Eyelink('Message', Eventflag(GenFlags.Response.str,task_name,run_number,trial,TaskStartTime)); % mark start time in file
                end
                respTime=firstPress(firstPress>0)-StimOnset;
                noresp=0;
                goodresp=1;
            end
        end
        
        % check for reaching time limit
        if noresp && GetSecs-TaskStartTime >= planned_onsets(trial) + maxtime
            noresp=0;
            keyPressed=badresp;
            respTime=maxtime;
        end
    end
    
    %-----------------------------------------------------------------
    % determine which choice to highlight
    %-----------------------------------------------------------------
    
    if goodresp==1 % if responded: add green rectangle around selected image
        if chose_rand<=0.5
            Screen('PutImage',w,Images{shuffle_stimlist1(trial)}, leftRect);
            Screen('PutImage',w,Images{shuffle_stimlist2(trial)}, rightRect);
        else
            Screen('PutImage',w,Images{shuffle_stimlist2(trial)}, leftRect);
            Screen('PutImage',w,Images{shuffle_stimlist1(trial)}, rightRect);
        end
        
        if keyPressed==leftstack
            Screen('FrameRect', w, green, leftRect, penWidth);
        elseif keyPressed==rightstack
            Screen('FrameRect', w, green, rightRect, penWidth);
        end
        
        CenterText(w,'+', white,0,0);
        Screen(w,'Flip');
        
    else % if did not respond: show text 'You must respond faster!'
        CenterText(w,sprintf('You must respond faster!') ,white,0,0);
        Screen(w,'Flip');
        if use_eyetracker
            %   Eyelink MSG: event, task, run, trial, runStart
            Eyelink('Message', Eventflag(GenFlags.RespondFaster.str,task_name,run_number,trial,TaskStartTime)); % mark start time in file
        end
    end
    
    %-----------------------------------------------------------------
    % show fixation ITI
    %-----------------------------------------------------------------
    CenterText(w,'+', white,0,0);
    Screen(w,'Flip',TaskStartTime + planned_onsets(trial) + respTime + aftertrialoutcome_dur - RefRate);
    if use_eyetracker
        %   Eyelink MSG: event, task, run, trial, runStart
        Eyelink('Message', Eventflag(GenFlags.Fixation.str,task_name,run_number,trial,TaskStartTime)); % mark start time in file
    end
    
    if goodresp ~= 1
        respTime=999000;
        % Colley ranking input - remove ties (no choices) from calculation.
        T(shuffle_stimlist1(trial),shuffle_stimlist2(trial))=0;
        T(shuffle_stimlist2(trial),shuffle_stimlist1(trial))=0;
        N(shuffle_stimlist1(trial),3)=N(shuffle_stimlist1(trial),3)-1;
        N(shuffle_stimlist2(trial),3)=N(shuffle_stimlist2(trial),3)-1;
        
    end
    
    %-----------------------------------------------------------------
    % write to output file
    %-----------------------------------------------------------------
    if ~demo
        if chose_rand<=0.5
            fprintf(fid1,'%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\n', subjectID, trial, StimOnset-TaskStartTime, char(leftname), char(rightname), shuffle_stimlist1(trial), shuffle_stimlist2(trial), keyPressed, out, respTime*1000);
        else
            fprintf(fid1,'%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\n', subjectID, trial, StimOnset-TaskStartTime, char(leftname), char(rightname), shuffle_stimlist2(trial), shuffle_stimlist1(trial), keyPressed,  out, respTime*1000);
        end
    end
    
    % add trial info to the Colley ranking mats
    if chose_rand<=0.5
        if keyPressed==leftstack
            N(shuffle_stimlist1(trial),1)=N(shuffle_stimlist1(trial),1)+1;
            N(shuffle_stimlist2(trial),2)=N(shuffle_stimlist2(trial),2)+1;
        elseif keyPressed==rightstack
            N(shuffle_stimlist2(trial),1)=N(shuffle_stimlist2(trial),1)+1;
            N(shuffle_stimlist1(trial),2)=N(shuffle_stimlist1(trial),2)+1;
        end
    else
        if keyPressed==rightstack
            N(shuffle_stimlist1(trial),1)=N(shuffle_stimlist1(trial),1)+1;
            N(shuffle_stimlist2(trial),2)=N(shuffle_stimlist2(trial),2)+1;
        elseif keyPressed==leftstack
            N(shuffle_stimlist2(trial),1)=N(shuffle_stimlist2(trial),1)+1;
            N(shuffle_stimlist1(trial),2)=N(shuffle_stimlist1(trial),2)+1;
        end
    end
    
    KbQueueFlush;
end % loop through trials

Postexperiment=GetSecs;
if use_eyetracker
    %   Eyelink MSG: event, task, run, trial, runStart
    Eyelink('Message', Eventflag(GenFlags.RunEnd.str,task_name,run_number,trial,TaskStartTime)); % mark start time in file
end

while GetSecs < Postexperiment+aftertaskfixation
    CenterText(w,'+', white,0,0);
    Screen(w,'Flip');
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

% display end of part message
WaitSecs(1);
Screen('FillRect', w, black);
if demo
    CenterText(w,sprintf('End of demo. Any questions?') ,white,0,-170);
else
    CenterText(w,sprintf('Thank you! Please call the experimenter.') ,white,0,-170);
end
Screen('Flip',w);
WaitSecs(5);
Screen('CloseAll');

if ~demo
    % Run Colley's ranking
    stimuli_ranking=colley(T,N);
    fid2=fopen([outPath,'/' subjectID '_ItemRankingResults_' timestamp '.txt'], 'a');
    fprintf(fid2,'Subject\tStimName\tStimNum\tRank\tWins\tLoses\tTotal\n');
    for j=1:number_of_stimuli
        fprintf(fid2,'%s\t%s\t%d\t%d\t%d\t%d\t%d\n', subjectID, char(stimname(j)), j, stimuli_ranking(j), N(j,1), N(j,2), N(j,3));
    end
    fclose(fid2);
    fclose(fid1);
    
    %   save data to a .mat file
    vars_in_workspace = who;
    vars2clear = vars_in_workspace(contains(vars_in_workspace,'Image'));
    for i = 1:length(vars2clear)
        clear(vars2clear{i});
    end
    save([output_filename,'.mat']); % save wprkspace
end

KbQueueFlush;
end
