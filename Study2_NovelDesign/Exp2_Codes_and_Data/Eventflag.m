%% This function creates formated flags for events during eyetracking experiments.
% takes input of Run and trial as number values
%requires ENUM- GenFlags 

function event= Eventflag(Event,Task,Run,Trial,runStart)

Run = sprintf('%03i',Run);
Trial = sprintf('%03i',Trial);

event=['flag_',Event,'_Task',Task,'_Run',Run,'_Trial',Trial,'_Time',num2str(GetSecs-runStart)];

end