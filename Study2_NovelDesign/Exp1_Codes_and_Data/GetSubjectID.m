function [subjectID] = GetSubjectID(exp_code)
% By: Tom Salomon, 2019
% This function uses the code 'whats_the_next_sub.m' to get a suggested
% code for the next participat. User should define the prefix (exp_code) to
% be used in the SubjectID either as an input or through input dialog box.

if nargin <1
    exp_code =[];
end
exp_code = char(exp_code);
okCode = ~isempty(exp_code);

while ~okCode
    exp_code = char(inputdlg('Experiment code not specified. Please enter desired prefix, e.g. "MRI_faces"','Experiment code',1));
    if ~isempty(exp_code) % input checkers
        okCode = true;
    else
        okCode = false;
    end
end

suggested_subID = whats_the_next_sub('read'); % get suggestion for the next participant code
okID = false;
while ~okID % make sure code is 3 digit numeric
    subjectID_num = str2double(inputdlg('Please insert subject code (3 digits only). 999 for debugging',...
        'subject code',1,{num2str(suggested_subID)}));
    if isnumeric(subjectID_num) && (~isempty(subjectID_num))% input checkers
        okID = subjectID_num<=999 && subjectID_num>=101;
    else
        okID = false;
    end
end
subjectID = sprintf('%s_%03i',exp_code, subjectID_num);
if subjectID_num<999 % Update next ID only if not 999 (debugging)
    next_sub_value = 1 + str2double(subjectID(end-2:end));
    whats_the_next_sub('update',next_sub_value);
end

end % end of function
