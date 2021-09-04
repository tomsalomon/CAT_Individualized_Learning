
function personal_details(subjectID, order, outputPath, sessionNum)
% function personalDetails(subjectID, order, mainPath)
%   This function gets a few personal details from the subject and saves it
%   to file named subjectID '_personalDetails' num2str(sessionNum) '.txt'


% get time and date
c = clock;
hr = sprintf('%02d', c(4));
min = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',min,'m'];

% open a txt file for the details
fid1 = fopen([outputPath '/' subjectID '_personalDetails' num2str(sessionNum) '_' timestamp '.txt'], 'a');
fprintf(fid1,'subjectID\torder\tdate\tgender(1-female, 2-male)\tage\tdominant hand (1-right, 2-left)\toccupation\n'); %write the header line

% ask the subject for the details
Gender = questdlg('Please select your gender:','Gender','Female','Male','Female');
while isempty(Gender)
    Gender = questdlg('Please select your gender:','Gender','Female','Male','Female');
end
if strcmp(Gender,'Male')
    Gender = 2;
else
    Gender = 1;
end

Age = inputdlg ('Please enter your age: ','Age',1);
while isempty(Age) || isempty(Age{1})
    Age = inputdlg ('Only integers between 18 and 40 are valid. Please enter your age: ','Age',1);
end
Age = cell2mat(Age);
Age = str2double(Age);
while mod(Age,1) ~= 0 || Age < 18 || Age > 40
    Age = inputdlg ('Only integers between 18 and 40 are valid. Please enter your age: ','Age',1);
    Age = cell2mat(Age);
    Age = str2double(Age);
end

DominantHand = questdlg('Please select your domoinant hand:','Dominant hand','Left','Right','Right');
while isempty(DominantHand)
    DominantHand = questdlg('Please select your domoinant hand:','Dominant hand','Left','Right','Right');
end
if strcmp(DominantHand,'Left')
    DominantHand = 2;
else
    DominantHand = 1;
end

Occupation = inputdlg('Please type your occupation (for example- a student for Psychology): ','Occupation',1);
Occupation = cell2mat(Occupation);

% Write details to file
fprintf(fid1,'%s\t%d\t%s\t%d\t%d\t%d\t%s\n', subjectID, order, timestamp, Gender, Age, DominantHand, Occupation);
fclose(fid1);
end


% 
% function personal_details (subjectID,order,mainPath)
% 
% if nargin<3
%     mainPath=pwd;
% end
% 
% % Get the date and time
% c = clock;
% hr = sprintf('%02d', c(4));
% min = sprintf('%02d', c(5));
% timestamp=[date,'_',hr,'h',min,'m'];
% time=[hr,'h',min,'m'];
% 
% % Ask subject's age, gender and education
% age = inputdlg ('Please enter your age: ','Age',1);
% gender = questdlg('Please select your gender:','Gender','Male','Female','Male');
% education = questdlg('Please select your eduactional status:','Education'...
%     ,'BA Student','MA/PhD Student','Else','BA Student');
% dom_hand = questdlg('What is your dominant hand?','Dominant Hand','Left','Right','Right');
% 
% fid=fopen([mainPath '/Output/' subjectID '_personal_details_' timestamp '.txt'], 'a');
% 
% fprintf(fid,'subjid\torder\tdate\ttime\tgender\tage\teducation\tdominant_hand\n'); %write the header line
% fprintf(fid,'%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\n', subjectID, order, date, time, gender, str2double(age{1}), education, dom_hand);
% fclose(fid);
% 
% end
