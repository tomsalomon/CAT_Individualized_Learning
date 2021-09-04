function sync_with_cloud(subjectID,export)
% this function copy data from the local experiment path to the cloud
% directory (export = 1, default) or vice versa (export = 0, import)

% define these paths

% export path - where data is saved in the experiment computer
export_path = './Output/';
% cloud path - where data is copied to/from on the cloud
cloud_path =  '~/Dropbox/experimentsOutput/TomSalomon/Boost_BayesianModel/';
% import path - where data is copied to in the analysis computer
import_path = './Output/';

% Experiment code (prefix before subject code)
ExpCode = 'BM';

if ~exist('export','var') % default - export data
    export = 1;
end

if ~exist('subjectID','var')
    subjectID = sprintf('%s_*',ExpCode);
end

if isnumeric(subjectID)
    subjectID = sprintf('%s_%03i',ExpCode,subjectID);
end

if export % export from experimental computer to cloud
    origin_dir_behave = export_path;
    target_dir_behave = cloud_path;
else % import from cloud to analysis computer
    origin_dir_behave = cloud_path;
    target_dir_behave = import_path;
end

sub_data_behave = dir([origin_dir_behave,'/*',subjectID,'*']);

for file_i = 1:numel(sub_data_behave)
    file_name = sub_data_behave(file_i).name;
    copyfile([origin_dir_behave,'/',file_name],[target_dir_behave,'/',file_name]);
end
