function next_sub = whats_the_next_sub(action,next_sub_value)
% By: Tom Salomon, 2019
% Subjects code tracking function. 
% Action input can be: 1 - read, 2 - update , 3 - create
% 

cloud_path = '~/Dropbox/experimentsOutput/TomSalomon/Boost_BayesianModel/';
next_sub_file = [cloud_path,'next_sub.mat'];

try
    load(next_sub_file,'next_sub');
catch
    next_sub = 1;
    if isfolder(cloud_path)
        action = 'create';
    else
        error('Cloud directory %s does not exist.\nTry redifining it first within the script',cloud_path);
    end
end

if ~exist('action','var')
    action = 'read';
end
if ~exist('next_sub_value','var')
    next_sub_value = 1 + next_sub;
end

if isnumeric(action)
possible_actions = {'read','update','create'};
action = possible_actions{action};
end

switch action
    case 'read'
    case 'update'
        next_sub = next_sub_value;
        save(next_sub_file,'next_sub')
                fprintf('\nUpdated the next_sub file');
    case 'create'
        next_sub = 1;
        save(next_sub_file,'next_sub')
        fprintf('\nCreated a next_sub file');
end
        fprintf('\nnext participant code should be %i\n',next_sub);

