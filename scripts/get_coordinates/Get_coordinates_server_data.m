
% this code gets and organizes coordinates from the server's data

clc; clear;

file_path = fullfile(fileparts(pwd),'ProcessedData','SubConds.xlsx');
table_data = readtable(file_path);
Subs = table_data.SubID; % this matches SubIDlist
StimLoc = table_data.StimLoc;

studydir = '/Users/liuq13/NODEAP/MRI';
SubIDlist = dir(fullfile(studydir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

wd = pwd; % get current working directory

%% load TMS coordinates in subjects' native space

StimCoord = zeros(nSubIDlist,3);
for subj = 1:nSubIDlist
    sn = SubIDlist(subj).name;
    subdir = fullfile(studydir,sn);
    sub_mat_file = fullfile(subdir,'stimulation_coordinates.mat');
    if ~exist(sub_mat_file,'file')
        disp(sn)
        continue;
    end
    load(sub_mat_file)
    if strcmp(StimLoc(subj),'Posterior')
        StimCoord(subj,:) = regions(2).StimCoord;
    end
    if strcmp(StimLoc(subj),'Anterior')
        StimCoord(subj,:) = regions(1).StimCoord;
    end
end

%%
varNames = {'x','y','z'};
Data = array2table(StimCoord, 'VariableNames', varNames);

% add subject ID
Data.SubID = Subs;
Data.StimLoc = StimLoc;

Data = Data(:, {'SubID','StimLoc','x','y','z'});
filename = fullfile(fileparts(pwd),'TMSvisualization','StimCoord_server_data.xlsx');
writetable(Data, filename)

