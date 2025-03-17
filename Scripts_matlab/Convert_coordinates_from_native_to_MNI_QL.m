
% This code is adapted from the same script used for SPEEDTMS,
% but here for NODEAP
% difference is that each person has only one set of stimulation
% coordinates, either targeting aOFC or pOFC

clc; clear; close all

HomeDir = '/Users/liuq13/Documents/Projects/NODEAP_data_analysis';
studydir = '/Volumes/X9Pro/NODEAP/MRI';

SubIDlist = dir(fullfile(studydir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

file_path = fullfile(HomeDir,'ProcessedData','SubConds.xlsx');
table_data = readtable(file_path); 
Subs = table_data.SubID; % this matches SubIDlist
StimLoc = table_data.StimLoc;

%% load TMS coordinates in subjects' native space

StimCoord = zeros(nSubIDlist,3);
for subj = 1:nSubIDlist
    sn = SubIDlist(subj).name;
    subdir = fullfile(studydir,sn);
    load(fullfile(subdir,'stimulation_coordinates.mat'))
    if strcmp(StimLoc(subj),'Posterior')
        StimCoord(subj,:) = regions(2).StimCoord;
    else
        StimCoord(subj,:) = regions(1).StimCoord;
    end
end

varNames = {'x','y','z'};
Data = array2table(StimCoord, 'VariableNames', varNames);

% add subject ID
Data.SubID = Subs;
Data.StimLoc = StimLoc;

Data = Data(:, {'SubID','StimLoc','x','y','z'});
filename = fullfile(fileparts(pwd),'TMSvisualization','StimCoord.xlsx');
writetable(Data, filename)

%% convert coordinates from subjects' native space to MNI space

for subj = 1:nSubIDlist
    
    sn = SubIDlist(subj).name;
    subdir = fullfile(studydir,sn);
 
    % find the .mat file and load it
    mat_path = fullfile(subdir,'nii_for_TMS/anat','*.mat');
    mat_file = dir(mat_path);
    mat_file_name = fullfile(mat_file.folder,mat_file.name);
    load(mat_file_name)

    % convert coordinates to MNI space 
    % 'Affine' comes from loading the mat file
    MNICoord(subj,:) = Affine(1:3,:) * [StimCoord(subj,:) 1]';
    
end

%% Organize coordinates into a table and save as txt file

dat = MNICoord;
dat = round(dat * 100) / 100;
dat(:,4) = 2 - strcmp(StimLoc,'Posterior');  % code for color
dat(:,5) = 1;  % code for size
dat = array2table(dat);
dat(:,6) = StimLoc;

% rename variable name and save
dat.Properties.VariableNames = ["x","y","z","color","size","label"];
filename = fullfile(HomeDir,'TMSvisualization','myNode.txt');
writetable(dat,filename,'Delimiter',' ')  

%% Use Daria's coordinates 
% with slight discrepancy than what I computed

filename = fullfile(HomeDir,'TMSvisualization','StimCoord_Daria.xlsx');
StimCoord = readtable(filename);
StimCoord = table2array(StimCoord(:,[3,4,5]));

for subj = 1:nSubIDlist
    
    sn = SubIDlist(subj).name;
    subdir = fullfile(studydir,sn);
 
    % find the .mat file and load it
    mat_path = fullfile(subdir,'nii_for_TMS/anat','*.mat');
    mat_file = dir(mat_path);
    mat_file_name = fullfile(mat_file.folder,mat_file.name);
    load(mat_file_name)

    % convert coordinates to MNI space 
    % 'Affine' comes from loading the mat file
    MNICoord(subj,:) = Affine(1:3,:) * [StimCoord(subj,:) 1]';
    
end

%% Organize coordinates into a table and save as txt file

dat = MNICoord;
dat = round(dat * 100) / 100;
dat(:,4) = 2 - strcmp(StimLoc,'Posterior');  % code for color (P=1, A=2)
dat(:,5) = 1;  % code for size
dat = array2table(dat);
dat(:,6) = strrep(Subs, 'NODEAP_', '');

% rename variable name and save
dat.Properties.VariableNames = ["x","y","z","color","size","label"];
filename = fullfile(HomeDir,'TMSvisualization','myNode_Daria.txt');
writetable(dat,filename,'Delimiter',' ')  





