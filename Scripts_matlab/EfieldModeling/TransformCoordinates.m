

% this code uses the subject2mni_coords function from SimNIBS to convert
% coordinates from native space to MNI space
% via the individual head models calculated by SimNIBS
% results show some deviation from those calculated from SPM (another
% folder) but visually not distinguishable

clear; clc; close all

userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Library/CloudStorage/Box-Box/NODEAP_data_analysis');
studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
EfieldDir = fullfile(studydir,'EfieldModeling');

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);
ROIdir = fullfile(studydir,'ConnectivityMasks');

% add SimNIBS toolbox
addpath('/Users/liuq13/Applications/SimNIBS-4.1/matlab_tools')

StimCoordFile = fullfile(HomeDir,'TMSvisualization','StimCoord_Daria.xlsx');
StimCoord = readtable(StimCoordFile);

file_path = fullfile(HomeDir,'ProcessedData','SubConds.xlsx');
table_data = readtable(file_path); 
Subs = table_data.SubID; % this matches SubIDlist
StimLoc = table_data.StimLoc;

%%
MNICoord = nan(nSubIDlist,3);
for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;
    fprintf('%s\n',SubID)
    
    coords_mni = table2array(StimCoord(subj,[3,4,5]));
    subdir = fullfile(EfieldDir,SubID,sprintf('m2m_%s',SubID));
    transformation_type = 'nonl';
    MNICoord(subj,:) = subject2mni_coords(coords_mni, subdir, transformation_type);
    
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
filename = fullfile(HomeDir,'Scripts_matlab','EfieldModeling','myNode_Daria_SimNIBS.txt');
writetable(dat,filename,'Delimiter',' ')  
