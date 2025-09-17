

% this code does post-hoc analysis on the motion parameters
% compare motion between cTBS and sham sessions

clc; clear all

%% setup
studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

MRIcount = readtable('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx',"ReadRowNames",true);

userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Library/CloudStorage/Box-Box/NODEAP_data_analysis');

%%

all_ctr_badvol = nan(nSubIDlist,n_rest_names); % # of bad vols across subs
all_fd = nan(nSubIDlist,n_rest_names); % total fd across vols 

for subj = 1:nSubIDlist

SubID = SubIDlist(subj).name;
SubDir = fullfile(MRIdir,SubID);
niidir = fullfile(SubDir, 'nifti');

% where to save nuisance regressors
NROutdir = fullfile(studydir,'NRegressor',SubID);

% load saved bad volume info
% 'badvol', 'ctr_badvol', 'ctr_NR'
load(fullfile(NROutdir,'BadVolumes.mat'));  
all_ctr_badvol(subj,:) = ctr_badvol;

for r = 1:n_rest_names % to loop over rest_names
    curr_rest = rest_names{r};
    if MRIcount{SubID,curr_rest}==0
        continue;
    end
    fprintf('%s of %s\n',curr_rest,SubID)
    
    if MRIcount{SubID,curr_rest}==0.5
        nscans = 205;
    else
        nscans = 310;
    end
    
    % find the motion parameters
    funcRunDir = fullfile(niidir,sprintf('%s_me',curr_rest));
    mpfilename = dir(fullfile(funcRunDir, 'rp_*.txt'));    %rp file contains 6 columns represnting positions  
    mp = load(fullfile(funcRunDir, mpfilename(1).name));          % mp motion parameters    
    mp_diff = [zeros(1,6); diff(mp)];
    
    translation_idx = [1,2,3]; % idx of translation difference parameters (x, y, z)
    rotation_idx = [4,5,6]; % idx of rotation difference parameters (roll, pitch, yaw) in radians
    
    head_radius = 50;  % Common assumption for head radius (in mm)

    % Translational FD (sum of absolute differences in x, y, z for each time point)
    translational_diff = sum(abs(mp_diff(:, translation_idx)), 2);  % Sum across the rows for translational parameters

    % Rotational FD (sum of absolute differences in roll, pitch, yaw, scaled by head radius)
    rotational_diff = head_radius * sum(abs(mp_diff(:, rotation_idx)), 2);  % Convert radians to mm and sum across the rows

    % Total FD is the sum of translational and rotational differences
    fd = translational_diff + rotational_diff;

    all_fd(subj,r) = sum(fd);
    
end



end

%%
% save this for later process
path = fullfile(HomeDir,'ProcessedData','motion_subs_sessions.mat');
save(path, 'all_ctr_badvol','all_fd'); 


