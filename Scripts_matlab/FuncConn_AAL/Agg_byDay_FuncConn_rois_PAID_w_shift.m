

%% aggregate correlation patterns from single sessions (without averaging)
% can use that to correlate with behavioral effects
% Day 1, Day 2, separately

clear; clc; close all

studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Library/CloudStorage/Box-Box/NODEAP_data_analysis');

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');
table_data = readtable(fullfile(HomeDir,'ProcessedData/SubConds.xlsx'));
StimLoc = table_data.StimLoc;

maskpath = fullfile(HomeDir,'Scripts_matlab/FuncConn_AAL');
aal_path = fullfile(maskpath,'TKsent','raal.nii');
aal_img = spm_vol(aal_path);
aal_data = spm_read_vols(aal_img);
num_rois = max(aal_data(:)); % Total number of ROIs
fprintf('Number of aal rois = %d\n',num_rois)

% load aal name file for labeling
aal_name_file = fullfile(maskpath,'TKsent','aal.mat');
load(aal_name_file)

% 4 TMS-related sphere ROIs
ROIpath = '/Volumes/X9Pro/NODEAP/ConnectivityMasks';
masks = {'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'
    };  
n_masks = length(masks);
mask_labels = {'aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'};

include_shift = true; % include temporal shifted gm values

%%
% aggregate correlation maps of each TMS condition
agg_func_conn_dir = fullfile(studydir, 'FuncConn_AAL_ROIs_PAID', 'agg_byDay');
if ~exist(agg_func_conn_dir,'dir')
    mkdir(agg_func_conn_dir)
end

%%
% read in correlation maps of each single session
rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);
dat = nan(num_rois,n_masks,n_rest_names,nSubIDlist);

for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;
    func_conn_dir = fullfile(studydir, 'FuncConn_AAL_ROIs_PAID', SubID); 
    for r = 1:n_rest_names
        
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,SubID)
        
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        
        if include_shift 
            mat_name = sprintf('conn_matrix_%s_w_shift_new.mat',curr_rest);
        else
            mat_name = sprintf('conn_matrix_%s_new.mat',curr_rest);
        end
        load(fullfile(func_conn_dir,mat_name))
        dat(:,:,r,subj) = correlation_matrix;
    end
end

%%
% save FC data for Day1 and Day2 separately
FCdat_Day2 = dat(:,:,[3,5,7],:);
FCdat_Day1 = dat(:,:,[2,4,6],:);

%% save FC data of each Day
path = '../../ProcessedData/FC_Day1_subs_sessions_w_shift_new.mat';
save(path, 'FCdat_Day1'); 

path = '../../ProcessedData/FC_Day2_subs_sessions_w_shift_new.mat';
save(path, 'FCdat_Day2');



