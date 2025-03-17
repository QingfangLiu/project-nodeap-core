

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
Order = table_data.StimOrder;
StimLoc = table_data.StimLoc;

nscans = 310; % 310 volumes for ME sequnce
nechos = 3;
TR = 1.5;
TEs = [14.60, 39.04, 63.48];
rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

maskpath = fullfile(HomeDir,'Scripts_matlab/FuncConn_AAL');
aal_path = fullfile(maskpath,'rAAL3v1.nii');
aal_img = spm_vol(aal_path);
aal_data = spm_read_vols(aal_img);
num_rois = max(aal_data(:)); % Total number of ROIs

%%
% 4 TMS-related sphere ROIs

ROIpath = '/Volumes/X9Pro/NODEAP/ConnectivityMasks';
masks = {'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'
    };  % masks in 3mm space
n_masks = length(masks);
mask_labels = {'aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'};

%%
% aggregate correlation maps of each TMS condition
all_corr_map_TMS = nan(num_rois,n_masks,2,nSubIDlist);

agg_func_conn_dir = fullfile(studydir, 'FuncConn_AAL_ROIs_PAID', 'agg');
if ~exist(agg_func_conn_dir,'dir')
    mkdir(agg_func_conn_dir)
end

for subj = 1:nSubIDlist
        
    SubID = SubIDlist(subj).name;
         
    func_conn_dir = fullfile(studydir, 'FuncConn_AAL_ROIs_PAID', SubID); 
    mat_name = fullfile(func_conn_dir,'sub_corr_maps.mat');
    load(mat_name)
    
    all_corr_map_TMS(:,:,:,subj) = mean_corr_map_TMS;
    
end
    
avg_corr_map_TMS = squeeze(mean(all_corr_map_TMS,4));

%%

fig = figure('Visible', 'off', 'Position', [100, 100, 600, 600]); % plot without showing
subplot(1,3,1)
imagesc(avg_corr_map_TMS(:,:,1));
colorbar; 
caxis([-1, 1]); 
title('FC (sham)');
xlabel('Masks');
ylabel('ROIs');

subplot(1,3,2)
imagesc(avg_corr_map_TMS(:,:,2));
colorbar; 
caxis([-1, 1]); 
title('FC (cTBS)');
xlabel('Masks');
ylabel('ROIs');

diff_map = avg_corr_map_TMS(:,:,2)-avg_corr_map_TMS(:,:,1);
subplot(1,3,3)
imagesc(diff_map);
colorbar; 
title('FC (cTBS-sham)');
xlabel('Masks');
ylabel('ROIs');

img_name = fullfile(agg_func_conn_dir,'avg_corr_map_TMS.bmp');
saveas(gcf, img_name);


%%
% simply use wilcoxon signed rank test to compare FC diff b/t sham and cTBS

x_sham = squeeze(all_corr_map_TMS(:,1,1,:));
x_cTBS = squeeze(all_corr_map_TMS(:,1,2,:));

for j = 1:n_masks
    
    fprintf('current mask %s\n\n',mask_labels{j})
    
for i = 1:num_rois
    
    x_sham_dat = x_sham(i,:);
    x_cTBS_dat = x_cTBS(i,:);
    
    if ~all(isnan(x_sham_dat))
        [p, h, stats] = signrank(x_sham_dat, x_cTBS_dat);
        if p < 0.05
            disp(i)
            disp(p)
            disp(stats)
        end
    end
    
end
end   

% found some trending difference, but no one would pass the multiple
% correction threshold


    
    