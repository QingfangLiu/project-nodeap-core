

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
% aggregate correlation maps of each TMS condition
all_corr_map_TMS = nan(num_rois,num_rois,2,nSubIDlist);

agg_func_conn_dir = fullfile(studydir, 'FuncConn_AAL_PAID', 'agg');
if ~exist(agg_func_conn_dir,'dir')
    mkdir(agg_func_conn_dir)
end

for subj = 1:nSubIDlist
        
    SubID = SubIDlist(subj).name;
         
    func_conn_dir = fullfile(studydir, 'FuncConn_AAL_PAID', SubID); 
    mat_name = fullfile(func_conn_dir,'sub_corr_maps.mat');
    load(mat_name)
    
    all_corr_map_TMS(:,:,:,subj) = mean_corr_map_TMS;
    
end
    
avg_corr_map_TMS = squeeze(mean(all_corr_map_TMS,4));

%%

fig = figure('Visible', 'off', 'Position', [100, 100, 1600, 400]); % plot without showing
subplot(1,3,1)
imagesc(avg_corr_map_TMS(:,:,1));
colorbar; 
caxis([-1, 1]); 
title('ROI-based Functional Connectivity Matrix (sham)');
xlabel('ROIs');
ylabel('ROIs');

subplot(1,3,2)
imagesc(avg_corr_map_TMS(:,:,2));
colorbar; 
caxis([-1, 1]); 
title('ROI-based Functional Connectivity Matrix (cTBS)');
xlabel('ROIs');
ylabel('ROIs');

diff_map = avg_corr_map_TMS(:,:,2)-avg_corr_map_TMS(:,:,1);
subplot(1,3,3)
imagesc(diff_map);
colorbar; 
title('ROI-based Functional Connectivity Matrix (cTBS-sham)');
xlabel('ROIs');
ylabel('ROIs');

img_name = fullfile(agg_func_conn_dir,'avg_corr_map_TMS.bmp');
saveas(gcf, img_name);


%%

% simply use wilcoxon signed rank test to compare FC diff b/t sham and cTBS
% should instead use data from each session of each subject

% ROI-28: OFC ant-R
x_sham = squeeze(all_corr_map_TMS(28,:,1,:));
x_cTBS = squeeze(all_corr_map_TMS(28,:,2,:));

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


%%
% ROI-30: OFC post-R

x_sham = squeeze(all_corr_map_TMS(30,:,1,:));
x_cTBS = squeeze(all_corr_map_TMS(30,:,2,:));

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


% for either ROI 28 or 30, there's some trending difference, but won't pass
% the multiple comparison threshold. 
% so not continue in this direction

%%

subs_corr_vals = nan(n_rest_names,num_rois,nSubIDlist);
use_roi = 30;

for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;   
    func_conn_dir = fullfile(studydir, 'FuncConn_AAL_PAID', SubID); 
    
    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,SubID)
            
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        
        mat_name = fullfile(func_conn_dir,sprintf('conn_matrix_%s.mat',curr_rest));
        load(mat_name)
        subs_corr_vals(r,:,subj) = correlation_matrix(use_roi,:);
        
    end
end

% save this for later process
path = fullfile(HomeDir,'ProcessedData',sprintf('subs_corr_vals_sessions_use_roi%d.mat',use_roi));
save(path, 'subs_corr_vals'); 






    
    
    