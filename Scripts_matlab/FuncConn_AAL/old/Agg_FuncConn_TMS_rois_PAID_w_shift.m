

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
aal_path = fullfile(maskpath,'rAAL3v1.nii');
aal_img = spm_vol(aal_path);
aal_data = spm_read_vols(aal_img);

num_rois = max(aal_data(:)); % Total number of ROIs

%%
% 4 TMS-related sphere ROIs in 2mm space
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
all_corr_map_TMS = nan(num_rois,n_masks,2,nSubIDlist);

agg_func_conn_dir = fullfile(studydir, 'FuncConn_AAL_ROIs_PAID', 'agg');
if ~exist(agg_func_conn_dir,'dir')
    mkdir(agg_func_conn_dir)
end

for subj = 1:nSubIDlist
        
    SubID = SubIDlist(subj).name;
    func_conn_dir = fullfile(studydir, 'FuncConn_AAL_ROIs_PAID', SubID);    
    if include_shift
        mat_name = fullfile(func_conn_dir,'sub_corr_maps_w_shift.mat');
    else
        mat_name = fullfile(func_conn_dir,'sub_corr_maps.mat');
    end
    load(mat_name)
    all_corr_map_TMS(:,:,:,subj) = mean_corr_map_TMS;
end
    
avg_corr_map_TMS = squeeze(mean(all_corr_map_TMS,4));

%%
min_corr = min(tanh(avg_corr_map_TMS),[],'all');
max_corr = max(tanh(avg_corr_map_TMS),[],'all');

figure('Visible', 'off', 'Position', [100, 100, 1200, 600]); % plot without showing
subplot(1,3,1)
imagesc(tanh(avg_corr_map_TMS(:,:,1)));
colorbar; 
caxis([min_corr, max_corr]); 
title('FC (sham)');
xlabel('Masks');
ylabel('ROIs');
xticks(1:n_masks);
xticklabels(mask_labels);  

subplot(1,3,2)
imagesc(tanh(avg_corr_map_TMS(:,:,2)));
colorbar; 
caxis([min_corr, max_corr]); 
title('FC (cTBS)');
xlabel('Masks');
ylabel('ROIs');
xticks(1:n_masks);
xticklabels(mask_labels);  

diff_map = avg_corr_map_TMS(:,:,2) - avg_corr_map_TMS(:,:,1);
subplot(1,3,3)
imagesc(diff_map);
colorbar; 
title('FC (cTBS-sham)');
xlabel('Masks');
ylabel('ROIs');
xticks(1:n_masks);
xticklabels(mask_labels);  

if include_shift
    img_name = fullfile(agg_func_conn_dir,'avg_corr_map_TMS_w_shift.bmp');
else
    img_name = fullfile(agg_func_conn_dir,'avg_corr_map_TMS.bmp');
end
saveas(gcf, img_name);


%%
% simply use wilcoxon signed rank test to compare FC diff b/t sham and cTBS
clc
disp('Using all 48 subjects')
t_stat_matrix = NaN(num_rois, n_masks);

for j = 1:n_masks
    fprintf('current mask %s\n\n',mask_labels{j})
for i = 1:num_rois
    x_sham_dat = squeeze(all_corr_map_TMS(i,j,1,:));
    x_cTBS_dat = squeeze(all_corr_map_TMS(i,j,2,:));
    
    if ~all(isnan(x_sham_dat))
        [H,P,CI,STATS] = ttest(x_cTBS_dat,x_sham_dat);
        if P < 0.01
            fprintf('with AAL roi %d\n',i)
            fprintf('p-value: %f',P)
            disp(STATS)
        end
        t_stat_matrix(i, j) = STATS.tstat;
    end
end
end   

% Plot the t-statistics as a heatmap
figure('Visible', 'off', 'Position', [100, 100, 1600, 600]); % plot without showing
subplot(1,4,1)
imagesc(t_stat_matrix);
colorbar;
title('T-stat Heatmap (positive value means greater cTBS value)');
xlabel('Masks');
ylabel('ROIs');
xticks(1:n_masks);
xticklabels(mask_labels);  
caxis([-3.5, 3.5]); 

% using only aOFC subjects
aOFC_subs = cellfun(@(x) isequal(x, 'Anterior'), StimLoc);
disp('Using only aOFC subjects')
t_stat_matrix = NaN(num_rois, n_masks);

for j = 1:n_masks    
    fprintf('current mask %s\n\n',mask_labels{j})
for i = 1:num_rois
    x_sham_dat = squeeze(all_corr_map_TMS(i,j,1,aOFC_subs));
    x_cTBS_dat = squeeze(all_corr_map_TMS(i,j,2,aOFC_subs));
    
    if ~all(isnan(x_sham_dat))
        [H,P,CI,STATS] = ttest(x_cTBS_dat,x_sham_dat);
        if P < 0.01
            fprintf('with AAL roi %d\n',i)
            fprintf('p-value: %f',P)
            disp(STATS)
        end
        t_stat_matrix(i, j) = STATS.tstat;
    end
    
end
end   

subplot(1,4,2)
imagesc(t_stat_matrix);
colorbar;
title('aOFC sub only');
xlabel('Masks');
ylabel('ROIs');
xticks(1:n_masks);
xticklabels(mask_labels);  
caxis([-3.5, 3.5]); 

% using only pOFC subjects
pOFC_subs = cellfun(@(x) isequal(x, 'Posterior'), StimLoc);
disp('Using only pOFC subjects')
t_stat_matrix = NaN(num_rois, n_masks);

for j = 1:n_masks    
    fprintf('current mask %s\n\n',mask_labels{j})
for i = 1:num_rois
    x_sham_dat = squeeze(all_corr_map_TMS(i,j,1,pOFC_subs));
    x_cTBS_dat = squeeze(all_corr_map_TMS(i,j,2,pOFC_subs));
    
    if ~all(isnan(x_sham_dat))
        [H,P,CI,STATS] = ttest(x_cTBS_dat,x_sham_dat);
        if P < 0.01
            fprintf('with AAL roi %d\n',i)
            fprintf('p-value: %f',P)
            disp(STATS)
        end
        t_stat_matrix(i, j) = STATS.tstat;
    end
    
end
end   

subplot(1,4,3)
imagesc(t_stat_matrix);
colorbar;
title('pOFC sub only');
xlabel('Masks');
ylabel('ROIs');
xticks(1:n_masks);
xticklabels(mask_labels);  
caxis([-3.5, 3.5]); 


% diff b/t aOFC and pOFC    
t_stat_matrix = NaN(num_rois, n_masks);

for j = 1:n_masks    
    fprintf('current mask %s\n\n',mask_labels{j})
for i = 1:num_rois
    
    x_sham_dat_aOFC = squeeze(all_corr_map_TMS(i,j,1,aOFC_subs));
    x_cTBS_dat_aOFC = squeeze(all_corr_map_TMS(i,j,2,aOFC_subs));
    x_cTBS_vs_sham_aOFC = x_cTBS_dat_aOFC - x_sham_dat_aOFC;
    
    x_sham_dat_pOFC = squeeze(all_corr_map_TMS(i,j,1,pOFC_subs));
    x_cTBS_dat_pOFC = squeeze(all_corr_map_TMS(i,j,2,pOFC_subs));
    x_cTBS_vs_sham_pOFC = x_cTBS_dat_pOFC - x_sham_dat_pOFC;
    
    if ~all(isnan(x_sham_dat_aOFC))
        [H,P,CI,STATS] = ttest2(x_cTBS_vs_sham_aOFC,x_cTBS_vs_sham_pOFC);
        if P < 0.01
            fprintf('with AAL roi %d\n',i)
            fprintf('p-value: %f',P)
            disp(STATS)
        end
        t_stat_matrix(i, j) = STATS.tstat;
    end
    
end
end   
    
subplot(1,4,4)
imagesc(t_stat_matrix);
colorbar;
title('aOFC - pOFC (cTBS-sham)');
xlabel('Masks');
ylabel('ROIs');
xticks(1:n_masks);
xticklabels(mask_labels);  


if include_shift
    img_name = fullfile(agg_func_conn_dir,'avg_corr_map_TMS_w_shift_stats.bmp');
else
    img_name = fullfile(agg_func_conn_dir,'avg_corr_map_TMS_stats.bmp');
end
saveas(gcf, img_name);


