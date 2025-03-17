

% modified from similar code
% but using the data from aal with smaller size of ROIs (116 ROIs, previously 170 rois)

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
agg_func_conn_dir = fullfile(studydir, 'FuncConn_AAL_ROIs_PAID', 'agg');
if ~exist(agg_func_conn_dir,'dir')
    mkdir(agg_func_conn_dir)
end

all_corr_map_TMS = nan(num_rois,n_masks,2,nSubIDlist);
for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;
    func_conn_dir = fullfile(studydir, 'FuncConn_AAL_ROIs_PAID', SubID); 
    if include_shift
        mat_name = fullfile(func_conn_dir,'sub_corr_maps_w_shift_new.mat');
    else
        mat_name = fullfile(func_conn_dir,'sub_corr_maps_new.mat');
    end
    load(mat_name)
    all_corr_map_TMS(:,:,:,subj) = mean_corr_map_TMS;
end
avg_corr_map_TMS = squeeze(mean(all_corr_map_TMS,4));  % mean matrix averaged across subjects, for each TMS condition

%%
% plot the average corr map of sham and cTBS separately
% also the diff of the two
min_corr = min(tanh(avg_corr_map_TMS),[],'all');
max_corr = max(tanh(avg_corr_map_TMS),[],'all');

figure('Visible', 'off', 'Position', [100, 100, 1400, 1600]); % plot without showing
subplot(1,3,1)
imagesc(tanh(avg_corr_map_TMS(:,:,1)));
colorbar; 
clim([min_corr, max_corr]); 
title('FC (sham)');
xticks(1:n_masks);
xticklabels(mask_labels);  
yticks(1:num_rois);
yticklabels(name); 

subplot(1,3,2)
imagesc(tanh(avg_corr_map_TMS(:,:,2)));
colorbar; 
clim([min_corr, max_corr]); 
title('FC (cTBS)');
xticks(1:n_masks);
xticklabels(mask_labels);
yticks(1:num_rois);
yticklabels(name); 

diff_map = avg_corr_map_TMS(:,:,2) - avg_corr_map_TMS(:,:,1);
subplot(1,3,3)
imagesc(diff_map);
colorbar; 
title('FC (cTBS-sham)');
xticks(1:n_masks);
xticklabels(mask_labels); 
yticks(1:num_rois);
yticklabels(name); 

if include_shift
    img_name = fullfile(agg_func_conn_dir,'avg_corr_map_TMS_w_shift_new.bmp');
else
    img_name = fullfile(agg_func_conn_dir,'avg_corr_map_TMS_new.bmp');
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
            fprintf('with AAL roi %s\n',name{i})
            fprintf('p-value: %f',P)
            disp(STATS)
        end
        t_stat_matrix(i, j) = STATS.tstat;
    end
end
end   

% Plot the t-statistics as a heatmap
figure('Visible', 'off', 'Position', [100, 100, 2400, 1600]); % plot without showing
subplot(1,4,1)
imagesc(t_stat_matrix);
colorbar;
title('T-stat Heatmap (positive value means greater cTBS value)');
xlabel('');
ylabel('');
xticks(1:n_masks);
xticklabels(mask_labels);  
caxis([-3.5, 3.5]); 
yticks(1:num_rois);
yticklabels(name);  

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
            fprintf('with AAL roi %s\n',name{i})
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
xlabel('');
ylabel('');
xticks(1:n_masks);
xticklabels(mask_labels);  
caxis([-3.5, 3.5]); 
yticks(1:num_rois);
yticklabels(name);  

t_stat_matrix_aOFC = t_stat_matrix;

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
            fprintf('with AAL roi %s\n',name{i})
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
xlabel('');
ylabel('');
xticks(1:n_masks);
xticklabels(mask_labels);  
caxis([-3.5, 3.5]); 
yticks(1:num_rois);
yticklabels(name);  

t_stat_matrix_pOFC = t_stat_matrix;

% diff b/t aOFC and pOFC    
t_stat_matrix = NaN(num_rois, n_masks);

for j = 1:n_masks    
    fprintf('current mask %s\n\n',mask_labels{j})
for i = 1:num_rois
    
    % simply take the diff b/t cTBS and sham
    x_sham_dat_aOFC = squeeze(all_corr_map_TMS(i,j,1,aOFC_subs));
    x_cTBS_dat_aOFC = squeeze(all_corr_map_TMS(i,j,2,aOFC_subs));
    x_cTBS_vs_sham_aOFC = x_cTBS_dat_aOFC - x_sham_dat_aOFC;
    
    x_sham_dat_pOFC = squeeze(all_corr_map_TMS(i,j,1,pOFC_subs));
    x_cTBS_dat_pOFC = squeeze(all_corr_map_TMS(i,j,2,pOFC_subs));
    x_cTBS_vs_sham_pOFC = x_cTBS_dat_pOFC - x_sham_dat_pOFC;
    
    if ~all(isnan(x_sham_dat_aOFC))
        [H,P,CI,STATS] = ttest2(x_cTBS_vs_sham_aOFC,x_cTBS_vs_sham_pOFC);
        if P < 0.01
            fprintf('with AAL roi %s\n',name{i})
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
xlabel('');
ylabel('');
xticks(1:n_masks);
xticklabels(mask_labels);  
yticks(1:num_rois);
yticklabels(name);  

if include_shift
    img_name = fullfile(agg_func_conn_dir,'avg_corr_map_TMS_w_shift_new_stats.bmp');
else
    img_name = fullfile(agg_func_conn_dir,'avg_corr_map_TMS_new_stats.bmp');
end
saveas(gcf, img_name);

%% 
% want to find ROIs based on each of aOFC and pOFC results
% on aOFC-stim subjects, look for ROIs with lower FC values in cTBS vs. sham
% on pOFC-stim subjects, look for ROIs with lower FC values in cTBS vs. sham

dat_use = nan(num_rois,2);
dat_use(:,1) = t_stat_matrix_aOFC(:,1);
dat_use(:,2) = t_stat_matrix_pOFC(:,3);

% look for ROIs with low negative values of each column
thresh = 0; % can change to other more constraint threshold value (i.e. more negative)
idx_from_aOFC = find(dat_use(:,1) < thresh);
idx_from_pOFC = find(dat_use(:,2) < thresh);
idx_common = intersect(idx_from_aOFC,idx_from_pOFC); % here I'm using common idx, not trying to include more
name(idx_common)                                     % take a quick look at which ROIs are included
n_roi_use = length(idx_common);

% then for each subject, get the pattern vectors from the cTBS session
corr_patterns = squeeze(all_corr_map_TMS(idx_common,[1,3],2,:));
corr_patterns = reshape(corr_patterns, n_roi_use*2, nSubIDlist);
corr_patterns = corr_patterns';
corr_patterns = zscore(corr_patterns);
StimLoc_idx = 2 - cellfun(@(x) isequal(x, 'Posterior'), StimLoc); % aOFC: 2; pOFC: 1

% try if you can decode aOFC/pOFC based on those correlation pattern vectors
% use libsvm toolbox to do leave-one-subject-out CV
% corr_patterns: nsub * nfeatures (2*n_roi_use)
% StimLoc_idx: nsub

addpath('/Users/liuq13/libsvm-3.25/matlab')
options = '-s 0 -t 2 -c 1e-1 -g 0.5';     % C-SVC with RBF kernel, C=1, gamma=0.5
predicted_labels = zeros(nSubIDlist, 1);

for i = 1:nSubIDlist
    test_data = corr_patterns(i, :);
    test_label = StimLoc_idx(i);
    
    training_data = corr_patterns([1:i-1, i+1:end], :);
    training_labels = StimLoc_idx([1:i-1, i+1:end]);
    
    model = svmtrain(training_labels, training_data, options);
    predicted_labels(i) = svmpredict(test_label, test_data, model);
end

accuracy = sum(predicted_labels == StimLoc_idx) / nSubIDlist * 100;
fprintf('\nLOOCV Accuracy: %.2f%%\n', accuracy);