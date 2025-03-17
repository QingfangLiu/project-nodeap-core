
% compute conn maps within the OFC and LPFC clusters, using each voxel as a
% seed
% use 3mm MNI space
% use filtered dat in global conn analysis

% this one calculates map with all gm voxels

clear; clc; close all
studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis');

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

table_data = readtable(fullfile(HomeDir,'ProcessedData/SubConds.xlsx'));
Order = table_data.StimOrder;

nscans = 310; % 310 volumes for ME sequnce
nechos = 3;
TR = 1.5;
TEs = [14.60, 39.04, 63.48];
rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

%%
maskpath = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis/Scripts_matlab/GlobalConn/mask');
masks = {'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'
    };  % masks in 3mm space
n_masks = length(masks);
mask_labels = {'aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'};

gm_nii = fullfile(maskpath,'gm_0.1_3mm.nii'); % gray matter mask
gm_dat = spm_read_vols(spm_vol(gm_nii));
gm_idx = find(gm_dat > 0);
n_gm_vox = length(gm_idx);

%% concatenante the ROIs
voxel_idx = [];
roi_idx = [];
for i = 1:n_masks
    [targetvol, ~] = spm_read_vols(spm_vol_nifti(fullfile(maskpath,masks{i})));
    tmploc = find(targetvol==1);  % loc idx in 53*63*52 space
    voxel_idx = [voxel_idx;tmploc];
    roi_idx = [roi_idx;repmat(i,length(tmploc),1)];
end

array_size = size(gm_dat);  % 53*63*52
[check_exist, roi_gm_linear_indices] = ismember(voxel_idx, gm_idx);
all(check_exist)
% check_exist: make sure all ROI indices are within the gray matter (as it
% should)
% roi_gm_linear_indices: indices of ROIs that are in the gray matter

all_idx = [voxel_idx,roi_gm_linear_indices,roi_idx];
nvox = size(all_idx,1);

% col1: voxel in 3mm MNI space
% col2: voxel in gray matter linear index
% col3: which roi: 1,2,3,4

%%
for subj = 1:nSubIDlist

    sn = SubIDlist(subj).name;
    
for r = 1:n_rest_names % to loop over rest_names
        
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,sn)
        
        func_conn_dir = fullfile(studydir, 'ConnMap_PAID_global', sn); % path to save the GLM result from PAID result
        if ~exist(func_conn_dir,'dir')
            mkdir(func_conn_dir)
        end
        
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        
        load(fullfile(studydir, 'GlobalConn_PAID', sn, curr_rest, 'tc_filtered_3mm.mat')); 
        
        % reduce the 'dat' to voxels we care about
        % dat = dat(:,all_idx(:,2));
        dat = zscore(dat); % zcore for convenient Pearson correlation calculation
                
%%
    corr_map = nan(nvox,n_gm_vox);
    parfor i = 1:nvox
        use_idx = all_idx(i,2);
        corr_map(i,:) = (dat(:,use_idx)' * dat)./ nscans; % Pearson correlation on z-scored values (with each voxel)
        if ~mod(i,50)
            fprintf('\t\t%d of %d\n', i, nvox);
        end
    end
    
    fig = figure('Visible', 'off', 'Position', [100, 100, 2000, 500]); % plot without showing
    imagesc(corr_map);
    colorbar; % Add a colorbar to show the scale
    caxis([-1, 1]); % Adjust these limits based on your data
    title(sprintf('Voxel-voxel connectivity map of %s',curr_rest));
    xlabel('');
    ylabel('');
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
    colormap(jet); 
    
    % Determine the boundaries for each group
    group_boundaries = find(diff(all_idx(:,3))) + 1;
    tick_positions = [mean([1, group_boundaries(1)]), ...
                  mean([group_boundaries(1), group_boundaries(2)]), ...
                  mean([group_boundaries(2), group_boundaries(3)]), ...
                  mean([group_boundaries(3), 349])];

    
    set(gca, ...
         'YTick', tick_positions, 'YTickLabel', mask_labels, ...
         'TickLength', [0 0], 'FontSize', 12, 'FontWeight', 'bold', 'XColor', 'k', 'YColor', 'k');
    
    img_name = fullfile(func_conn_dir,sprintf('conn_map_%s.bmp',curr_rest));
    saveas(gcf, img_name);
    
    % save this corr mat 
    mat_name = fullfile(func_conn_dir,sprintf('conn_map_%s.mat',curr_rest));
    save(mat_name,'corr_map')
       
    
end
end

%%
% to calculate the mean of the conn maps

for subj = 1:nSubIDlist
    sn = SubIDlist(subj).name;
    func_conn_dir = fullfile(studydir, 'ConnMap_PAID_global', sn); 
    all_corr_map = cell(n_rest_names,1);
    
    % coding sham and cTBS as 1 and 2, respectively
    % null as 3
    tmp_order = Order(subj);
    switch tmp_order
        case 123
            TMS_types = [3,2,1,1,2,1,1];            
        case 132
            TMS_types = [3,2,1,1,1,1,2];       
        case 213
            TMS_types = [3,1,2,2,1,1,1];            
        case 231
            TMS_types = [3,1,2,1,1,2,1];            
        case 312
            TMS_types = [3,1,1,2,1,1,2];               
        case 321
            TMS_types = [3,1,1,1,2,2,1];
    end
    TMS_types(count_table(subj,:)<1) = nan;

    for r = 1:n_rest_names 
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,sn)
        if count_table(subj,r)<1
            disp('Not enough data, put nan values ...')
            all_corr_map{r,1} = nan(nvox,n_gm_vox);
            continue;
        end
        mat_name = fullfile(func_conn_dir,sprintf('conn_map_%s.mat',curr_rest));
        load(mat_name)
        all_corr_map{r,1} = corr_map;
    end
    
    indices_sham = find(TMS_types == 1);  
    indices_cTBS = find(TMS_types == 2); 
        
    all_maps_cat = cat(3, all_corr_map{:});
    mean_corr_map = squeeze(nanmean(all_maps_cat,2)); % avg across gm voxels
    
    mean_corr_map_TMS(:,1) = nanmean(mean_corr_map(:,indices_sham),2);
    mean_corr_map_TMS(:,2) = nanmean(mean_corr_map(:,indices_cTBS),2);
    
    abs_mean_corr_map = squeeze(nanmean(abs(all_maps_cat),2)); % avg abs map across gm voxels
    
    abs_mean_corr_map_TMS(:,1) = nanmean(abs_mean_corr_map(:,indices_sham),2);
    abs_mean_corr_map_TMS(:,2) = nanmean(abs_mean_corr_map(:,indices_cTBS),2);
    
    % save both mean corr maps 
    mat_name = fullfile(func_conn_dir,'mean_corr_map_TMS.mat');
    save(mat_name,'mean_corr_map_TMS','abs_mean_corr_map_TMS')
    
    % later added: also save the mean map across sessions for each subject
    % so later I can check across-subj consistency of the FC patterns
    sub_corr_map = squeeze(nanmean(all_maps_cat,3)); % avg across sessions
    
    fig = figure('Visible', 'off', 'Position', [100, 100, 2000, 500]); % plot without showing
    imagesc(corr_map);
    colorbar; % Add a colorbar to show the scale
    caxis([-1, 1]); % Adjust these limits based on your data
    title('Voxel-voxel connectivity map');
    xlabel('');
    ylabel('');
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
    colormap(jet); 
    
    % Determine the boundaries for each group
    group_boundaries = find(diff(all_idx(:,3))) + 1;
    tick_positions = [mean([1, group_boundaries(1)]), ...
                  mean([group_boundaries(1), group_boundaries(2)]), ...
                  mean([group_boundaries(2), group_boundaries(3)]), ...
                  mean([group_boundaries(3), 349])];

    
    set(gca, ...
         'YTick', tick_positions, 'YTickLabel', mask_labels, ...
         'TickLength', [0 0], 'FontSize', 12, 'FontWeight', 'bold', 'XColor', 'k', 'YColor', 'k');
    
    img_name = fullfile(func_conn_dir,'sub_corr_map.bmp');
    saveas(gcf, img_name);
    mat_name = fullfile(func_conn_dir,'sub_corr_map.mat');
    save(mat_name,'sub_corr_map')
    
end


%% organize corr map across subjects
for subj = 1:nSubIDlist
    sn = SubIDlist(subj).name;
    func_conn_dir = fullfile(studydir, 'ConnMap_PAID_global', sn);
    mat_name = fullfile(func_conn_dir,'mean_corr_map_TMS.mat');
    load(mat_name)
    for k = 1:4
        loc = all_idx(:,3)==k;
        subs_corr_map_sham(subj,k) = mean(mean_corr_map_TMS(loc,1));
        subs_corr_map_cTBS(subj,k) = mean(mean_corr_map_TMS(loc,2));
        abs_subs_corr_map_sham(subj,k) = mean(abs_mean_corr_map_TMS(loc,1));
        abs_subs_corr_map_cTBS(subj,k) = mean(abs_mean_corr_map_TMS(loc,2));
    end
end

agg_dir = fullfile(studydir, 'ConnMap_PAID_global','Agg');
writematrix(subs_corr_map_sham, fullfile(agg_dir,'subs_corr_map_sham.csv'));
writematrix(subs_corr_map_cTBS, fullfile(agg_dir,'subs_corr_map_cTBS.csv'));
writematrix(abs_subs_corr_map_sham, fullfile(agg_dir,'abs_subs_corr_map_sham.csv'));
writematrix(abs_subs_corr_map_cTBS, fullfile(agg_dir,'abs_subs_corr_map_cTBS.csv'));

%%
% to look at FC reliability pattern
% by calculating correlations between each session and the average of else
% do this for each voxel of interest

for subj = 1:nSubIDlist
    sn = SubIDlist(subj).name;
    func_conn_dir = fullfile(studydir, 'ConnMap_PAID_global', sn); 
    all_corr_map = cell(n_rest_names,1);
    
    % coding sham and cTBS as 1 and 2, respectively
    % null as 3
    tmp_order = Order(subj);
    switch tmp_order
        case 123
            TMS_types = [3,2,1,1,2,1,1];            
        case 132
            TMS_types = [3,2,1,1,1,1,2];       
        case 213
            TMS_types = [3,1,2,2,1,1,1];            
        case 231
            TMS_types = [3,1,2,1,1,2,1];            
        case 312
            TMS_types = [3,1,1,2,1,1,2];               
        case 321
            TMS_types = [3,1,1,1,2,2,1];
    end
    TMS_types(count_table(subj,:)<1) = nan;

    for r = 1:n_rest_names 
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,sn)
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        mat_name = fullfile(func_conn_dir,sprintf('conn_map_%s.mat',curr_rest));
        load(mat_name)
        all_corr_map{r,1} = corr_map;
    end

    correlation = nan(nvox,n_rest_names);
    for r = 1:n_rest_names  
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,sn)
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        % map of this run
        curr_map = all_corr_map{r,1};
        % mean map of other runs
        all_other_maps = all_corr_map;
        all_other_maps{r} = [];
        all_other_maps_cat = cat(3, all_other_maps{:});
        mean_corr_map = nanmean(all_other_maps_cat,3);

        parfor i = 1:nvox
            correlation(i,r) = corr(mean_corr_map(i, :)', curr_map(i, :)');
        end
    end

    avg_correlation = nanmean(correlation,2);
    
    % agg based on sham or cTBS sessions
    avg_correlation_cTBS = nanmean(correlation(:,TMS_types==2),2);
    avg_correlation_sham = nanmean(correlation(:,TMS_types==1),2);

    % save voxel-wise FC reliability as a mat file
    mat_name = fullfile(func_conn_dir,'FC_reliability.mat');
    save(mat_name,'correlation','avg_correlation','all_idx',...
    'avg_correlation_cTBS','avg_correlation_sham')

end

% organize func conn reliability together
for subj = 1:nSubIDlist
    sn = SubIDlist(subj).name;
    func_conn_dir = fullfile(studydir, 'ConnMap_PAID_global', sn);
    mat_name = fullfile(func_conn_dir,'FC_reliability.mat');
    load(mat_name)
    for k = 1:4
        subs_correlation(subj,k) = mean(avg_correlation(all_idx(:,3)==k));
        subs_correlation_cTBS(subj,k) = mean(avg_correlation_cTBS(all_idx(:,3)==k));
        subs_correlation_sham(subj,k) = mean(avg_correlation_sham(all_idx(:,3)==k));
    end
end

figure;
for i = 1:4
    subplot(2, 2, i);
    histogram(subs_correlation(:, i));
    title(sprintf('Histogram of %s', mask_labels{i}));
    xlabel('Functional conn reliability');
    ylabel('Frequency');
end

agg_dir = fullfile(studydir, 'ConnMap_PAID_global','Agg');
if ~exist(agg_dir,'dir')
    mkdir(agg_dir)
end
img_name = fullfile(agg_dir,'FuncConn_reliability.bmp');
saveas(gcf, img_name);
writematrix(subs_correlation, fullfile(agg_dir,'subs_correlation.csv'));
writematrix(subs_correlation_cTBS, fullfile(agg_dir,'subs_correlation_cTBS.csv'));
writematrix(subs_correlation_sham, fullfile(agg_dir,'subs_correlation_sham.csv'));

%%
% want to compare cTBS vs sham in terms of func conn reliability
% for each subject, from 4 sham sessions, sample 2 sessions; also 2 cTBS
% sessions, one Day 0 session
% simi to bootstrap idea

for subj = 1:nSubIDlist
    sn = SubIDlist(subj).name;
    func_conn_dir = fullfile(studydir, 'ConnMap_PAID_global', sn); 
    
    % coding sham and cTBS as 1 and 2, respectively
    % null as 3
    tmp_order = Order(subj);
    switch tmp_order
        case 123
            TMS_types = [3,2,1,1,2,1,1];            
        case 132
            TMS_types = [3,2,1,1,1,1,2];       
        case 213
            TMS_types = [3,1,2,2,1,1,1];            
        case 231
            TMS_types = [3,1,2,1,1,2,1];            
        case 312
            TMS_types = [3,1,1,2,1,1,2];               
        case 321
            TMS_types = [3,1,1,1,2,2,1];
    end
    
    all_corr_map = [];
    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,sn)
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        mat_name = fullfile(func_conn_dir,sprintf('conn_map_%s.mat',curr_rest));
        load(mat_name)
        all_corr_map{r,1} = corr_map;
    end
    
    TMS_types(count_table(subj,:)<1) = nan;
    
    if length(find(TMS_types == 2))<2
        disp('Not enough number of cTBS sessions, use one session to compute mean')
        n_sample = 1;
    else
        n_sample = 2;
    end
    
    niter = 1e2;
    parfor iter = 1:niter
        
        fprintf('Working on iter %d\n',iter)
        indices_sham = find(TMS_types == 1);  
        sampled_indices_sham = randsample(indices_sham, n_sample);
        indices_cTBS = find(TMS_types == 2);  
        sampled_indices_cTBS = randsample(indices_cTBS, n_sample);
        
        if count_table(subj,1)<1
            indices_use = [sampled_indices_cTBS,sampled_indices_sham];
        else
            indices_use = [sampled_indices_cTBS,sampled_indices_sham,1]; % all session indices to use
        end
        
        use_corr_maps = all_corr_map(indices_use);
        all_corr_maps = cat(3, use_corr_maps{:});
        avg_use_map = mean(all_corr_maps, 3);
    
    % compute corr with sham runs
    correlation_sham = nan(nvox,n_sample);
    for k = 1:n_sample
        tmp_map = all_corr_map{sampled_indices_sham(k)};
        for i = 1:nvox
            correlation_sham(i,k) = corr(avg_use_map(i, :)', tmp_map(i, :)');
        end
    end
    corr_sham = mean(correlation_sham,2);
    
    % compute corr with cTBS runs
    correlation_cTBS = nan(nvox,n_sample);
    for k = 1:n_sample
        tmp_map = all_corr_map{indices_cTBS(k)};
        for i = 1:nvox
            correlation_cTBS(i,k) = corr(avg_use_map(i, :)', tmp_map(i, :)');
        end
    end
    corr_cTBS = mean(correlation_cTBS,2);
    out_corr{iter} = [corr_sham,corr_cTBS];
    end
    
    % avg across iterations
    all_out_corr = cat(3, out_corr{:});
    avg_out_corr = mean(all_out_corr, 3);
    
    % save this corr mat 
    mat_name = fullfile(func_conn_dir,'FC_corr_sham_cTBS.mat');
    save(mat_name,'avg_out_corr','TMS_types','all_idx')
    
end

%%
% organize func conn sham vs cTBS together
for subj = 1:nSubIDlist
    sn = SubIDlist(subj).name;
    func_conn_dir = fullfile(studydir, 'ConnMap_PAID_global', sn);
    mat_name = fullfile(func_conn_dir,'FC_corr_sham_cTBS.mat');
    load(mat_name)
    for k = 1:4
        subs_correlation_sham(subj,k) = mean(avg_out_corr(all_idx(:,3)==k,1));
        subs_correlation_cTBS(subj,k) = mean(avg_out_corr(all_idx(:,3)==k,2));
    end
end

agg_dir = fullfile(studydir, 'ConnMap_PAID_global','Agg');
writematrix(subs_correlation_cTBS, fullfile(agg_dir,'subs_correlation_cTBS_bs.csv'));
writematrix(subs_correlation_sham, fullfile(agg_dir,'subs_correlation_sham_bs.csv'));



%% look at across-subject consistency 
% using pairwise correlation on average subject-wise corr map
% load mean pattern from each subject
agg_corr_map = nan(nSubIDlist,n_gm_vox,n_masks);
for subj = 1:nSubIDlist
    sn = SubIDlist(subj).name;
    func_conn_dir = fullfile(studydir, 'ConnMap_PAID_global', sn);
    mat_name = fullfile(func_conn_dir,'sub_corr_map.mat');
    load(mat_name)
    for k = 1:n_masks
        loc = all_idx(:,3)==k;
        agg_corr_map(subj,:,k) = nanmean(sub_corr_map(loc,:),1);
    end
end

% compute and plot pairwise correlation
for k = 1:n_masks
    
conn_matrix = agg_corr_map(:,:,k);
cor_matrix = corr(conn_matrix');

% Extract the upper triangular part of the matrix, excluding the diagonal
upper_tri = triu(cor_matrix, 1);  % 1 to exclude the diagonal
non_zero_elements = upper_tri(upper_tri ~= 0);   % Get the non-zero elements (actual correlations)

mean_upper_tri = mean(non_zero_elements);
sd_upper_tri = std(non_zero_elements);

figure('Position', [100, 100, 800, 800]);  % Adjust size as needed
imagesc(cor_matrix);
colorbar;
title('Pairwise Correlation Between Subjects');
xlabel('Subject');
ylabel('Subject');
axis square;
colormap(jet);  % Change the colormap if needed
caxis([-1 1]);  % Set color axis for correlation values
subtitle(['Mean Pairwise Correlation: ', num2str(mean_upper_tri, '%.4f'), ', SD: ', num2str(sd_upper_tri, '%.4f')]);

agg_dir = fullfile(studydir, 'ConnMap_PAID_global','Agg');
img_name = fullfile(agg_dir,sprintf('PairwiseCorr_%s.bmp',mask_labels{k}));
saveas(gcf, img_name);

end

close all

%% also correlate subject-wise corr pattern with overall mean 

% load overall conn map
out_dir = fullfile(studydir, 'ConnMap_PAID_global', 'Agg');
mat_name = fullfile(out_dir,'avg_conn_map.mat');
load(mat_name)
    
out_corr = nan(nSubIDlist,n_masks);
for k = 1:n_masks
    tmp_loc = all_idx(:,3)==k;
    tmp_avg_corr = mean(averageMatrix(tmp_loc,:),1);
for subj = 1:nSubIDlist
    out_corr(subj,k) = corr(agg_corr_map(subj,:,k)',tmp_avg_corr');
end
end

figure('Position', [100, 100, 800, 800]);  
for i = 1:4
    subplot(2, 2, i);  % Create a subplot (2 rows, 2 columns, i-th position)
    histogram(out_corr(:, i));  % Plot histogram for the i-th column
    title(['Correlation of subjective FC pattern with group mean:', mask_labels{i}]);
    xlabel('Value');
    ylabel('Frequency');
end

agg_dir = fullfile(studydir, 'ConnMap_PAID_global','Agg');
img_name = fullfile(agg_dir,'Corr_sub_overall.bmp');
saveas(gcf, img_name);


