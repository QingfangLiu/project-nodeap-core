
% compute conn maps across AAL ROIs
% do ROI-based connectivity

% with the sphere ROIs (masks) from TMS experiments

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
% gm, wm, csf all in 2mm space

gm_nii = fullfile(maskpath,'gm_0.1_2mm.nii'); % gray matter mask
gm_dat = spm_read_vols(spm_vol(gm_nii));
gm_idx = find(gm_dat > 0);

wm_nii = fullfile(maskpath,'wm_0.9_2mm.nii'); % white matter mask
wm_dat = spm_read_vols(spm_vol(wm_nii));
wm_idx = find(wm_dat > 0);

csf_nii = fullfile(maskpath,'csf_0.9_2mm.nii'); % CSF mask
csf_dat = spm_read_vols(spm_vol(csf_nii));
csf_idx = find(csf_dat > 0);

%%
% 4 TMS-related sphere ROIs

ROIpath = '/Volumes/X9Pro/NODEAP/ConnectivityMasks';
masks = {'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'
    };  % masks in 2mm space
n_masks = length(masks);
mask_labels = {'aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'};


%% 
% whether or not to include shift
include_shift = true;

%%
for subj = 1 %:nSubIDlist
       
    SubID = SubIDlist(subj).name;
    SubDir = fullfile(MRIdir,SubID);
     
for r = 1 %:n_rest_names 
        
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,SubID)
        
        func_conn_dir = fullfile(studydir, 'FuncConn_AAL_ROIs_PAID', SubID); 
        
        if include_shift
            mat_name = fullfile(func_conn_dir,sprintf('conn_matrix_%s_w_shift.mat',curr_rest));
        else
            mat_name = fullfile(func_conn_dir,sprintf('conn_matrix_%s.mat',curr_rest));
        end
               
        if ~exist(func_conn_dir,'dir')
            mkdir(func_conn_dir)
        end
        
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        
        if exist(mat_name,'file') % if this run has been done
            continue;
        end
        
        path = fullfile (SubDir, 'nifti', 'functional', curr_rest);
        n = dir(fullfile(path, 's6w2*.nii')); % data with 2mm voxel size, smoothed
        
        fmri_data = spm_read_vols(spm_vol(fullfile(path, n.name)));
        
        gm_mean = zeros(nscans,1);
        wm_mean = zeros(nscans,1);
        csf_mean = zeros(nscans,1);
        
        for i = 1:nscans
            tmp = fmri_data(:,:,:,i);
            gm_mean(i) = mean(tmp(gm_idx)); % gray matter: take the mean
            wm_mean(i) = mean(tmp(wm_idx)); % white matter: take the mean
            csf_mean(i) = mean(tmp(csf_idx)); % CSF: take the mean
        end
        
        gm_mean = gm_mean';
              
        mreg = load(fullfile(studydir, 'NRegressor',SubID, sprintf('nuisance_regressors_%s.txt',curr_rest)));
        cleaned_time_series = zeros(nscans, num_rois); % To store cleaned time series
        
        for roi = 1:num_rois
            
            disp(roi)
            roi_mask = (aal_data == roi); % Logical mask for the current ROI
            roi_ts = zeros(nscans, 1); 

            for t = 1:nscans
                roi_ts_t = fmri_data(:,:,:,t);                   % 3D volume for the current time point
                roi_ts(t) = mean(roi_ts_t(roi_mask), 'omitnan'); % Calculate mean ignoring NaNs
            end

            if include_shift   % also filter out gm with temporal shift
                max_shift = 5;  % how much to shift in total
                shifted_matrix = zeros(2 * max_shift + 1, nscans);
                for shift = -max_shift : max_shift
                    if shift < 0
                        shifted_matrix(shift + max_shift + 1, :) = [gm_mean(abs(shift) + 1:end), zeros(1, abs(shift))];
                    elseif shift > 0
                        shifted_matrix(shift + max_shift + 1, :) = [zeros(1, shift), gm_mean(1:end-shift)];
                    else
                        shifted_matrix(shift + max_shift + 1, :) = gm_mean;
                    end
                end
                mreg = [zscore([mreg, shifted_matrix', wm_mean, csf_mean, [1:nscans]']), ones(nscans,1)];
            else
                mreg = [zscore([mreg, gm_mean, wm_mean, csf_mean, [1:nscans]']), ones(nscans,1)]; % add mean(gm), mean(wm), mean(csf), drift, and constant
            end
            
            b = inv(mreg'*mreg)*(mreg'*roi_ts);
            cleaned_time_series(:, roi) = roi_ts - mreg*b;
            
        end

        ROI_time_series = zeros(nscans, n_masks); % To store time series of each ROI
        for i = 1:n_masks
            [targetvol, ~] = spm_read_vols(spm_vol_nifti(fullfile(ROIpath,masks{i})));
            tmploc = (targetvol==1); 
            roi_ts = zeros(nscans, 1); 

            for t = 1:nscans
                roi_ts_t = fmri_data(:,:,:,t); % 3D volume for the current time point
                roi_ts(t) = mean(roi_ts_t(tmploc), 'omitnan'); % Calculate mean ignoring NaNs
            end
            
            if include_shift   % also filter out gm with temporal shift
                max_shift = 5;  % how much to shift in total
                shifted_matrix = zeros(2 * max_shift + 1, nscans);
                for shift = -max_shift : max_shift
                    if shift < 0
                        shifted_matrix(shift + max_shift + 1, :) = [gm_mean(abs(shift) + 1:end), zeros(1, abs(shift))];
                    elseif shift > 0
                        shifted_matrix(shift + max_shift + 1, :) = [zeros(1, shift), gm_mean(1:end-shift)];
                    else
                        shifted_matrix(shift + max_shift + 1, :) = gm_mean;
                    end
                end
                mreg = [zscore([mreg, shifted_matrix', wm_mean, csf_mean, [1:nscans]']), ones(nscans,1)];
            else
                mreg = [zscore([mreg, gm_mean, wm_mean, csf_mean, [1:nscans]']), ones(nscans,1)]; % add mean(gm), mean(wm), mean(csf), drift, and constant
            end
            
            b = inv(mreg'*mreg)*(mreg'*roi_ts);
            ROI_time_series(:, i) = roi_ts - mreg * b;
            
        end
        
        correlation_matrix = zeros(num_rois, n_masks);
        for i = 1:num_rois
            for j = 1:n_masks
                correlation_matrix(i, j) = corr(cleaned_time_series(:, i), ROI_time_series(:, j));
            end
        end

        % save this corr mat 
        save(mat_name,'correlation_matrix')
       
end
end



%% to calculate the mean of the conn maps

for subj = 1:nSubIDlist
        
    SubID = SubIDlist(subj).name;
        
    func_conn_dir = fullfile(studydir, 'FuncConn_AAL_ROIs_PAID', SubID); 
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
        fprintf('%s of %s\n',curr_rest,SubID)
        if count_table(subj,r)<1
            disp('Not enough data, put nan values ...')
            all_corr_map{r,1} = nan(num_rois,n_masks);
            continue;
        end
        mat_name = fullfile(func_conn_dir,sprintf('conn_matrix_%s.mat',curr_rest));
        load(mat_name)
        all_corr_map{r,1} = correlation_matrix;
    end
    
    indices_sham = find(TMS_types == 1);  
    indices_cTBS = find(TMS_types == 2); 
        
    all_maps_cat = cat(3, all_corr_map{:});
    
    mean_corr_map_TMS = nan(num_rois,n_masks,2);
    mean_corr_map_TMS(:,:,1) = nanmean(all_maps_cat(:,:,indices_sham),3);
    mean_corr_map_TMS(:,:,2) = nanmean(all_maps_cat(:,:,indices_cTBS),3);
    sub_corr_map = squeeze(nanmean(all_maps_cat,3));  % mean corr map across sessions 
       
    % plot subject-wise average
    fig = figure('Visible', 'off'); 
    subplot(1,3,1)
    imagesc(sub_corr_map);
    colorbar; 
    caxis([-1, 1]); 
    title('ROI-based Functional Connectivity');
    xlabel('Masks');
    ylabel('ROIs');
    
    subplot(1,3,2)
    imagesc(mean_corr_map_TMS(:,:,1));
    colorbar; 
    caxis([-1, 1]); 
    title('FC (sham)');
    xlabel('Masks');
    ylabel('ROIs');
    
    subplot(1,3,3)
    imagesc(mean_corr_map_TMS(:,:,2));
    colorbar; 
    caxis([-1, 1]); 
    title('FC (cTBS)');
    xlabel('Masks');
    ylabel('ROIs');
    
    if include_shift    
        img_name = fullfile(func_conn_dir,'sub_corr_map_w_shift.bmp');
        mat_name = fullfile(func_conn_dir,'sub_corr_maps_w_shift.mat');
    else
        img_name = fullfile(func_conn_dir,'sub_corr_map.bmp');  
        mat_name = fullfile(func_conn_dir,'sub_corr_maps.mat');   
    end
    
    saveas(gcf, img_name);
    save(mat_name,'mean_corr_map_TMS','sub_corr_map')
    
end




