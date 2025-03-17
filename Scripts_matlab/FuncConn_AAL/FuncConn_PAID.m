
% compute conn maps across AAL ROIs
% do ROI-based connectivity

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
for subj = 1:nSubIDlist
       
    SubID = SubIDlist(subj).name;
    SubDir = fullfile(MRIdir,SubID);
     
for r = 1:n_rest_names 
        
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,SubID)
        
        func_conn_dir = fullfile(studydir, 'FuncConn_AAL_PAID', SubID); 
        mat_name = fullfile(func_conn_dir,sprintf('conn_matrix_%s.mat',curr_rest));
        
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
              
        mreg = load(fullfile(studydir, 'NRegressor',SubID, sprintf('nuisance_regressors_%s.txt',curr_rest)));
        cleaned_time_series = zeros(nscans, num_rois); % To store cleaned time series
        
        % Loop through each ROI to extract and clean time series
        for roi = 1:num_rois
            % Create a mask for the current ROI
            roi_mask = (aal_data == roi); % Logical mask for the current ROI
            roi_ts = zeros(nscans, 1); 

            for t = 1:nscans
                roi_ts_t = fmri_data(:,:,:,t); % 3D volume for the current time point
                roi_ts(t) = mean(roi_ts_t(roi_mask), 'omitnan'); % Calculate mean ignoring NaNs
            end

            % Regress out nuisance regressors (from the previously defined nuisance_regressors)
            X = [ones(nscans, 1), mreg, zscore([gm_mean, wm_mean, csf_mean, [1:nscans]'])]; % Design matrix with intercept
            b = (X' * X) \ (X' * roi_ts); % Calculate regression coefficients

            % Compute the residuals (cleaned time series)
            cleaned_time_series(:, roi) = roi_ts - X * b; % Remove the effects of nuisance regressors
        end

        % Compute the correlation matrix
        correlation_matrix = corr(cleaned_time_series);
        
        fig = figure('Visible', 'off'); % plot without showing
        imagesc(correlation_matrix);
        colorbar;
        caxis([-1, 1]); % Adjust these limits based on your data
        title('ROI-based Functional Connectivity Matrix');
        xlabel('ROIs');
        ylabel('ROIs');
        axis square;

        %%
        img_name = fullfile(func_conn_dir,sprintf('conn_matrix_%s.bmp',curr_rest));
        saveas(gcf, img_name);

        % save this corr mat 
        save(mat_name,'correlation_matrix')
       
end
end

%
%% to calculate the mean of the conn maps

for subj = 1:nSubIDlist
        
    SubID = SubIDlist(subj).name;
        
    func_conn_dir = fullfile(studydir, 'FuncConn_AAL_PAID', SubID); 
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
            all_corr_map{r,1} = nan(num_rois,num_rois);
            continue;
        end
        mat_name = fullfile(func_conn_dir,sprintf('conn_matrix_%s.mat',curr_rest));
        load(mat_name)
        all_corr_map{r,1} = correlation_matrix;
    end
    
    indices_sham = find(TMS_types == 1);  
    indices_cTBS = find(TMS_types == 2); 
        
    all_maps_cat = cat(3, all_corr_map{:});
    
    mean_corr_map_TMS = nan(num_rois,num_rois,2);
    mean_corr_map_TMS(:,:,1) = nanmean(all_maps_cat(:,:,indices_sham),3);
    mean_corr_map_TMS(:,:,2) = nanmean(all_maps_cat(:,:,indices_cTBS),3);
    sub_corr_map = squeeze(nanmean(all_maps_cat,3));  % mean corr map across sessions 
       
    % plot subject-wise average
    fig = figure('Visible', 'off'); 
    imagesc(sub_corr_map);
    colorbar; 
    caxis([-1, 1]); 
    title('ROI-based Functional Connectivity Matrix');
    xlabel('ROIs');
    ylabel('ROIs');
    axis square;
    img_name = fullfile(func_conn_dir,'sub_corr_map.bmp');
    saveas(gcf, img_name);
    
    mat_name = fullfile(func_conn_dir,'sub_corr_maps.mat');
    save(mat_name,'mean_corr_map_TMS','sub_corr_map')
    
end


