
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
rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

maskpath = fullfile(HomeDir,'Scripts_matlab/FuncConn_AAL');
aal_path = fullfile(maskpath,'rAAL3v1.nii');
aal_img = spm_vol(aal_path);
aal_data = spm_read_vols(aal_img);
num_rois = max(aal_data(:)); % Total number of ROIs: 170

gm_nii = fullfile(maskpath,'gm_0.1_2mm.nii'); % gray matter mask
gm_dat = spm_read_vols(spm_vol(gm_nii));
gm_idx = find(gm_dat > 0);
n_gm_voxels = length(gm_idx);

include_shift = true; % include temporal shifted gm values

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
for subj = 1:nSubIDlist
       
    SubID = SubIDlist(subj).name;
         
for r = 1:n_rest_names 
        
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,SubID)
        func_conn_dir = fullfile(studydir, 'FuncConn_AAL_ROIs_PAID', SubID); 
        
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        
        if include_shift
            mat_name = fullfile(func_conn_dir,sprintf('conn_matrix_%s_w_shift.mat',curr_rest));
        else
            mat_name = fullfile(func_conn_dir,sprintf('conn_matrix_%s.mat',curr_rest));
        end
        
%         if exist(mat_name,'file') % if this run has been done
%             continue;
%         end
        
        % load data matrix 'dat': nscans * n_gm_voxels
        if include_shift
            tc_name = 'tc_filtered_2mm_w_shift.mat';
        else
            tc_name = 'tc_filtered_2mm.mat';
        end
        
        load(fullfile(studydir, 'FuncConn_AAL_ROIs_PAID', SubID, curr_rest, tc_name));
         
        LOC_roi = cell(num_rois,1); % index of roi within gray matter voxels
        for roi = 1:num_rois
            roi_mask = find(aal_data == roi); 
            [~,LOC_roi{roi,1}] = ismember(roi_mask,gm_idx);
        end
        
        LOC_mask = cell(n_masks,1); % index of mask within gray matter voxels
        for i = 1:n_masks
            [targetvol, ~] = spm_read_vols(spm_vol_nifti(fullfile(ROIpath,masks{i})));
            tmploc = find(targetvol==1); 
            [~,LOC_mask{i,1}] = ismember(tmploc,gm_idx);
        end
        
        correlation_matrix = zeros(num_rois, n_masks);
        for i = 1:num_rois
            for j = 1:n_masks
                use_LOC_roi = LOC_roi{i}(LOC_roi{i}~=0);
                use_roi_dat = mean(dat(:,use_LOC_roi),2);
                
                use_LOC_mask = LOC_mask{j}(LOC_mask{j}~=0);
                use_mask_dat = mean(dat(:,use_LOC_mask),2);
                correlation_matrix(i, j) = corr(use_roi_dat, use_mask_dat);
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
                
        if include_shift
            mat_name = fullfile(func_conn_dir,sprintf('conn_matrix_%s_w_shift.mat',curr_rest));
        else
            mat_name = fullfile(func_conn_dir,sprintf('conn_matrix_%s.mat',curr_rest));
        end
        load(mat_name)
        all_corr_map{r,1} = atanh(correlation_matrix); % fisher z-transformation using atanh
    end
    
    all_maps_cat = cat(3, all_corr_map{:});
    sub_corr_map = squeeze(nanmean(all_maps_cat,3));  % mean corr map across sessions 
    
    indices_sham = find(TMS_types == 1);  
    indices_cTBS = find(TMS_types == 2); 
        
    mean_corr_map_TMS = nan(num_rois,n_masks,2);
    mean_corr_map_TMS(:,:,1) = nanmean(all_maps_cat(:,:,indices_sham),3);
    mean_corr_map_TMS(:,:,2) = nanmean(all_maps_cat(:,:,indices_cTBS),3);
    
    if include_shift
        mat_name = fullfile(func_conn_dir,'sub_corr_maps_w_shift.mat');
    else
        mat_name = fullfile(func_conn_dir,'sub_corr_maps.mat');
    end
    
    save(mat_name,'mean_corr_map_TMS','sub_corr_map')
    
end



