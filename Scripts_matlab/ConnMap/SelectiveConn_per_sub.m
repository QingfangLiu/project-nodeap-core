
% this analysis selects voxels first based on connectivity map
% then evalulates ROIs' conn with those selected voxels
% to see whether and how the overall conn with those selected voxels changes with TMS

% different from another script
% voxel selection was based on subject-level results

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

n_voxels_vec = 100:50:500;  % an array of voxel number to use


%%
maskpath = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis/Scripts_matlab/GlobalConn/mask');
masks = {'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'
    };  % masks in 3mm space
n_masks = length(masks);

labels = {'aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'};

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

%% organize the conn map across sess and subj
% as a 7*48 big cell array
% this takes a bit time to load all

all_corr_map = cell(n_rest_names,nSubIDlist);
for subj = 1:nSubIDlist
    sn = SubIDlist(subj).name;
    func_conn_dir = fullfile(studydir, 'ConnMap_PAID_global', sn); 
    
    for r = 1:n_rest_names 
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,sn)
        if count_table(subj,r)<1  % only including complete run
            disp('Not enough data, put nan values ...')
            all_corr_map{r,subj} = nan(nvox,n_gm_vox);
            continue;
        end
        mat_name = fullfile(func_conn_dir,sprintf('conn_map_%s.mat',curr_rest));
        load(mat_name)
        all_corr_map{r,subj} = corr_map;
    end
     
end


%% decide the voxels to use
for z = 1:length(n_voxels_vec)
    
    n_top_voxels = n_voxels_vec(z);  
    fprintf('Top %d voxels\n',n_top_voxels)
    
    mean_corr_map_TMS = nan(nSubIDlist,4,2);
    abs_mean_corr_map_TMS = nan(nSubIDlist,4,2);
    
for subj = 1:nSubIDlist
    sn = SubIDlist(subj).name;
        
    func_conn_dir = fullfile(studydir, 'ConnMap_PAID_global', sn); 
    use_corr_map = all_corr_map(:,subj);   % map for this subject
    all_maps_cat = cat(3, use_corr_map{:});
    
    sumMatrix = zeros(nvox, n_gm_vox); 
    countMatrix = zeros(nvox, n_gm_vox);
    for i = 1:numel(use_corr_map)
        tmp = isnan(use_corr_map{i});       % a matrix showing if the entry is nan
        curr_map = use_corr_map{i};
        curr_map(tmp) = 0;                  % make nan entry as zero
        sumMatrix = sumMatrix + curr_map;   % Add the element to the sum
        countMatrix = countMatrix + double(~tmp);
    end
    averageMatrix = sumMatrix ./ countMatrix;
    top_idx = nan(4,n_top_voxels);
    for k = 1:4
        loc = all_idx(:,3)==k;
        conn_vec = mean(averageMatrix(loc,:),1);   % conn vector to this ROI
        [~, sortedIndices] = sort(conn_vec, 'descend');
        top_idx(k,:) = sortedIndices(1:n_top_voxels);   % in gray matter linear index
    end

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
    indices_sham = find(TMS_types == 1);  
    indices_cTBS = find(TMS_types == 2); 
        
    % look at each ROI
    for k = 1:4
        loc = find(all_idx(:,3)==k);
        curr_map = all_maps_cat(loc,top_idx(k,:),:);
        
        mean_corr_map_TMS(subj,k,1) = nanmean(curr_map(:,:,indices_sham),'all');
        mean_corr_map_TMS(subj,k,2) = nanmean(curr_map(:,:,indices_cTBS),'all');
    
        abs_curr_map = abs(curr_map); 
        abs_mean_corr_map_TMS(subj,k,1) = nanmean(abs_curr_map(:,:,indices_sham),'all');
        abs_mean_corr_map_TMS(subj,k,2) = nanmean(abs_curr_map(:,:,indices_cTBS),'all');
    end    
end

    % save the output for this # voxel
    out_dir = fullfile(studydir, 'ConnMap_PAID_global', 'Agg');
    mat_name = fullfile(out_dir,sprintf('FC_corr_sham_cTBS_top%d_per_sub.mat',n_top_voxels));
    save(mat_name,'mean_corr_map_TMS','abs_mean_corr_map_TMS') 

end



    