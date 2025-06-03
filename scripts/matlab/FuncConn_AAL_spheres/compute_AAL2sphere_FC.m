
% ========================================================================
% Compute ROI-to-sphere functional connectivity maps using preprocessed fMRI data
%
% Description:
% This script calculates Pearson correlation coefficients between each AAL ROI 
% (gray matter-restricted) and each of 4 TMS-related spherical ROIs (in MNI 2mm space)
% across multiple resting-state fMRI sessions and subjects.
%
% Steps:
% 1. Load AAL atlas and gray matter mask.
% 2. Precompute voxel indices for all AAL ROIs and sphere masks within gray matter.
% 3. Iterate through valid subject-session pairs.
% 4. For each, load the filtered time series, extract mean time series per ROI/mask.
% 5. Compute and save correlation matrices (ROI × sphere).
%
% Parallelization:
% Parallelized using `parfor` to speed up subject/session processing.
%
% Outputs:
% correlation_matrix — ROI (116) × sphere (4) matrix for each subject/session
% ========================================================================

clear; clc; close all

studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

maskpath = fullfile(getenv('HOME'),'NODEAP_scripts','Scripts_matlab/FuncConn_AAL/masks');
MRIcount = readtable('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx',"ReadRowNames",true);

gm_idx  = find(spm_read_vols(spm_vol(fullfile(maskpath, 'gm_0.1_2mm.nii')))  > 0);

aal_atlas_3D = spm_read_vols(spm_vol('./aal116/raal.nii'));  % size = [X Y Z]
aal_labels_gm = aal_atlas_3D(gm_idx);  % label of each voxel in gm
roi_labels = unique(aal_labels_gm); % unique labels (should be 0-116)
roi_labels(roi_labels == 0) = [];  % remove background 
nROIs = length(roi_labels); % 116 ROIs

sphere_path = '/Volumes/X9Pro/NODEAP/ConnectivityMasks';
spheres = {'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'
    };  
n_spheres = length(spheres);


%%
% Build job list
taskList = {};
count = 0;
for subj = 1:length(SubIDlist)
    SubID = SubIDlist(subj).name;
    for r = 1:length(rest_names)
        curr_rest = rest_names{r};
        val = MRIcount{SubID, curr_rest};
        if val ~= 1, continue; end  % ignore missing or half-missing sessions
        count = count + 1;
        taskList(count,:) = {SubID, curr_rest};
    end
end

%%
% Start parallel pool if needed
if isempty(gcp('nocreate'))
    parpool('local');
end

%%
% ROI voxel indices
LOC_roi = cell(nROIs,1); 
for roi = 1:nROIs
    roi_voxels = find(aal_labels_gm == roi_labels(roi));
    [~, LOC_roi{roi}] = ismember(roi_voxels, gm_idx);
end

% Sphere mask voxel indices
LOC_sphere = cell(n_spheres,1); 
for i = 1:n_spheres
    targetvol = spm_read_vols(spm_vol_nifti(fullfile(sphere_path,spheres{i})));
    tmploc = find(targetvol==1); 
    [~, LOC_sphere{i}] = ismember(tmploc,gm_idx);
end
   
%%

for t = 1:size(taskList,1)
    SubID = taskList{t,1};
    curr_rest = taskList{t,2};
    
    try
        save_dir = fullfile(studydir, 'FuncConn_AAL', SubID, curr_rest);
        save_path = fullfile(save_dir, 'conn_matrix_w_sphere.mat');
        if exist(save_path,'file'), continue; end
        
        fprintf('Computing aal-sphere correlation for %s of %s\n',curr_rest,SubID)
        
        % load filtered data
        dat_dir = fullfile(studydir, 'FuncConn_AAL', SubID, curr_rest);
        dat_path = fullfile(dat_dir,'tc_filtered_2mm.mat');
        
        if ~exist(dat_path, 'file')
            warning('Missing data file: %s', dat_path);
            continue;
        end
        dat_struct = load(dat_path);
        dat = dat_struct.dat_filtered;
        
        % Compute correlation matrix
        correlation_matrix = zeros(nROIs, n_spheres);
        for i = 1:nROIs
            roi_idx = LOC_roi{i}(LOC_roi{i} ~= 0);
            roi_ts = mean(dat(:, roi_idx), 2);

            for j = 1:n_spheres
                sphere_idx = LOC_sphere{j}(LOC_sphere{j} ~= 0);
                mask_ts = mean(dat(:, sphere_idx), 2);
                correlation_matrix(i, j) = corr(roi_ts, mask_ts);
            end
        end

        if ~exist(save_dir, 'dir'), mkdir(save_dir); end
        s = struct("correlation_matrix", correlation_matrix);
        save(save_path, "-fromstruct", s);  

    catch ME
        fprintf('Error with %s %s: %s\n', SubID, curr_rest, ME.message);
    end
end







