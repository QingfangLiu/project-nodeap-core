% ========================================================================
% Purpose
%   For each available subject/session, compute Pearson correlations
%   between each AAL116 ROI (GM-restricted) and each of 4 OFC spheres
%   (aOFC seed/target, pOFC seed/target) using pre-filtered GM time series.
%
% Inputs (folder structure assumed)
%   dat_folder/<SubID>/nifti/functional/<Sess>/        (source fMRI lives here)
%   dat_folder/FuncConn_AAL/<SubID>/<Sess>/tc_filtered_2mm.mat
%       (contains dat_filtered: [T x N_gm] GM-only time series)
%   project_folder/atlas_masks/tissue_masks/{gm_0.1_2mm.nii, wm_0.9_2mm.nii, csf_0.9_2mm.nii}
%   project_folder/atlas_masks/atlas_aal116/raal.nii         (AAL116 resliced to your grid)
%   project_folder/atlas_masks/ofc_connectivity_masks/*.nii (4 spheres, binary)
%   count_xlsx (rows = subjects in disk order, cols = sessions in rest_names)
%
% Output (per subject/session)
%   dat_folder/FuncConn_AAL/<SubID>/<Sess>/conn_matrix_w_sphere.mat
%       - correlation_matrix : [nROIs x 4] ROI-to-sphere Pearson r
%       - roi_labels         : AAL label integers used (length == nROIs)
%       - sphere_order       : 4 mask filenames (order used)
%
% Key choices
%   - ROIs: derived from AAL labels intersected with GM mask.
%   - Spheres: mapped to GM index-space; mean-TS averaged over GM-intersected voxels.
%   - Simple availability gate: if count_table(subj,r) ~= 1, skip.
%
% Dependencies
%   SPM12 (spm_vol, spm_read_vols), Statistics Toolbox (corr)
% ========================================================================

clear; clc; close all

%% ---------------------- Paths & constants (edit) ------------------------
dat_folder     = '/Volumes/X9Pro/NODEAP/MRI';
study_folder   = '/Volumes/X9Pro/NODEAP';
project_folder = '/Users/liuq13/project-nodeap-core';
count_xlsx     = '/Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx';

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};

maskpath_tissue = fullfile(project_folder, 'atlas_masks', 'tissue_masks');
aal_path        = fullfile(project_folder, 'atlas_masks', 'atlas_aal116', 'raal.nii');
sphere_path     = fullfile(project_folder, 'atlas_masks', 'ofc_connectivity_masks');

sphere_files = { ...
    'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'};
n_spheres = numel(sphere_files);

%% ---------------------- Subjects from disk ------------------------------
SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);
SubIDs    = string({SubIDlist.name})';

%% ---------------------- Availability (simple heuristic) -----------------
count_table_raw = readmatrix(count_xlsx);
% If >=8 columns, assume col1 is ID and sessions are cols 2..8; else 1..7
if size(count_table_raw,2) >= 8
    count_table = count_table_raw(:, 2:8);
else
    count_table = count_table_raw(:, 1:7);
end

%% ---------------------- Load masks & AAL (GM labels) --------------------
Vgm    = spm_vol(fullfile(maskpath_tissue, 'gm_0.1_2mm.nii'));
gm_idx = find(spm_read_vols(Vgm) > 0);

Va      = spm_vol(aal_path);
aal_3D  = spm_read_vols(Va);              % [X Y Z]
aal_gm  = aal_3D(gm_idx);                 % labels only at GM voxels
roi_labels = unique(aal_gm);
roi_labels(roi_labels == 0) = [];         % drop background
nROIs = numel(roi_labels);

%% ---------------------- Map sphere voxels into GM index-space -----------
LOC_sphere = cell(n_spheres,1);
for i = 1:n_spheres
    Vi   = spm_vol(fullfile(sphere_path, sphere_files{i}));
    Yi   = spm_read_vols(Vi);
    vox  = find(Yi == 1);                 % sphere voxels in whole volume
    [tf, loc] = ismember(vox, gm_idx);    % map to columns of dat_filtered
    LOC_sphere{i} = loc(tf);              % GM-intersected sphere indices
end

%% ---------------------- Build task list ---------------------------------
taskList = {};
k = 0;
for subj = 1:numel(SubIDlist)
    for r = 1:numel(rest_names)
        if count_table(subj, r) ~= 1, continue; end
        k = k + 1;
        taskList(k,:) = {SubIDlist(subj).name, rest_names{r}};
    end
end

%% ---------------------- Parallel (optional) -----------------------------
if isempty(gcp('nocreate')), parpool('local'); end

%% ---------------------- Main loop: ROI-to-sphere correlations -----------
parfor t = 1:size(taskList,1)
    SubID     = taskList{t,1};
    curr_rest = taskList{t,2};

    try
        save_dir  = fullfile(study_folder, 'FuncConn_AAL', SubID, curr_rest);
        save_path = fullfile(save_dir, 'conn_matrix_w_sphere.mat');
        if exist(save_path, 'file'), continue; end

        fprintf('ROI-to-sphere FC: %s / %s\n', SubID, curr_rest);

        % Load filtered GM time series
        ts_path = fullfile(save_dir, 'tc_filtered_2mm.mat');
        if ~exist(ts_path, 'file')
            warning('Missing filtered time series: %s', ts_path);
            continue;
        end
        S = load(ts_path);
        if ~isfield(S, 'dat_filtered')
            warning('dat_filtered missing in %s', ts_path);
            continue;
        end
        dat = S.dat_filtered;     % [T x N_gm]

        T = size(dat,1);
        correlation_matrix = nan(nROIs, n_spheres);

        % Precompute ROI logical masks across GM columns
        % (logical over GM voxels: aal_gm == roi_label)
        for iroi = 1:nROIs
            lbl     = roi_labels(iroi);
            vox_roi = (aal_gm == lbl);                   % 1 x N_gm logical
            if ~any(vox_roi), continue; end
            roi_ts = mean(dat(:, vox_roi), 2, 'omitnan');% [T x 1]

            for jsph = 1:n_spheres
                idx_sph = LOC_sphere{jsph};
                if isempty(idx_sph), continue; end
                sph_ts = mean(dat(:, idx_sph), 2, 'omitnan');
                r = corr(roi_ts, sph_ts, 'rows', 'pairwise');
                correlation_matrix(iroi, jsph) = r;
            end
        end

        if ~exist(save_dir, 'dir'), mkdir(save_dir); end
        
        s = struct("correlation_matrix", correlation_matrix, ...
                   "sphere_order", {sphere_files});
        save(save_path, "-fromstruct", s);

    catch ME
        fprintf('Error (ROI-sphere): %s %s -> %s\n', SubID, curr_rest, ME.message);
    end
end