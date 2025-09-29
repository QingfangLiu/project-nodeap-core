
% ========================================================================
%
% Purpose:
%   Compute resting-state functional connectivity (FC) between four
%   right-hemisphere OFC spheres (aOFC seed/target, pOFC seed/target)
%   across all subjects and sessions that have available data.
%
% Inputs:
%   - dat_folder: folder containing subject subfolders 'NODEAP*'
%                 with derivatives:
%                 FuncConn_AAL/<SubID>/<Sess>/tc_filtered_2mm.mat
%                 (variable: dat_filtered [T x Nvox_in_GM])
%   - project_folder: repo folder containing masks under:
%                 atlas_masks/tissue_masks/gm_0.1_2mm.nii
%                 atlas_masks/ofc_connectivity_masks/*.nii
%   - count_xlsx: availability matrix (rows = subjects, cols = sessions)
%                 value == 1 means session available
%
% Outputs:
%   For each valid subject-session pair:
%     FuncConn_AAL/<SubID>/<Sess>/conn_matrix_4spheres.mat
%       - correlation_matrix : 4x4 Pearson correlation matrix
%       - sphere_order       : names of sphere masks in order
%
% Dependencies:
%   SPM12 (spm_vol, spm_read_vols, corr)
% ========================================================================

clear; clc; close all

%% ---------------------- Paths & constants -------------------------------
study_folder   = '/Volumes/X9Pro/NODEAP';
dat_folder     = '/Volumes/X9Pro/NODEAP/MRI';
project_folder = '/Users/liuq13/project-nodeap-core';
count_xlsx     = '/Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx';

maskpath_tissue = fullfile(project_folder, 'atlas_masks', 'tissue_masks');
maskpath_sphere = fullfile(project_folder, 'atlas_masks', 'ofc_connectivity_masks');

% Subject list (disk order)
SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);
SubIDs    = string({SubIDlist.name})';
nSub      = numel(SubIDs);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};

maskpath = fullfile(project_folder,'atlas_masks','tissue_masks');

% Availability table (rows must align with SubID order above or contain ID col)
count_table_raw = readmatrix(count_xlsx);
% Heuristic: if table has >=8 cols, assume col 1 is an ID and sessions are cols 2..8
if size(count_table_raw,2) >= 8
    count_table = count_table_raw(:, 2:8);
else
    count_table = count_table_raw(:, 1:7);
end

%% ---------------------- Load GM mask ------------------------------------
Vgm   = spm_vol(fullfile(maskpath_tissue, 'gm_0.1_2mm.nii'));
Ygm   = spm_read_vols(Vgm);
gm_idx = find(Ygm > 0);

sphere_files = { ...
    'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'};
n_spheres = numel(sphere_files);

%% ---------------------- Sphere voxel indices (within GM) ----------------
LOC_sphere = cell(n_spheres,1);
for i = 1:n_spheres
    Vi   = spm_vol(fullfile(maskpath_sphere, sphere_files{i}));
    Yi   = spm_read_vols(Vi);
    vox  = find(Yi == 1);
    [~, loc] = ismember(vox, gm_idx);
    LOC_sphere{i} = loc(loc > 0); % keep only voxels inside GM
end

%% ---------------------- Build subject-session task list -----------------
taskList = {};
count = 0;
for subj = 1:nSub
    SubID = SubIDlist(subj).name;
    for r = 1:numel(rest_names)
        if count_table(subj,r) ~= 1, continue; end  % ignore missing sessions
        count = count + 1;
        taskList(count,:) = {SubID, rest_names{r}};
    end
end

%% ---------------------- Main loop: FC computation -----------------------
for t = 1:size(taskList,1)
    SubID     = taskList{t,1};
    curr_rest = taskList{t,2};

    try
        
        fprintf('Computing sphere-sphere FC for %s of %s\n', curr_rest, SubID);

        dat_dir  = fullfile(study_folder, 'FuncConn_AAL', SubID, curr_rest);
        dat_path = fullfile(dat_dir, 'tc_filtered_2mm.mat');
        save_path = fullfile(dat_dir, 'conn_matrix_4spheres.mat');
        if exist(save_path, 'file'), continue; end
        if ~exist(dat_path, 'file')
            warning('Missing data file: %s', dat_path);
            continue;
        end

        S = load(dat_path);  % expects dat_filtered [T x Ngm]
        if ~isfield(S, 'dat_filtered')
            warning('dat_filtered variable not found in %s', dat_path);
            continue;
        end
        dat = S.dat_filtered;

        % Extract mean time series for each sphere
        Tlen   = size(dat,1);
        ts_all = zeros(Tlen, n_spheres);
        for i = 1:n_spheres
            idx = LOC_sphere{i};
            if isempty(idx)
                warning('Sphere %s has no voxels in GM.', sphere_files{i});
                ts_all(:,i) = NaN;
            else
                ts_all(:,i) = mean(dat(:, idx), 2, 'omitnan');
            end
        end

        % Compute correlation matrix
        correlation_matrix = corr(ts_all, 'rows', 'pairwise');

        SOUT.correlation_matrix = correlation_matrix;
        SOUT.sphere_order       = sphere_files;
        save(save_path, '-struct', 'SOUT');

    catch ME
        fprintf('Error with %s %s: %s\n', SubID, curr_rest, ME.message);
    end
end







