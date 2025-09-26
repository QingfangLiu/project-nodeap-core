
% ========================================================================
% Purpose
%   1) Extract GM-only, nuisance-regressed voxel time series from each
%      available subject/session.
%   2) Average those voxels into AAL116 ROI signals and compute ROI-wise FC.
%
% Inputs (edit the paths in the section below)
%   - dat_folder:
%       <dat_folder>/<SubID>/nifti/functional/<Sess>/s6w2*.nii   (4D norm/smoothed)
%   - project_folder:
%       atlas_masks/tissue_masks/{gm_0.1_2mm.nii, wm_0.9_2mm.nii, csf_0.9_2mm.nii}
%       atlas_masks/atlas_aal116/raal.nii
%   - count_xlsx: availability matrix (rows = SubIDs sorted by disk order,
%       cols = sessions in rest_names). Value==1 means a valid session.
%
% Outputs (per subject/session)
%   <study_folder>/FuncConn_AAL/<SubID>/<Sess>/tc_filtered_2mm.mat   (dat_filtered)
%   <study_folder>/FuncConn_AAL/<SubID>/<Sess>/conn_matrix.mat       (FC_matrix)
%
% Key choices
%   - GM/WM/CSF masks: thresholds already baked into filenames (0.1, 0.9, 0.9).
%   - Nuisance regressors: loads .../NRegressor/<SubID>/nuisance_regressors_<Sess>.txt,
%     then augments with GM/WM/CSF means, linear drift, z-scores, plus intercept.
%   - Availability: simple check `if count_table(subj,r) ~= 1, continue; end`.
%
% Dependencies
%   - SPM12 on MATLAB path (spm_vol, spm_read_vols)
%   - Statistics Toolbox (zscore)
% ========================================================================

clc; clear;

%% ------------------- Paths & constants (edit these) ---------------------
dat_folder     = '/Volumes/X9Pro/NODEAP/MRI';
study_folder   = '/Volumes/X9Pro/NODEAP';
project_folder = '/Users/liuq13/project-nodeap-core';
count_xlsx     = '/Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx';

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
N_EXPECTED_SCANS = 310;       % set to [] to skip the scan-count check

maskpath = fullfile(project_folder, 'atlas_masks', 'tissue_masks');
aal_path = fullfile(project_folder, 'atlas_masks', 'atlas_aal116', 'raal.nii');

%% ------------------- Discover subjects by disk order --------------------
SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);
SubIDs    = string({SubIDlist.name})';

%% ------------------- Load availability (simple heuristic) ---------------
count_table_raw = readmatrix(count_xlsx);
% If >=8 columns, assume col1 is ID and sessions are cols 2..8, else 1..7
if size(count_table_raw,2) >= 8
    count_table = count_table_raw(:, 2:8);
else
    count_table = count_table_raw(:, 1:7);
end

%% ------------------- Load masks & AAL atlas -----------------------------
Vgm = spm_vol(fullfile(maskpath, 'gm_0.1_2mm.nii'));
Vwm = spm_vol(fullfile(maskpath, 'wm_0.9_2mm.nii'));
Vcs = spm_vol(fullfile(maskpath, 'csf_0.9_2mm.nii'));
gm_idx  = find(spm_read_vols(Vgm) > 0);
wm_idx  = find(spm_read_vols(Vwm) > 0);
csf_idx = find(spm_read_vols(Vcs) > 0);

% AAL (resliced) atlas and GM-restricted labels
Va      = spm_vol(aal_path);
aal_3D  = spm_read_vols(Va);                 % [X Y Z]
aal_gm  = aal_3D(gm_idx);                    % labels for GM voxels only
roi_labels = unique(aal_gm);
roi_labels(roi_labels == 0) = [];            % drop background
nROIs = numel(roi_labels);

%% ------------------- Build task list from availability ------------------
taskList = {}; k = 0;
for subj = 1:numel(SubIDlist)
    SubID = SubIDlist(subj).name;
    for r = 1:numel(rest_names)
        if count_table(subj, r) ~= 1, continue; end
        k = k + 1;
        taskList{k,1} = SubID;
        taskList{k,2} = rest_names{r};
    end
end

%% ------------------- Start parallel pool (if not already) --------------
if isempty(gcp('nocreate')), parpool('local'); end
%%

%% ------------------- Part 1: Filtered GM voxel time series --------------
parfor t = 1:size(taskList,1)
    SubID     = taskList{t,1};
    curr_rest = taskList{t,2};

    try
        subj_dir   = fullfile(dat_folder, SubID);
        func_dir   = fullfile(subj_dir, 'nifti', 'functional', curr_rest);
        
        save_dir   = fullfile(study_folder, 'FuncConn_AAL', SubID, curr_rest);
        save_path  = fullfile(save_dir, 'tc_filtered_2mm.mat');

        if exist(save_path, 'file'), continue; end

        vols = dir(fullfile(func_dir, 's6w2*.nii'));
        if isempty(vols), continue; end

        fprintf('Extracting filtered GM time series: %s / %s\n', SubID, curr_rest);

        % Read 4D normalized/smoothed data
        Vimgs     = spm_vol(fullfile(func_dir, vols(1).name));
        data_4D   = spm_read_vols(Vimgs);         % [X Y Z T]
        nscans    = size(data_4D, 4);

        if ~isempty(N_EXPECTED_SCANS) && nscans ~= N_EXPECTED_SCANS
            warning('Scan count mismatch %s %s: expected %d, found %d', ...
                    SubID, curr_rest, N_EXPECTED_SCANS, nscans);
            continue;
        end

        % Reshape to [T x Vox] and restrict to GM
        dat     = reshape(data_4D, [], nscans)';  % [T x Vox]
        dat_gm  = dat(:, gm_idx);

        % Tissue means & drift
        wm_mean = mean(dat(:, wm_idx), 2);
        csf_mean= mean(dat(:, csf_idx), 2);
        gm_mean = mean(dat_gm, 2);
        drift   = (1:nscans)';

        % Nuisance regressors (external)
        nreg_dir  = fullfile(study_folder, 'NRegressor', SubID);
        nreg_file = fullfile(nreg_dir, sprintf('nuisance_regressors_%s.txt', curr_rest));
        if ~exist(nreg_file, 'file')
            warning('Missing nuisance file: %s', nreg_file);
            continue;
        end
        mreg = load(nreg_file);                   % [T x R]
        mreg = [mreg, gm_mean, wm_mean, csf_mean, drift];
        mreg = zscore(mreg);
        mreg = [mreg, ones(nscans,1)];            % intercept

        % Regress out nuisance from GM voxels
        B            = mreg \ dat_gm;
        dat_filtered = dat_gm - mreg * B;         % [T x Ngm]

        % Save
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end
        SOUT.dat_filtered = dat_filtered;
        save(save_path, '-struct', 'SOUT');

    catch ME
        fprintf('Error (filter): %s %s -> %s\n', SubID, curr_rest, ME.message);
    end
end


%% ------------------- Part 2: ROI signals & FC matrices ------------------
parfor t = 1:size(taskList,1)
    SubID     = taskList{t,1};
    curr_rest = taskList{t,2};

    try
        save_dir  = fullfile(study_folder, 'FuncConn_AAL', SubID, curr_rest);
        ts_path   = fullfile(save_dir, 'tc_filtered_2mm.mat');
        fc_path   = fullfile(save_dir, 'conn_matrix.mat');
        if exist(fc_path, 'file'), continue; end
        if ~exist(ts_path, 'file'), continue; end

        fprintf('Computing AAL FC: %s / %s\n', SubID, curr_rest);

        S            = load(ts_path);                 % dat_filtered [T x Ngm]
        dat_filtered = S.dat_filtered;
        nscans       = size(dat_filtered, 1);

        % Build AAL ROI signals by averaging GM voxels with each label
        ROI_signals = zeros(nscans, nROIs);
        for i = 1:nROIs
            lbl = roi_labels(i);
            vox = (aal_gm == lbl);                   % logical over GM voxels
            if any(vox)
                ROI_signals(:, i) = mean(dat_filtered(:, vox), 2, 'omitnan');
            else
                ROI_signals(:, i) = NaN;
            end
        end

        % ROI x ROI Pearson correlation (pairwise rows to handle any NaNs)
        FC_matrix = corr(ROI_signals, 'rows', 'pairwise');

        SOUT.FC_matrix = FC_matrix;
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end
        save(fc_path, '-struct', 'SOUT');

    catch ME
        fprintf('Error (FC): %s %s -> %s\n', SubID, curr_rest, ME.message);
    end
end


