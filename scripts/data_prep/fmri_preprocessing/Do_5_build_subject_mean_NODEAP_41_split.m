
%% ========================================================================
%  Script: build_subject_mean_NODEAP_41_split.m
%
%  Why separate?
%    Subject NODEAP_41 has S3D2 with different voxel dims (e.g., 104×104×78),
%    so we compute means in two blocks:
%      • Block 1: sessions D0..S3D1  ->  functional/mean_fvol_1.nii
%      • Block 2: session  S3D2      ->  functional/mean_fvol_2.nii
%
%  Inputs
%    <dat_folder>/NODEAP_41/nifti/functional/<Sess>/{fvol_###.nii | fvol_4d.nii}
%    Availability table (row order must match dir listing):
%      /Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx
%
%  Outputs
%    <dat_folder>/NODEAP_41/nifti/functional/mean_fvol_1.nii
%    <dat_folder>/NODEAP_41/nifti/functional/mean_fvol_2.nii
% ========================================================================

clear; clc;

% ---------- Config ----------
dat_folder  = '/Volumes/X9Pro/NODEAP/MRI';
rest_names  = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
count_xlsx  = '/Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx';

% SPM setup
spm('Defaults','fMRI');
global defaults
defaults.stats.maxmem = 2 ^ 34;   % for 32 GB RAM
defaults.stats.resmem = true;

% ---------- Locate NODEAP_41 ----------
SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);
names     = {SubIDlist.name};
idx41     = find(strcmp(names, 'NODEAP_41'), 1);
if isempty(idx41)
    error('NODEAP_41 not found under %s', dat_folder);
end
SubID  = 'NODEAP_41';
SubDir = fullfile(dat_folder, SubID);
niidir = fullfile(SubDir, 'nifti');
func_root = fullfile(niidir, 'functional');

% Read availability matrix (rows=subjects, cols=rest_names)
count_table = readmatrix(count_xlsx);

% ---------- Block 1: D0..S3D1 ----------
block1_sessions = 1:6;   % rest_names{1..6}
mri_data = [];
past_volumes = 0;

for r = block1_sessions
    % If count table exists and marks missing, skip
    if size(count_table,1) >= idx41 && size(count_table,2) >= r
        if count_table(idx41, r) == 0
            fprintf('%s — %s: not enough data, skip.\n', SubID, rest_names{r});
            continue;
        }
    end

    curr_rest    = rest_names{r};
    curr_funcdir = fullfile(func_root, curr_rest);
    fprintf('%s — %s\n', SubID, curr_rest);

    % Ensure a 4D fvol exists (merge 3D → 4D if needed)
    func_4d = ensure_4d(curr_funcdir);
    if isempty(func_4d), continue; end

    % Append session to stack
    V = spm_vol(func_4d);
    for i = 1:numel(V)
        mri_data(:,:,:,i + past_volumes) = spm_read_vols(V(i)); 
    end
    past_volumes = past_volumes + numel(V);
end

% Write mean for block 1 
if past_volumes > 0
    mean_hdr       = V(1);  % use last session's header
    if ~exist(func_root,'dir'), mkdir(func_root); end
    mean_hdr.fname = fullfile(func_root, 'mean_fvol_1.nii');
    spm_write_vol(mean_hdr, mean(mri_data,4));
    fprintf('  ↳ Wrote: %s\n', mean_hdr.fname);
else
    fprintf('No data loaded for D0..S3D1 — mean_fvol_1.nii not written.\n');
end

% ---------- Block 2: S3D2 only ----------
r = 7;  % S3D2
if size(count_table,1) >= idx41 && size(count_table,2) >= r && count_table(idx41, r) ~= 0
    curr_rest    = rest_names{r};
    curr_funcdir = fullfile(func_root, curr_rest);
    fprintf('%s — %s (separate dims)\n', SubID, curr_rest);

    func_4d = ensure_4d(curr_funcdir);
    if ~isempty(func_4d)
        V = spm_vol(func_4d);
        for i = 1:numel(V)
            mri_data2(:,:,:,i) = spm_read_vols(V(i)); 
        end
        mean_hdr2       = V(1);
        mean_hdr2.fname = fullfile(func_root, 'mean_fvol_2.nii');
        spm_write_vol(mean_hdr2, mean(mri_data2,4));
        fprintf('  ↳ Wrote: %s\n', mean_hdr2.fname);
    end
else
    fprintf('S3D2 missing or marked absent — mean_fvol_2.nii not written.\n');
end

fprintf('Done in %.1f sec\n', toc);

% ---------- Helper: ensure a 4D fvol exists, return its path -------------
function func_4d = ensure_4d(curr_funcdir)
    func_4d = fullfile(curr_funcdir, 'fvol_4d.nii');
    if exist(func_4d,'file')
        return;
    end
    if ~exist(curr_funcdir,'dir')
        func_4d = '';
        return;
    end
    files3D = dir(fullfile(curr_funcdir, 'fvol_*.nii'));
    if isempty(files3D)
        func_4d = '';
        return;
    end
    fprintf('  ↳ Merging 3D fvol_###.nii → 4D at %s\n', curr_funcdir);
    func_files = arrayfun(@(x) fullfile(curr_funcdir, x.name), files3D, 'UniformOutput', false);
    spm_file_merge(func_files, func_4d, 0);
    % remove 3D inputs (your original behavior)
    for k = 1:numel(func_files), delete(func_files{k}); end
end
