%% ========================================================================
%  Script: build_subject_mean_fMRI.m
%
%  Purpose
%    - For each subject, concatenate all available session-level combined
%      ME fMRI volumes (fvol_###.nii) into a single 4D per session (if needed),
%      then stack across sessions and write the across-session mean image:
%        nifti/functional/mean_fvol.nii
%
%  Inputs (directory layout)
%    <dat_folder>/<SubID>/nifti/functional/<Sess>/fvol_###.nii
%    count table: /Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx
%      (rows = subjects; cols = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'})
%      nonzero entries indicate the session exists (310 or 205 vols)
%
%  Output (per subject)
%    <dat_folder>/<SubID>/nifti/functional/mean_fvol.nii
%
%  Notes
%    - Skips subject NODEAP_41 due to different image dimension.
%    - If mean_fvol.nii already exists for a subject, that subject is skipped.
%    - If <Sess>/fvol_4d.nii doesn’t exist, merges 3D fvol_###.nii → fvol_4d.nii
%      and deletes the 3D fvol_###.nii afterward (as in your original code).
% ========================================================================

tic
clear; clc;

% Data root 
dat_folder = '/Volumes/X9Pro/NODEAP/MRI';

% Subjects
SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);
nSub      = numel(SubIDlist);

% Read MRI count table
count_xlsx = '/Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx';
count_table = readmatrix(count_xlsx);   % falls back to xlsread-like behavior

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

spm('Defaults','fMRI')
global defaults
defaults.stats.maxmem = 2 ^ 34; % based on RAM = 32GB
defaults.stats.resmem = true; % allow temp files to be kept in memory


%%
for subj = 1:nSubIDlist

SubID = SubIDlist(subj).name;

if strcmp(SubID,'NODEAP_41')
    continue;
end

SubDir = fullfile(HomeDir,SubID);
niidir = fullfile(SubDir, 'nifti');

% if the mean has been created
if exist(fullfile(niidir,'functional', 'mean_fvol.nii'),'file')
    continue;
end

past_volumes = 0;

for r = 1:n_rest_names % to loop over rest_names
    if count_table(subj,r+1)==0
        disp('Not enough data, skipping this iteration ...')
        continue;
    end
    
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubID)
    curr_func_dir = fullfile(niidir,'functional',curr_rest);
    func_4d = fullfile(curr_func_dir, 'fvol_4d.nii');
    
    if ~exist(func_4d,'file')
        disp('Working on 3D to 4D ...')
        func_files_list = dir(fullfile(curr_func_dir,'fvol*.nii'));
        func_files = arrayfun(@(x) fullfile(curr_func_dir, x.name), func_files_list, 'UniformOutput', false);
        spm_file_merge(func_files, func_4d, 0);  % Merge 3D volumes into a 4D volume
        for i = 1:numel(func_files)
            delete(func_files{i});
        end
        func_4d = fullfile(curr_func_dir, 'fvol_4d.nii');
    end
    
    V = spm_vol(func_4d);  % Use SPM to read the NIfTI file
    for i = 1:numel(V)
        mri_data(:,:,:,i + past_volumes) = spm_read_vols(V(i));
    end
    past_volumes = past_volumes + numel(V);

end

    mean_volume_file = V(1);  % Use the header information from the last session    
    mean_volume_file.fname = fullfile(niidir,'functional', 'mean_fvol.nii');
    spm_write_vol(mean_volume_file, mean(mri_data,4));

end

