
%% ========================================================================
%  Script: realign_rest_e1_all_sessions.m
%
%  Purpose
%    - Collect echo-1 (e1) NIfTI volumes across rest sessions for each subject
%    - Run SPM realign (estimate only) on the concatenated session lists
%    - Handle subject/session idiosyncrasies (e.g., NODEAP_44 short run)
%
%  Assumed layout (same as convert script)
%    dat_folder/
%      NODEAP_XX/
%        ├─ nifti/
%        │   └─ <Sess>_me/     % contains <Sess>_rest..._e1.nii (and e2/e3)
%        └─ MRIjobs/
%
%  Notes
%    - Default nscans = 310; NODEAP_44/S1D1 = 205
%    - For NODEAP_41, S3D2 is realigned separately from the other sessions
%    - This script performs Estimate only—it does not create r*.nii resliced images and does not modify NIfTI headers on disk.
%    - SPM writes per-image .mat files storing the estimated rigid-body transforms and an rp_*.txt file 
%    - (6 motion parameters per volume) for each session.
% ========================================================================

%% 0) Environment --------------------------------------------------------------
clear; clc;

dat_folder   = '/Volumes/X9Pro/NODEAP/MRI';

%% 1) Discover subjects ----------------------------------------------------------
SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);       % keep directories only
nSub      = numel(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

% Initialize SPM
spm fmri;
spm_jobman('initcfg');

%% 1) Loop over subjects -------------------------------------------------------
for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;
    SubDir = fullfile(dat_folder,SubID);
    niidir = fullfile(SubDir, 'nifti');
    jobdir = fullfile(SubDir, 'MRIjobs');
    
    fprintf('\n================ %s ================\n', SubID);
    
    % Collect echo-1 filenames across available sessions
    rctr = 0; % count rest sessions (in case of missing session)
    filename = []; % all e1 scans across all sessions

    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
        medir = fullfile(niidir, sprintf('%s_me',curr_rest)); % path to me data

        % Expected #volumes
        if strcmp(SubID,'NODEAP_44') && strcmp(curr_rest,'S1D1')
            nscans = 205;
        else
            nscans = 310;
        end
    
        % Skip if session not present
        if ~exist(medir,'dir')
            fprintf('  ↳ Missing NIfTI dir (skip): %s\n', medir);
            continue;
        end

        rctr = rctr + 1; % to account for potential missing sessions
        
        % organize echo 1 images across sessions
        n_e1 = dir(fullfile(medir, sprintf('%s_rest*_e1.nii',curr_rest)));
        for i = 1:nscans
            filename{rctr}{i,1} = sprintf('%s,%01d', fullfile(medir, n_e1(1).name),i);
        end

    end
    
%% 2) Realign (estimate) setup --------------------------------------------
matlabbatch = [];

if strcmp(SubID,'NODEAP_41') % for NODEAP_41: S3D2 session is realigned separately
    matlabbatch{1}.spm.spatial.realign.estimate.data = filename(1:6);
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 0; % register to first to be faster
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';

    matlabbatch{2}.spm.spatial.realign.estimate.data = filename(7);
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.sep = 4;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.rtm = 0; % register to first to be faster
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.interp = 2;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.weight = '';
else
    matlabbatch{1}.spm.spatial.realign.estimate.data = filename;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 0; % register to first to be faster
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';
end

% save and run job
fname = fullfile(jobdir, 'me_realignment_e1_all_sessions.mat');
save(fname, 'matlabbatch');
spm_jobman('run', matlabbatch);

end

