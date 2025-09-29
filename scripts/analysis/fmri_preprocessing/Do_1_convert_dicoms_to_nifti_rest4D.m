
%% ========================================================================
%  Script: convert_dicoms_to_nifti_rest4D.m
%  Purpose
%    - Convert DICOMs (T1 + multi-echo resting-state) to NIfTI using dcm2niix
%    - Initialize SPM12 (fMRI)
%    - Collect echo-1 (e1) file lists per session for downstream 4D handling
%
%  Inputs (folder layout assumed)
%    /Volumes/X9Pro/NODEAP/MRI/<SubID>/dicoms/
%      ├─ D0_T1/                 % T1 DICOMs
%      └─ <Sess>_rest/           % Resting-state DICOMs for each session in rest_names
%
%  Outputs (created if missing)
%    /Volumes/X9Pro/NODEAP/MRI/<SubID>/nifti/anat/          % T1 NIfTI
%    /Volumes/X9Pro/NODEAP/MRI/<SubID>/nifti/<Sess>_me/     % Rest NIfTI (multi-echo)
%
%  Requirements
%    - SPM12 on MATLAB path
%    - dcm2niix installed (set absolute path below)
%
%  Notes
%    - Special case: NODEAP_44 at session S1D1 has 205 volumes; others default to 310
%    - Sessions without DICOMs are skipped gracefully
% ========================================================================

%% 0) Environment ----------------------------------------------------------------
clear; clc;

addpath('/Users/liuq13/Applications/spm12/');
spm fmri;  % SPM initialization

dat_folder   = '/Volumes/X9Pro/NODEAP/MRI';
dcm2niipath  = '/Applications/MRIcroGL.app/Contents/Resources/dcm2niix';

% Multi-echo metadata (kept here in case you need them downstream)
TEs    = [14.60, 39.04, 63.48];
nechos = 3;

% Rest sessions to process (order matters for output list)
rest_names    = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names  = numel(rest_names);

% Default number of volumes per resting session
default_nscans = 310;


%% 1) Discover subjects ----------------------------------------------------------
SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);       % keep directories only
nSub      = numel(SubIDlist);

fprintf('Found %d subject folders under %s\n', nSub, dat_folder);

%% 2) Loop over subjects ---------------------------------------------------------
for s = 1:nSub
    SubID = SubIDlist(s).name;
    SubDir = fullfile(dat_folder, SubID);
    niidir = fullfile(SubDir, 'nifti');
    jobdir = fullfile(SubDir, 'MRIjobs');

    if ~exist(jobdir, 'dir'); mkdir(jobdir); end
    if ~exist(niidir, 'dir'); mkdir(niidir); end

    fprintf('\n================ %s ================\n', SubID);

    %% 2a) Convert T1 ------------------------------------------------------------
    anatpath = fullfile(SubDir, 'dicoms', 'D0_T1');
    outdir_anat = fullfile(niidir, 'anat');

    if exist(anatpath, 'dir')
        if ~exist(outdir_anat, 'dir')
            mkdir(outdir_anat);
            run_dcm2niix(dcm2niipath, outdir_anat, anatpath);
        else
            fprintf('T1 already converted: %s\n', outdir_anat);
        end
    else
        fprintf('No T1 DICOMs found for %s at %s (skipping T1)\n', SubID, anatpath);
    end

    %% 2b) Convert resting-state sessions ---------------------------------------
    rctr = 0;                      % counts *converted/available* rest sessions
    filename = {};                 % echo-1 file lists per available session

    for r = 1:n_rest_names
        curr_rest = rest_names{r};
        func_dicoms_path = fullfile(SubDir, 'dicoms', sprintf('%s_rest', curr_rest));
        medir            = fullfile(niidir, sprintf('%s_me', curr_rest));

        fprintf('%s: processing %s ...\n', SubID, curr_rest);

        %  Skip if no DICOMs present
        if ~exist(func_dicoms_path, 'dir') || numel(dir(func_dicoms_path)) <= 2
            fprintf('  ↳ No DICOMs: %s (skipped)\n', func_dicoms_path);
            continue;
        end

        % Determine expected number of scans for this subject/session
        nscans = get_nscans(SubID, curr_rest, default_nscans);

        % Convert if needed
        if ~exist(medir, 'dir')
            mkdir(medir);
            run_dcm2niix(dcm2niipath, medir, func_dicoms_path);
        else
            fprintf('  ↳ NIfTI already exists: %s\n', medir);
        end

        % Build echo-1 filename list for this session
        e1_list = dir(fullfile(medir, sprintf('%s_rest*_e1.nii', curr_rest)));
        if isempty(e1_list)
            fprintf('  ↳ WARNING: no e1 NIfTI found matching %s_rest*_e1.nii\n', curr_rest);
            continue;
        end

        rctr = rctr + 1;              % increment only for sessions that exist
        filename{rctr} = cell(nscans, 1);
        e1_nii = fullfile(medir, e1_list(1).name);

        for i = 1:nscans
            filename{rctr}{i,1} = sprintf('%s,%d', e1_nii, i);
        end

        fprintf('  ↳ Session %s ready (e1 volumes listed: %d)\n', curr_rest, nscans);
    end

    % If you need to persist 'filename', save it here per subject, e.g.:
    % save(fullfile(jobdir, sprintf('%s_e1_lists.mat', SubID)), 'filename');

end

%% ========================== Helper Functions ==================================
function run_dcm2niix(dcm2niipath, outdir, in_dicom_dir)
%RUN_DCM2NIix Thin wrapper to call dcm2niix
    cmd = sprintf('"%s" -o "%s" "%s"', dcm2niipath, outdir, in_dicom_dir);
    fprintf('  ↳ Converting with dcm2niix:\n    %s\n', cmd);
    status = system(cmd);
    if status ~= 0
        warning('  ↳ dcm2niix failed (status=%d) for: %s', status, in_dicom_dir);
    end
end

function n = get_nscans(SubID, sessName, default_n)
%GET_NSCANS Returns expected #volumes per session, with subject-specific overrides.
    if strcmp(SubID, 'NODEAP_44') && strcmp(sessName, 'S1D1')
        n = 205;   % Special-case cutoff
    else
        n = default_n;
    end
end
