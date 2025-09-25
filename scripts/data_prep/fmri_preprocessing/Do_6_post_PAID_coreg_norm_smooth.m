
%% ========================================================================
%  Script: post_PAID_coreg_norm_smooth.m
%
%  Purpose
%    Continue rt-fMRI preprocessing on PAID-combined data:
%      (1) Coregister mean EPI to T1, and apply transform to session 4D EPIs
%      (2) Normalize EPIs to MNI space and write out 3 mm & 2 mm versions
%      (3) Smooth both normalized series at 6 mm FWHM
%
%  Assumes (per subject)
%    <dat_folder>/<SubID>/nifti/
%       ├─ anat/D0_T1*.nii                  % anatomical T1 (native)
%       ├─ functional/mean_fvol.nii         % across-session mean of fvol
%       └─ functional/<Sess>/fvol_4d.nii    % PAID-combined 4D per session
%    And a subject×session availability table:
%       /Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx
%
%  Outputs (per subject)
%    Coreg: updates transforms (no new files written)
%    Normalization (3 mm):  w3fvol_4d.nii under functional/<Sess>/
%    Smoothing (on 3 mm):   s6w3fvol_4d.nii    (6 mm)
%    Normalization (2 mm):  w2fvol_4d.nii
%    Smoothing (on 2 mm):   s6w2fvol_4d.nii    (6 mm)
%
%  Notes
%    - Skips NODEAP_41 here (handled separately due to S3D2 dimension mismatch).
%    - T1 normalization first estimates deformation (y_*.nii) and writes 1 mm T1.
%    - Coregister = Estimate only: writes/updates SPM transform files; no new images.
% ========================================================================


%% 0) Setup ----------------------------------------------------------------
clear; clc; tic;

% Data root (fallback allowed)
dat_folder = '/Volumes/X9Pro/NODEAP/MRI';

% Subject list
SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);
nSub      = numel(SubIDlist);

% Availability table (rows=subjects matching SubIDlist order; cols=rest_names)
count_xlsx  = '/Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx';
count_table = readmatrix(count_xlsx);

% Sessions (order matters)
rest_names   = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = numel(rest_names);


spm('Defaults','fMRI')
global defaults
defaults.stats.maxmem = 2 ^ 34; % based on RAM = 32GB
defaults.stats.resmem = true; % allow temp files to be kept in memory

%% 1) Loop over subjects ---------------------------------------------------
for subj = 1:nSubIDlist

SubID = SubIDlist(subj).name;
SubDir = fullfile(dat_folder,SubID);
niidir = fullfile(SubDir, 'nifti');
jobdir = fullfile(SubDir, 'MRIjobs');

if strcmp(SubID, 'NODEAP_41'); continue; end  % handled separately
    fprintf('\n================ %s ================\n', SubID);

% --- Inputs: anatomical + mean EPI
aname = dir(fullfile(niidir,'anat','D0_T1*.nii'));
afile = fullfile(niidir,'anat',aname(1).name);
func_mean = fullfile(niidir,'functional', 'mean_fvol.nii');

if isempty(dir(fullfile(niidir,'anat','w*.nii')))
    % Normalize T1 to generate warping parameters
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {afile};
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {afile};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm =  {fullfile(spm('dir'),'tpm','TPM.nii')};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;

    % Run preprocessing
    job_name = fullfile(jobdir, 'Normalization.mat');
    save(job_name, 'matlabbatch');
    spm_jobman('run', matlabbatch)
    clear matlabbatch
end

% --- Collect available session 4D files for this subject

func_4d = {};    
rctr = 0; % count rest sessions (in case of missing session)

for r = 1:n_rest_names % to loop over rest_names
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubID)
    curr_func_dir = fullfile(niidir,'functional',curr_rest);
    
    % skip if either directory does not exist
    if count_table(subj,r)==0
        disp('Not enough data, skipping this iteration ...')
        continue;
    end
    rctr = rctr + 1; % to account for potential missing sessions
    func_4d{rctr,1} = fullfile(curr_func_dir, 'fvol_4d.nii');
end
    
%% 2) Coregister mean EPI → T1, bring sessions along -------------------
    % Estimate only (no reslicing): updates transforms for func_4d
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {afile};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {func_mean};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = func_4d;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    % Normalise (write) to 3 mm
    def_field = fullfile(niidir,'anat', sprintf('y_%s', aname(1).name));
    matlabbatch{2}.spm.spatial.normalise.write.subj.def = {def_field};
    matlabbatch{2}.spm.spatial.normalise.write.subj.resample = func_4d;
    matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w3';

    % Smooth (6 mm) the 3 mm outputs
    matlabbatch{3}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{3}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{3}.spm.spatial.smooth.dtype = 0;
    matlabbatch{3}.spm.spatial.smooth.im = 0;
    matlabbatch{3}.spm.spatial.smooth.prefix = 's6';

    % Run preprocessing
    job_name = fullfile(jobdir, 'Post_PAID_Processing1.mat');
    save(job_name, 'matlabbatch');
    spm_jobman('run', matlabbatch)
    clear matlabbatch
    
     %% 3) Normalise & smooth again at 2 mm ---------------------------------
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {def_field};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = func_4d;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w2';

    % Smooth EPI
    matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{2}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{2}.spm.spatial.smooth.dtype = 0;
    matlabbatch{2}.spm.spatial.smooth.im = 0;
    matlabbatch{2}.spm.spatial.smooth.prefix = 's6';

    % Run preprocessing
    job_name = fullfile(jobdir, 'Post_PAID_Processing2.mat');
    save(job_name, 'matlabbatch');
    spm_jobman('run', matlabbatch)
    clear matlabbatch

end
    
