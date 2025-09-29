
%% ========================================================================
%  Script: post_PAID_coreg_norm_smooth_NODEAP_41.m
%
%  Purpose
%    Handle NODEAP_41 specially because S3D2 has different image dims.
%    Steps:
%      (1) Coregister mean EPIs to T1, separately for blocks:
%            • Block 1: D0..S3D1  (mean_fvol_1.nii, sessions 1–6)
%            • Block 2: S3D2      (mean_fvol_2.nii, session 7)
%      (2) Normalize all session 4D EPIs to MNI at 3 mm + smooth 6 mm
%      (3) Normalize again at 2 mm + smooth 6 mm
%
%  Assumes (per subject)
%    <dat_folder>/NODEAP_41/nifti/
%       ├─ anat/D0_T1*.nii
%       ├─ functional/mean_fvol_1.nii     (from D0..S3D1)
%       ├─ functional/mean_fvol_2.nii     (from S3D2 only)
%       └─ functional/<Sess>/fvol_4d.nii  (PAID-combined 4D per session)
%
%  Outputs (under functional/<Sess>/)
%    w3fvol_4d.nii,  s6w3fvol_4d.nii
%    w2fvol_4d.nii,  s6w2fvol_4d.nii
% ========================================================================

clear; clc; 

% ------------------ Common config ------------------
dat_folder = '/Volumes/X9Pro/NODEAP/MRI';

rest_names   = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};

% SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');
global defaults
defaults.stats.maxmem = 2 ^ 34;   % for ~32 GB RAM
defaults.stats.resmem = true;

% ------------------ Locate subject -----------------
SubID = 'NODEAP_41';
SubDir = fullfile(dat_folder, SubID);
niidir = fullfile(SubDir, 'nifti');
jobdir = fullfile(SubDir, 'MRIjobs');
if ~exist(jobdir,'dir'), mkdir(jobdir); end
    
% Inputs
% get anatomical image
aname = dir(fullfile(niidir,'anat','D0_T1*.nii'));
afile = fullfile(niidir,'anat',aname(1).name);
% mean func vol across sessions
func_mean_1 = fullfile(niidir,'functional','mean_fvol_1.nii'); % D0..S3D1
func_mean_2 = fullfile(niidir,'functional','mean_fvol_2.nii'); % S3D2

% If T1 not normalized yet, estimate deformation (y_*.nii) and write 1 mm T1
if isempty(dir(fullfile(niidir,'anat','w*.nii')))
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

% ------------------ Build session lists ------------------
% Block 1 (sessions 1..6): D0..S3D1

func_4d_1 = {};  
for r = 1:6
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubID)
    curr_func_dir = fullfile(niidir,'functional',curr_rest);
    func_4d_1{r,1} = fullfile(curr_func_dir, 'fvol_4d.nii');
end

% Block 2 (session 7): S3D2
func_4d = func_4d_1;
for r = 7
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubID)
    curr_func_dir = fullfile(niidir,'functional',curr_rest);
    func_4d_2 = fullfile(curr_func_dir, 'fvol_4d.nii');
    func_4d{r,1} = func_4d_2;
end

% ------------------ Coregistration (estimate only) ------------------
matlabbatch = [];

if ~isempty(func_4d_1) && exist(func_mean_1,'file')
    % Block 1: mean_fvol_1 → T1, apply to sessions 1..6
    matlabbatch{end+1}.spm.spatial.coreg.estimate.ref    = {afile};
    matlabbatch{end}.spm.spatial.coreg.estimate.source   = {func_mean_1};
    matlabbatch{end}.spm.spatial.coreg.estimate.other    = func_4d_1;
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
end

if ~isempty(func_4d_2) && exist(func_mean_2,'file')
    % Block 2: mean_fvol_2 → T1, apply to session 7
    matlabbatch{end+1}.spm.spatial.coreg.estimate.ref    = {afile};
    matlabbatch{end}.spm.spatial.coreg.estimate.source   = {func_mean_2};
    matlabbatch{end}.spm.spatial.coreg.estimate.other    = func_4d_2;
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
end

% Deformation field from T1 step
def_field = fullfile(niidir,'anat',['y_' anat_list(1).name]);

% ------------------ Normalize @ 3 mm + smooth 6 mm ------------------
matlabbatch{end+1}.spm.spatial.normalise.write.subj.def      = {def_field};
matlabbatch{end}.spm.spatial.normalise.write.subj.resample   = func_4d_all;
matlabbatch{end}.spm.spatial.normalise.write.woptions.bb     = [-78 -112 -70; 78 76 85];
matlabbatch{end}.spm.spatial.normalise.write.woptions.vox    = [3 3 3];
matlabbatch{end}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{end}.spm.spatial.normalise.write.woptions.prefix = 'w3';

matlabbatch{end+1}.spm.spatial.smooth.data(1) = cfg_dep( ...
    'Normalise: Write: Normalised Images (Subj 1)', ...
    substruct('.', 'val', '{}', {numel(matlabbatch)}, '.', 'val', '{}', {1}, '.', 'val', '{}', {1}, '.', 'val', '{}', {1}), ...
    substruct('()', {1}, '.', 'files'));
matlabbatch{end}.spm.spatial.smooth.fwhm   = [6 6 6];
matlabbatch{end}.spm.spatial.smooth.dtype  = 0;
matlabbatch{end}.spm.spatial.smooth.im     = 0;
matlabbatch{end}.spm.spatial.smooth.prefix = 's6';

save(fullfile(jobdir,'Post_PAID_Processing_NODEAP41_3mm.mat'),'matlabbatch');
spm_jobman('run', matlabbatch);
clear matlabbatch

% ------------------ Normalize @ 2 mm + smooth 6 mm ------------------
matlabbatch = [];
matlabbatch{1}.spm.spatial.normalise.write.subj.def      = {def_field};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = func_4d_all;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb   = [-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox  = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w2';

matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep( ...
    'Normalise: Write: Normalised Images (Subj 1)', ...
    substruct('.', 'val', '{}', {1}, '.', 'val', '{}', {1}, '.', 'val', '{}', {1}, '.', 'val', '{}', {1}), ...
    substruct('()', {1}, '.', 'files'));
matlabbatch{2}.spm.spatial.smooth.fwhm   = [6 6 6];
matlabbatch{2}.spm.spatial.smooth.dtype  = 0;
matlabbatch{2}.spm.spatial.smooth.im     = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's6';

save(fullfile(jobdir,'Post_PAID_Processing_NODEAP41_2mm.mat'),'matlabbatch');
spm_jobman('run', matlabbatch);
clear matlabbatch