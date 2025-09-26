% ========================================================================
% Purpose
%   1) Reslice SPM12's TPM.nii (probability maps) to match a reference
%      normalized 2 mm image.
%   2) Threshold the resliced TPM to create GM/WM/CSF binary masks.
%
% Output files (written to tissue mask folder):
%   - r2mmTPM.nii             (TPM resliced to the reference grid)
%   - gm_0.1_2mm.nii          (GM > 0.10)
%   - wm_0.9_2mm.nii          (WM > 0.90)
%   - csf_0.9_2mm.nii         (CSF > 0.90)
%
% Notes
%   - Uses spline interpolation (interp=4) appropriate for probability maps.
%   - Thresholds can be adjusted below (GM_T, WM_T, CSF_T).
%   - Requires SPM12 on MATLAB path.
% ========================================================================

clear; clc;

addpath('/Users/liuq13/Applications/spm12/');

%% ---------------------- Edit these paths if needed ----------------------
project_folder = '/Users/liuq13/project-nodeap-core';
dat_folder     = '/Volumes/X9Pro/NODEAP/MRI';

ref_sub  = 'NODEAP_06';
ref_sess = 'D0';
ref_file = 'w2fvol_4d.nii';                 % normalized 2 mm 4D image

tissue_dir   = fullfile(project_folder, 'atlas_masks', 'tissue_masks');
sourceTPM    = fullfile(tissue_dir, 'TPM.nii'); % copy of SPM TPM.nii placed here
% (If you prefer, set sourceTPM = fullfile(spm('Dir'),'tpm','TPM.nii');)

%% ---------------------- Thresholds --------------------------------------
GM_T  = 0.10;   % GM threshold
WM_T  = 0.90;   % WM threshold
CSF_T = 0.90;   % CSF threshold

%% ---------------------- Setup & checks ----------------------------------
refImg_file = fullfile(dat_folder, ref_sub, 'nifti', 'functional', ref_sess, ref_file);
refImg      = [refImg_file ',1'];  % use first volume as reference

assert(exist(sourceTPM, 'file')==2,  'TPM not found: %s', sourceTPM);
assert(exist(refImg_file,'file')==2, 'Reference image not found: %s', refImg_file);

if ~exist(tissue_dir, 'dir'); mkdir(tissue_dir); end

spm('defaults','fmri'); spm_jobman('initcfg');

%% ---------------------- 1) Reslice TPM to 2 mm grid ---------------------
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.write.ref                 = {refImg};
matlabbatch{1}.spm.spatial.coreg.write.source              = {sourceTPM};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp     = 4;      % 4th-degree spline
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap       = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask       = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix     = 'r2mm';

spm_jobman('run', matlabbatch);

rTPM_path = fullfile(tissue_dir, 'r2mmTPM.nii');
fprintf('Resliced TPM written: %s\n', rTPM_path);

%% ---------------------- 2) Threshold resliced TPM -----------------------
hdr  = spm_vol(rTPM_path);   % 4D: [GM, WM, CSF, ...]
nhdr = hdr(1);               % template header for writing masks

% GM mask
gm  = spm_read_vols(hdr(1));
mask = zeros(size(gm)); mask(gm > GM_T) = 1;
nhdr.fname = fullfile(tissue_dir, sprintf('gm_%g_2mm.nii', GM_T));
spm_write_vol(nhdr, mask);

% WM mask
wm  = spm_read_vols(hdr(2));
mask = zeros(size(wm)); mask(wm > WM_T) = 1;
nhdr.fname = fullfile(tissue_dir, sprintf('wm_%g_2mm.nii', WM_T));
spm_write_vol(nhdr, mask);

% CSF mask
csf = spm_read_vols(hdr(3));
mask = zeros(size(csf)); mask(csf > CSF_T) = 1;
nhdr.fname = fullfile(tissue_dir, sprintf('csf_%g_2mm.nii', CSF_T));
spm_write_vol(nhdr, mask);

fprintf('Masks written to %s\n', tissue_dir);
