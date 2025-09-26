
% ========================================================================
% Purpose
%   Reslice the AAL116 atlas (2 mm MNI) to match a reference fMRI image
%   (e.g., your normalized 2x2x2 mm 4D data). Uses nearest-neighbor
%   interpolation to preserve ROI labels.
%
% Inputs (edit these three):
%   project_folder : repo root containing atlas_masks/atlas_aal116/aal.nii
%   dat_folder     : root of NODEAP MRI derivatives
%   ref_sub, ref_sess, ref_file : a subject/session that already has a
%                                 normalized 2 mm image to serve as the
%                                 spatial reference (first vol used).
%
% Output
%   A resliced atlas named 'raal.nii' saved next to the source atlas.
%
% Requirements
%   SPM12 on path. (This script calls spm('defaults','fmri') and initcfg.)
% ========================================================================

clear; clc;

addpath('/Users/liuq13/Applications/spm12/');

project_folder = '/Users/liuq13/project-nodeap-core';
dat_folder     = '/Volumes/X9Pro/NODEAP/MRI';

ref_sub  = 'NODEAP_06';
ref_sess = 'D0';
ref_file = 'w2fvol_4d.nii';     % normalized 4D image (2 mm); first vol used
%% ------------------------------------------------------------------------

% Paths
atlas_dir   = fullfile(project_folder, 'atlas_masks', 'atlas_aal116');
sourceImg   = fullfile(atlas_dir, 'aal.nii');
refImg_file = fullfile(dat_folder, ref_sub, 'nifti', 'functional', ref_sess, ref_file);
refImg      = [refImg_file ',1'];   % use vol 1 as reference

% Checks
assert(exist(sourceImg, 'file')==2,  'Source atlas not found: %s', sourceImg);
assert(exist(refImg_file,'file')==2, 'Reference image not found: %s', refImg_file);

% Init SPM
spm('defaults','fmri'); spm_jobman('initcfg');

% Batch: Coregister (Write) -> reslice source to reference grid
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.write.ref                 = {refImg};
matlabbatch{1}.spm.spatial.coreg.write.source              = {sourceImg};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp     = 0;        % NN to keep labels
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap       = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask       = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix     = 'r';

% Run
spm_jobman('run', matlabbatch);

% Inform user
fprintf('Done. Resliced atlas written as %s\n', fullfile(atlas_dir,'raal.nii'));


