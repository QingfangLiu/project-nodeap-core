
% this code reslices aal atlas that TK sent to my data space
% this AAL has 116 ROIs
% MNI on 2mm 

clear; clc;

run('../start.m');
sourceImg = './aal116/aal.nii';  

% reslice AAL image to match current data size & dim
refImg = '/Volumes/X9Pro/NODEAP/MRI/NODEAP_06/nifti/functional/D0/w2fvol_4d.nii,1'; % Normalized data in 2*2*2 mm space

% reslice
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.write.ref = {refImg};
matlabbatch{1}.spm.spatial.coreg.write.source = {sourceImg};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0; % use nearest neighbor
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0]; % no wrap
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

% run job
spm_jobman('run', matlabbatch)

% a file 'raal.nii' should be generated 

