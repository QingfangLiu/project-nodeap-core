
% this code reslices aal atlas that TK sent to my data space
% MNI on 2mm 

clear; clc;
sourceImg = './TKsent/aal.nii';  

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

% check file dimension match
% reference image
y1 = spm_read_vols(spm_vol(refImg));
size(y1)

newfile = './TKsent/raal.nii'; 
y2 = spm_read_vols(spm_vol(newfile));
size(y2)
