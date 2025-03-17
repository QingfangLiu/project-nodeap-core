
% copied from SPEEDTMS code
% obtain masks in 3mm space for global connectedness analysis

clear; clc;

%% First: reslice TPM to match current data size & dim

sourceImg = 'TPM.nii'; % copied from /Applications/spm12/tmp/TPM.nii
refImg = '/Volumes/X9Pro/NODEAP/MRI/NODEAP_06_MRI/tedana_output/D0/aic/s6w3desc-optcomDenoised_bold.nii,1'; % Normalized data in 3*3*3 size

% reslice
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.write.ref = {refImg};
matlabbatch{1}.spm.spatial.coreg.write.source = {sourceImg};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4; % use 4th degree spline
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0]; % no wrap
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r3mm';

% run job
spm_jobman('run', matlabbatch)

% a file 'rTPM.nii' should be generated 

%% Second: threshold rTMP to get masks

hdr = spm_vol('r3mmTPM.nii');
nhdr = hdr(1);

% gray matter
x = spm_read_vols(hdr(1)); 
y = zeros(size(x));
y(x > 0.1) = 1; % thereshold at 10%
nhdr.fname = 'gm_0.1_3mm.nii';
spm_write_vol(nhdr,y);

% white matter
x = spm_read_vols(hdr(2)); 
y = zeros(size(x));
y(x > 0.9) = 1; % thereshold at 90%
nhdr.fname = 'wm_0.9_3mm.nii';
spm_write_vol(nhdr,y);

% CSF
x = spm_read_vols(hdr(3)); 
y = zeros(size(x));
y(x > 0.9) = 1; % thereshold at 90%
nhdr.fname = 'csf_0.9_3mm.nii';
spm_write_vol(nhdr,y);


%% Check size

% reference image
y1 = spm_read_vols(spm_vol(refImg));
size(y1)

% gray matter
y2 = spm_read_vols(spm_vol('gm_0.1_3mm.nii'));
size(y2)

% gray matter should have the same size as the reference image






