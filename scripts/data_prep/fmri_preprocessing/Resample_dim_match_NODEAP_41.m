

% reslice images from 104*104*78 dimension to 104*104*60
% for NODEAP_41, S3D2 data

functdir = '/Volumes/X9Pro/NODEAP/MRI/NODEAP_41/nii/functional/S3D2';
fnames = dir(fullfile(functdir,'f*.nii')); 
ffile = fullfile(functdir,fnames(1).name);

ffiles = spm_select('expand', ffile); % expand the volumes from a 4D nifti file
ffiles = cellstr(ffiles); 

clear matlabbatch
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {'/Volumes/X9Pro/NODEAP/MRI/NODEAP_41/nii/functional/S3D1/fvol_4d.nii,1'}; % dim: 104*104*60
matlabbatch{1}.spm.spatial.coreg.estwrite.source = ffiles(1); % dim: 104*104*78
matlabbatch{1}.spm.spatial.coreg.estwrite.other = ffiles;
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

spm_jobman('run', matlabbatch);
