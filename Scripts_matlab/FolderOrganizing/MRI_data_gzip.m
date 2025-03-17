

% This is to zip nii files to nii.gz to save some disk space
% do this to multi-echo nifti data

clear; clc;
HomeDir = '/Volumes/X9Pro/NODEAP/MRI';

SubIDlist = dir(fullfile(HomeDir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);
count_table = zeros(nSubIDlist,n_rest_names);

for subj = 1:nSubIDlist

SubID = SubIDlist(subj).name;
SubDir = fullfile(HomeDir,SubID);
niidir = fullfile(SubDir, 'nii');

for r = 1:n_rest_names % to loop over rest_names
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubID)
    functdir = fullfile(niidir, sprintf('%s_me',curr_rest)); 
    if exist(functdir,'dir') && numel(dir(fullfile(functdir, '*.nii')))>0
        cd(functdir)
        disp('Working on gzip ...')
        gzip *.nii
        disp('Working on deleting ...')
        delete *.nii
    end
end
end