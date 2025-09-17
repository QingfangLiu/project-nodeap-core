
% This code continues processing rt-fMRI data on PAID data.
% steps include: coregistration, normalization and smoothing
% normalization writes functional EPIs in 3mm 
% note: normalization of anatomical image has been done before

% before everything: the 3D volumes are merged to 4D volume to save some
% space

% redo this for NODEAP_41 S3D2 only
% this session went through another coregister-reslice step


tic
clear; clc;
HomeDir = '/Volumes/X9Pro/NODEAP/MRI';

nscans = 310; % 310 volumes for ME sequnce
SubIDlist = dir(fullfile(HomeDir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);
spm_figure('CreateWin','Graphics'); % print output without the entire GUI

spm('Defaults','fMRI')
global defaults
defaults.stats.maxmem = 2 ^ 34; % based on RAM = 32GB
defaults.stats.resmem = true; % allow temp files to be kept in memory

%%
% for subj = 1:nSubIDlist

% SubID = SubIDlist(subj).name;

SubID = 'NODEAP_41';

SubDir = fullfile(HomeDir,SubID);
niidir = fullfile(SubDir, 'nii');
jobdir = fullfile(SubDir, 'jobs');

% get anatomical image
aname = dir(fullfile(niidir,'anat','D0_T1*.nii'));
afile = fullfile(niidir,'anat',aname(1).name);
    
% for r = 1:n_rest_names % to loop over rest_names
    r = 7;
    
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubDir)
    
    % get mean EPI as the source file
    curr_me_dir = fullfile(niidir,sprintf('%s_me',curr_rest));
    curr_func_dir = fullfile(SubDir,'nii','functional',curr_rest);
    
    % unzip needed nii files
    meanfile_name = dir(fullfile(curr_me_dir,'mean*.nii'));
    if isempty(meanfile_name)
        mean_gz = dir(fullfile(curr_me_dir, 'mean*.nii.gz'));
        mean_gz_filename = fullfile(curr_me_dir, mean_gz(1).name);
        gunzip(mean_gz_filename);
        meanfile_name = dir(fullfile(curr_me_dir,'mean*.nii'));
    end
    meanfile = fullfile(curr_me_dir,meanfile_name(1).name);
    func_4d = fullfile(curr_func_dir, 'rfvol_4d.nii');
    
if isempty(dir(fullfile(curr_func_dir,'s6rf*.nii')))
    
    matlabbatch = [];
    % Smooth EPI
    matlabbatch{1}.spm.spatial.smooth.data(1) = {func_4d};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's6';

    % Run preprocessing
    spm_jobman('run', matlabbatch)
    clear matlabbatch
    
end

% end

% end

toc

    
