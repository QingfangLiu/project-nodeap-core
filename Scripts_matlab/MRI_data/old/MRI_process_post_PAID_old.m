
% This code continues processing rt-fMRI data on PAID data.
% steps include: coregistration, normalization and smoothing
% normalization writes functional EPIs in 3mm 
% note: normalization of anatomical image has been done before

% before everything: the 3D volumes are merged to 4D volume to save some
% space

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
for subj = 1:nSubIDlist

SubID = SubIDlist(subj).name;
SubDir = fullfile(HomeDir,SubID);
niidir = fullfile(SubDir, 'nii');
jobdir = fullfile(SubDir, 'jobs');

% get anatomical image
aname = dir(fullfile(niidir,'anat','D0_T1*.nii'));
afile = fullfile(niidir,'anat',aname(1).name);
    
for r = 1:n_rest_names % to loop over rest_names
    
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubDir)
    
    % get mean EPI as the source file
    curr_me_dir = fullfile(niidir,sprintf('%s_me',curr_rest));
    curr_func_dir = fullfile(SubDir,'nii','functional',curr_rest);
    
    % skip if either directory does not exist
    if ~ (exist(curr_me_dir,'dir') && exist(curr_func_dir,'dir'))
        disp('Either me dir or func dir not existing, skipping this iteration ...')
        fprintf('\n')
        continue;
    end
    
    % skip if both exist
    if ~isempty(dir(fullfile(curr_func_dir,'s6w3f*.nii'))) && ...
       ~isempty(dir(fullfile(curr_func_dir,'s6w2f*.nii'))) && ...
       ~isempty(dir(fullfile(curr_func_dir,'s6f*.nii')))
        disp('Processing has been done, skipping this iteration ...')
        continue;
    end
        
    % unzip needed nii files
    meanfile_name = dir(fullfile(curr_me_dir,'mean*.nii'));
    if isempty(meanfile_name)
        mean_gz = dir(fullfile(curr_me_dir, 'mean*.nii.gz'));
        mean_gz_filename = fullfile(curr_me_dir, mean_gz(1).name);
        gunzip(mean_gz_filename);
        meanfile_name = dir(fullfile(curr_me_dir,'mean*.nii'));
    end
    meanfile = fullfile(curr_me_dir,meanfile_name(1).name);
    
    func_4d = fullfile(curr_func_dir, 'fvol_4d.nii');
    if ~exist(func_4d,'file')
        func_files_list = dir(fullfile(curr_func_dir,'fvol*.nii'));
        func_files = arrayfun(@(x) fullfile(curr_func_dir, x.name), func_files_list, 'UniformOutput', false);
        spm_file_merge(func_files, func_4d, 0);  % Merge 3D volumes into a 4D volume
        for i = 1:numel(func_files)
            delete(func_files{i});
        end
    end
    
if isempty(dir(fullfile(curr_func_dir,'s6w3f*.nii')))
    
    matlabbatch = [];
    % Coregister mean realigned EPI to T1, and bring along the realigned EPIs 
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {afile};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {meanfile};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {func_4d};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    % Normalise write
    fname = fullfile(niidir,'anat', sprintf('y_%s', aname(1).name));
    matlabbatch{2}.spm.spatial.normalise.write.subj.def = {fname};
    matlabbatch{2}.spm.spatial.normalise.write.subj.resample = {func_4d};
    matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w3';

    % Smooth EPI
    matlabbatch{3}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{3}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{3}.spm.spatial.smooth.dtype = 0;
    matlabbatch{3}.spm.spatial.smooth.im = 0;
    matlabbatch{3}.spm.spatial.smooth.prefix = 's6';

    % Run preprocessing
    job_name = fullfile(jobdir, sprintf('Post_PAID_Processing_%s.mat',curr_rest));
    save(job_name, 'matlabbatch');
    spm_jobman('run', matlabbatch)
    clear matlabbatch

end

if isempty(dir(fullfile(curr_func_dir,'s6f*.nii')))
    
    matlabbatch = [];
    % Smooth EPI
    matlabbatch{1}.spm.spatial.smooth.data(1) = {func_4d};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's6';

    % Run preprocessing
    job_name = fullfile(jobdir, sprintf('Post_PAID_Processing2_%s.mat',curr_rest));
    save(job_name, 'matlabbatch');
    spm_jobman('run', matlabbatch)
    clear matlabbatch
    
end

if isempty(dir(fullfile(curr_func_dir,'s6w2f*.nii')))
    
    matlabbatch = [];
    % Coregister mean realigned EPI to T1, and bring along the realigned EPIs 
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {afile};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {meanfile};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {func_4d};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    % Normalise write
    fname = fullfile(niidir,'anat', sprintf('y_%s', aname(1).name));
    matlabbatch{2}.spm.spatial.normalise.write.subj.def = {fname};
    matlabbatch{2}.spm.spatial.normalise.write.subj.resample = {func_4d};
    matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w2';

    % Smooth EPI
    matlabbatch{3}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{3}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{3}.spm.spatial.smooth.dtype = 0;
    matlabbatch{3}.spm.spatial.smooth.im = 0;
    matlabbatch{3}.spm.spatial.smooth.prefix = 's6';

    % Run preprocessing
    job_name = fullfile(jobdir, sprintf('Post_PAID_Processing3_%s.mat',curr_rest));
    save(job_name, 'matlabbatch');
    spm_jobman('run', matlabbatch)
    clear matlabbatch

end

% move the spm ps file to current tedana dir
spm_ps_file = dir(fullfile(pwd,'spm_*.ps'));
if ~isempty(spm_ps_file)
    movefile(fullfile(pwd,spm_ps_file.name),curr_func_dir)
end

end

end

toc

    
