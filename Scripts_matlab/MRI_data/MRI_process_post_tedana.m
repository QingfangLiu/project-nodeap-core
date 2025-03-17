
% This code continues processing rt-fMRI data on optimally combined and denoised data
% from tedana output 'desc-optcomDenoised_bold.nii.gz'
% steps include: coregistration, normalization and smoothing
% normalization writes functional EPIs in 2mm and 3mm 
% smoothing was done with fwhm 6mm
% later added one with only spatial smoothing
% output includes 's6w2*.nii' 's6w3*.nii' 's6desc*.nii'

tic
clear; clc;
HomeDir = '/Volumes/X9Pro/NODEAP/MRI';

nscans = 310; % 310 volumes for ME sequnce
SubIDlist = dir(fullfile(HomeDir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);
%spm fmri
spm_figure('CreateWin','Graphics'); % print output without the entire GUI

%%
for subj = 1:nSubIDlist

SubID = SubIDlist(subj).name;
SubDir = fullfile(HomeDir,SubID);
niidir = fullfile(SubDir, 'nii');
jobdir = fullfile(SubDir, 'jobs');

if ~exist(jobdir,'dir') 
   mkdir(jobdir)
end

% get anatomical image
aname = dir(fullfile(niidir,'anat','D0_T1*.nii'));
afile = fullfile(niidir,'anat',aname(1).name);
    
for r = 1:n_rest_names % to loop over rest_names
    
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubDir)
    
    % get mean EPI as the source file
    curr_me_dir = fullfile(niidir,sprintf('%s_me',curr_rest));
    curr_tedana_dir = fullfile(SubDir,'tedana_output',curr_rest,'aic');
    
    % skip if either directory does not exist
    if ~ (exist(curr_me_dir,'dir') && exist(curr_tedana_dir,'dir'))
        disp('Either me dir or tedana dir not existing, skipping this iteration ...')
        continue;
    end
    
    % skip if all three types of processed files exist
    if ~isempty(dir(fullfile(curr_tedana_dir,'s6w2*.nii'))) && ...
        ~isempty(dir(fullfile(curr_tedana_dir,'s6w3*.nii'))) && ...
        ~isempty(dir(fullfile(curr_tedana_dir,'s6desc*.nii')))
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
    
    func_files = fullfile(curr_tedana_dir,'desc-optcomDenoised_bold.nii');
    if ~exist(func_files,'file')
        func_files_gz = fullfile(curr_tedana_dir,'desc-optcomDenoised_bold.nii.gz');
        gunzip(func_files_gz)
    end
    
if isempty(dir(fullfile(curr_tedana_dir,'s6w2*.nii')))
    
    matlabbatch = [];
    % Coregister mean realigned EPI to T1, and bring along the realigned EPIs 
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {afile};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {meanfile};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {func_files};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    % Normalize T1 to generate warping parameters
    matlabbatch{2}.spm.spatial.normalise.estwrite.subj.vol = {afile};
    matlabbatch{2}.spm.spatial.normalise.estwrite.subj.resample = {afile};
    matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.tpm =  {fullfile(spm('dir'),'tpm','TPM.nii')};
    matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{2}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{2}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1];
    matlabbatch{2}.spm.spatial.normalise.estwrite.woptions.interp = 4;

    % Normalise write
    fname = fullfile(niidir,'anat', sprintf('y_%s', aname(1).name));
    matlabbatch{3}.spm.spatial.normalise.write.subj.def = {fname};
    matlabbatch{3}.spm.spatial.normalise.write.subj.resample = {func_files};
    matlabbatch{3}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{3}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{3}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{3}.spm.spatial.normalise.write.woptions.prefix = 'w2';

    % Smooth EPI
    matlabbatch{4}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{4}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{4}.spm.spatial.smooth.dtype = 0;
    matlabbatch{4}.spm.spatial.smooth.im = 0;
    matlabbatch{4}.spm.spatial.smooth.prefix = 's6';

    % Run preprocessing
    job_name = fullfile(jobdir, sprintf('Post_tedana_Processing_%s.mat',curr_rest));
    save(job_name, 'matlabbatch');
    spm_jobman('run', matlabbatch)
    clear matlabbatch

end

if isempty(dir(fullfile(curr_tedana_dir,'s6w3*.nii')))
    matlabbatch = [];
    % Normalise write
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {fname};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {func_files};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w3';

    % Smooth EPI
    matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{2}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{2}.spm.spatial.smooth.dtype = 0;
    matlabbatch{2}.spm.spatial.smooth.im = 0;
    matlabbatch{2}.spm.spatial.smooth.prefix = 's6';

    % Run preprocessing
    job_name = fullfile(jobdir, sprintf('Post_tedana_Processing2_%s.mat',curr_rest));
    save(job_name, 'matlabbatch');
    spm_jobman('run', matlabbatch)
    clear matlabbatch
end

% no coregistration or normalization, only spatial smoothing
if isempty(dir(fullfile(curr_tedana_dir,'s6desc*.nii')))
    
    matlabbatch = [];
    % Smooth EPI
    matlabbatch{1}.spm.spatial.smooth.data(1) = {func_files};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's6';

    % Run preprocessing
    job_name = fullfile(jobdir, sprintf('Post_tedana_Processing3_%s.mat',curr_rest));
    save(job_name, 'matlabbatch');
    spm_jobman('run', matlabbatch)
    clear matlabbatch

end

% move the spm ps file to current tedana dir
spm_ps_file = dir(fullfile(pwd,'spm_*.ps'));
if ~isempty(spm_ps_file)
    movefile(fullfile(pwd,spm_ps_file.name),curr_tedana_dir)
end

end

end

toc

    
