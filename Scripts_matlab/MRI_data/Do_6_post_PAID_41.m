
% This code continues processing rt-fMRI data on PAID data.
% steps include: coregistration, normalization and smoothing
% normalization writes functional EPIs in 3mm 

clear; clc;
HomeDir = '/Volumes/X9Pro/NODEAP/MRI';
[~, scriptName, ~] = fileparts(mfilename('fullpath'));

SubIDlist = dir(fullfile(HomeDir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');
userHome = getenv('HOME');
if strcmp(userHome,'/Users/liuq13')
    ghostscriptPath = '/opt/homebrew/bin/gs'; 
else
    ghostscriptPath = '/usr/local/bin/gs'; 
end

if ~exist(ghostscriptPath)
    error('The ghostscript folder does not exist: %s', ghostscriptPath);
end

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);
% spm fmri

spm('Defaults','fMRI')
global defaults
defaults.stats.maxmem = 2 ^ 34; % based on RAM = 32GB
defaults.stats.resmem = true; % allow temp files to be kept in memory

%%
for subj = 18

    tic
SubID = SubIDlist(subj).name;
SubDir = fullfile(HomeDir,SubID);
niidir = fullfile(SubDir, 'nifti');
jobdir = fullfile(SubDir, 'MRIjobs');

output_file = sprintf('output_%s_%s.txt', scriptName, string(datetime));
diary(output_file); % keep notes of outputs

% get anatomical image
aname = dir(fullfile(niidir,'anat','D0_T1*.nii'));
afile = fullfile(niidir,'anat',aname(1).name);
% mean func vol across sessions
func_mean_1 = fullfile(niidir,'functional', 'mean_fvol_1.nii');
func_mean_2 = fullfile(niidir,'functional', 'mean_fvol_2.nii');

if isempty(dir(fullfile(niidir,'anat','w*.nii')))
    % Normalize T1 to generate warping parameters
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {afile};
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {afile};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm =  {fullfile(spm('dir'),'tpm','TPM.nii')};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;

    % Run preprocessing
    job_name = fullfile(jobdir, 'Normalization.mat');
    save(job_name, 'matlabbatch');
    spm_jobman('run', matlabbatch)
    clear matlabbatch
end

func_4d_1 = {};  

for r = 1:6
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubID)
    curr_func_dir = fullfile(niidir,'functional',curr_rest);
    func_4d_1{r,1} = fullfile(curr_func_dir, 'fvol_4d.nii');
end

func_4d = func_4d_1;

for r = 7
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubID)
    curr_func_dir = fullfile(niidir,'functional',curr_rest);
    func_4d_2 = fullfile(curr_func_dir, 'fvol_4d.nii');
    func_4d{r,1} = func_4d_2;
end

%     matlabbatch = [];
%     % Coregister mean realigned EPI to T1, and bring along the realigned EPIs 
%     matlabbatch{1}.spm.spatial.coreg.estimate.ref = {afile};
%     matlabbatch{1}.spm.spatial.coreg.estimate.source = {func_mean_1};
%     matlabbatch{1}.spm.spatial.coreg.estimate.other = func_4d_1;
%     matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
%     matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
%     matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
%     matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
%     
%     % Coregister mean realigned EPI to T1, and bring along the realigned EPIs 
%     matlabbatch{2}.spm.spatial.coreg.estimate.ref = {afile};
%     matlabbatch{2}.spm.spatial.coreg.estimate.source = {func_mean_2};
%     matlabbatch{2}.spm.spatial.coreg.estimate.other = {func_4d_2};
%     matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
%     matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
%     matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
%     matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
% 
%     % Normalise write
%     fname = fullfile(niidir,'anat', sprintf('y_%s', aname(1).name));
%     matlabbatch{3}.spm.spatial.normalise.write.subj.def = {fname};
%     matlabbatch{3}.spm.spatial.normalise.write.subj.resample = func_4d;
%     matlabbatch{3}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
%     matlabbatch{3}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
%     matlabbatch{3}.spm.spatial.normalise.write.woptions.interp = 4;
%     matlabbatch{3}.spm.spatial.normalise.write.woptions.prefix = 'w3';
% 
%     % Smooth EPI
%     matlabbatch{4}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
%     matlabbatch{4}.spm.spatial.smooth.fwhm = [6 6 6];
%     matlabbatch{4}.spm.spatial.smooth.dtype = 0;
%     matlabbatch{4}.spm.spatial.smooth.im = 0;
%     matlabbatch{4}.spm.spatial.smooth.prefix = 's6';
% 
%     % Run preprocessing
%     job_name = fullfile(jobdir, 'Post_PAID_Processing1.mat');
%     save(job_name, 'matlabbatch');
%     spm_jobman('run', matlabbatch)
%     clear matlabbatch

    matlabbatch = [];
    % Normalise write
    fname = fullfile(niidir,'anat', sprintf('y_%s', aname(1).name));
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {fname};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = func_4d;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w2';

    % Smooth EPI
    matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{2}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{2}.spm.spatial.smooth.dtype = 0;
    matlabbatch{2}.spm.spatial.smooth.im = 0;
    matlabbatch{2}.spm.spatial.smooth.prefix = 's6';

    % Run preprocessing
    job_name = fullfile(jobdir, 'Post_PAID_Processing2.mat');
    save(job_name, 'matlabbatch');
    spm_jobman('run', matlabbatch)
    clear matlabbatch

% move the spm ps file to current tedana dir
spm_ps_file = dir(fullfile(pwd,'spm_*.ps'));
if ~isempty(spm_ps_file)
    inputPS = spm_ps_file.name;
    outputPDF = sprintf('spm_post_PAID_%s.pdf',string(datetime));
    command = sprintf('%s -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -sOutputFile="%s" "%s"', ghostscriptPath, outputPDF, inputPS);
    system(command);
    movefile(fullfile(pwd,inputPS),jobdir)
    movefile(fullfile(pwd,outputPDF),jobdir)
end

toc
diary off;
movefile(fullfile(pwd,output_file),jobdir)

end
    
