
% conduct functional conn analysis using data from PAID and other
% following analyses
% do this for rs-fMRI data after each TMS
% using data in MNI space 

clear; clc; close all

studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

FCdir = fullfile(studydir,'FuncConn_PAID_MNI');

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

nscans = 310; % 310 volumes for ME sequnce
nechos = 3;
TR = 1.5;
TEs = [14.60, 39.04, 63.48];
rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

FC_models = {'aOFC_seed','aOFC_stim','pOFC_seed','pOFC_stim'};
n_FC_models = length(FC_models);
masks = {'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'};

userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis');
maskpath = fullfile(HomeDir,'Scripts_matlab/FuncConn/mask');

spm('Defaults','fMRI')
global defaults
defaults.stats.maxmem = 2 ^ 34; % based on RAM = 32GB
defaults.stats.resmem = true; % allow temp files to be kept in memory

gm_mask = fullfile(maskpath,'gm_0.1_2mm.nii');

%%
for subj = 1:nSubIDlist

    sn = SubIDlist(subj).name;
    subdir = fullfile(MRIdir,sn);
    niidir = fullfile(subdir, 'nifti');
    
    jobdir = fullfile(FCdir, sn, 'FCjobs');
    if ~exist(jobdir,'dir')
        mkdir(jobdir)
    end
    
for r = 1:n_rest_names % to loop over rest_names
        
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,sn)
        functdir = fullfile(niidir,'functional',curr_rest); % path to rsfRMI images
        func_conn_dir = fullfile(FCdir, sn, curr_rest); % path to save the GLM result from PAID result
        curr_me_dir = fullfile(niidir,sprintf('%s_me',curr_rest));
    
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
                
%% SET up GLM for each seed

% get resting state data
fnames = dir(fullfile(functdir,'s6w2f*.nii')); % smoothed data in MNI space

ffile = fullfile(functdir,fnames(1).name);
tmpvol = spm_read_vols(spm_vol(ffile)); % 4D array, 79*95*79*310, this step takes a bit of time

% reorganize the 4D array so spm reads it
ffiles = spm_select('expand', ffile); % expand the volumes from a 4D nifti file
ffiles = cellstr(ffiles); 

% use the pre-calculated nuisance regressor
nfile = fullfile(studydir, 'NRegressor',sn, sprintf('nuisance_regressors_%s.txt',curr_rest));

for i=1:n_FC_models
    
    [targetvol, ~] = spm_read_vols(spm_vol_nifti(fullfile(studydir,'ConnectivityMasks', masks{i})));
    targetidx = find(targetvol==1);
    
    % make GLM folder
    modeldir = fullfile(func_conn_dir,sprintf('GLM_%s',FC_models{i}));
    if ~exist(modeldir,'dir')
        mkdir(modeldir);
    end
    
    % continue to next if this has been run
    if exist(fullfile(modeldir,'spmT_0001.nii'),'file')
        continue;
    end
    
    % make vectors of mean target region activity
    target = zeros(nscans,1);
    for ii = 1:nscans
        dd = tmpvol(:,:,:,ii);
        target(ii) = nanmean(dd(targetidx));
    end
    
    matlabbatch = [];
    matlabbatch{1}.spm.stats.fmri_spec.dir = {modeldir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = ffiles;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress.name = FC_models{i};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress.val = target; %  timeseries of mean target region activity
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {nfile};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 180;  % change this from default 128 to 180
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mask = {gm_mask};
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.1; % change from default 0.8 to 0.1
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'none';  % change from default 'AR(1)' to 'none'
        
    % estimate  model
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % do +1 contrast (t-test) on the target regressor
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = FC_models{i};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
    matlabbatch{3}.spm.stats.con.delete = 0;
    
    % run job
    job_name = fullfile(jobdir, sprintf('FuncConn_PAID_%s_%s.mat',curr_rest,FC_models{i}));
    save(job_name, 'matlabbatch');
    spm_jobman('run', matlabbatch)
    clear matlabbatch
     
end
    
end
end
