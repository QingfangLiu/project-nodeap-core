
% redo the func conn analysis by using "correct" motion parameter as
% regressors
% results saved in GLM*_updated folders
% goal is to see if changing this will affect func conn values

clear 
close all

studydir = '/Volumes/X9Pro/NODEAP/MRI';
SubIDlist = dir(fullfile(studydir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

nscans = 310; % 310 volumes for ME sequnce
nechos = 3;
TR = 1.5;
TEs = [14.60, 39.04, 63.48];

%%
for subj = 2:nSubIDlist

    sn = SubIDlist(subj).name;

% Right anterior OFC
regions(1).seedMNI = [34, 54, - 14];
regions(1).seedName = 'aOFC_seed_right';
regions(1).targetMNI = [44, 28, 38];
regions(1).targetName = 'aOFC_target_right';

% Right posterior OFC
regions(2).seedMNI = [28, 38, -16];
regions(2).seedName = 'pOFC_seed_right';
regions(2).targetMNI = [46, 38, 14];
regions(2).targetName = 'pOFC_target_right';

subdir = fullfile(studydir,sn);
niidir = fullfile(subdir, 'nii_for_TMS');
medir = fullfile(niidir, 'day0_me'); % path to multi-echo data
functdir = fullfile(niidir, 'functional'); % path to rsfRMI images
unnormdir = fullfile(subdir, 'UnnormedMasks');  %path to native space masks for each subejct

%% SET up GLM for each seed

% get resting state data
fnames = dir(fullfile(functdir,'s6rf*.nii'));
tmpvol = spm_read_vols(spm_vol(fullfile(functdir, fnames(1).name)));
tmpvol = zeros([size(tmpvol), length(fnames)]);
ffiles = [];
for i = 1:length(fnames)
    ffiles = [ffiles; {fullfile(functdir, fnames(i).name)}];
    tmpvol(:,:,:,i) = spm_read_vols(spm_vol_nifti(ffiles{i}));
end

% using motion parameter from me analysis
nname = dir(fullfile(medir,'rp*.txt'));
nfile = fullfile(medir, nname(1).name); % realignment parameters

for i=1:length(regions)
    matlabbatch = [];
    
    [seedvol, seedXYZ] = spm_read_vols(spm_vol_nifti(fullfile(unnormdir,sprintf('rw%s.nii', regions(i).seedName))));
    seedidx = find(seedvol==1);
    
    % make GLM folder (_updated)
    modeldir = fullfile(niidir,sprintf('GLM_%s_updated',regions(i).seedName));
    
    if ~exist(modeldir,'dir')
        mkdir(modeldir);
    end
    
    % make vectors of mean seed region activity
    seed = zeros(nscans,1);
    for ii = 1:nscans
        dd = tmpvol(:,:,:,ii);
        seed(ii) = nanmean(dd(seedidx));
    end
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {modeldir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = ffiles;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress.name = regions(i).seedName;
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress.val = seed; %  timeseries of mean seed region activity
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {nfile};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    % estimate  model
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % do +1 contrast (t-test) on the seed region regressor
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = regions(i).seedName;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
    matlabbatch{3}.spm.stats.con.delete = 0;
    
    % run job
    spm_jobman('run', matlabbatch)
    
end

end
