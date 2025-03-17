
% conduct functional conn analysis using data from PAID and other
% following analyses
% do this for rs-fMRI data after each TMS
% using data in native space (similar to our processing Day 0 data)

clear
clc
close all

studydir = '/Volumes/X9Pro/NODEAP/MRI';
SubIDlist = dir(fullfile(studydir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

nscans = 310; % 310 volumes for ME sequnce
nechos = 3;
TR = 1.5;
TEs = [14.60, 39.04, 63.48];
rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

%%
for subj = 47

    sn = SubIDlist(subj).name;
    subdir = fullfile(studydir,sn);
    jobdir = fullfile(subdir, 'jobs');
    niidir = fullfile(subdir, 'nii');
    unnormdir = fullfile(subdir, 'UnnormedMasks');  %path to native space masks for each subejct

for r = 3:n_rest_names % to loop over rest_names
        
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,subdir)
        functdir = fullfile(niidir,'functional',curr_rest); % path to rsfRMI images
        func_conn_dir = fullfile(subdir, 'FuncConn_PAID',curr_rest); % path to save the GLM result from PAID result
        curr_me_dir = fullfile(niidir,sprintf('%s_me',curr_rest));
    
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
    
        % skip if it has been run
%         if exist(fullfile(func_conn_dir,'stimulation_coordinates.mat'),'file')
%             disp('Processing has been done, skipping this iteration ...')
%             continue;
%         end
    
clear regions
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
            
%% SET up GLM for each seed

% get resting state data
% fnames = dir(fullfile(functdir,'s6f*.nii')); % smoothed data in native space
% ffile = fullfile(functdir,fnames(1).name);
% tmpvol = spm_read_vols(spm_vol(ffile)); % 4D array, 104*104*60*310, this step takes a bit of time

% reorganize the 4D array so spm reads it
% ffiles = spm_select('expand', ffile); % expand the volumes from a 4D nifti file
% ffiles = cellstr(ffiles); 
% 
% % realignment parameters (from multi-echo folder)
% nname = dir(fullfile(curr_me_dir,'rp*.txt'));
% nfile = fullfile(curr_me_dir, nname(1).name); 

for i=1:length(regions)
    
    [seedvol, seedXYZ] = spm_read_vols(spm_vol_nifti(fullfile(unnormdir,sprintf('rw%s.nii', regions(i).seedName))));
    seedidx = find(seedvol==1);
    
    % make GLM folder
    modeldir = fullfile(func_conn_dir,sprintf('GLM_%s',regions(i).seedName));
    
    if ~exist(modeldir,'dir')
        mkdir(modeldir);
    end
    
    % make vectors of mean seed region activity
%     seed = zeros(nscans,1);
%     for ii = 1:nscans
%         dd = tmpvol(:,:,:,ii);
%         seed(ii) = nanmean(dd(seedidx));
%     end
    
    matlabbatch = [];
%     matlabbatch{1}.spm.stats.fmri_spec.dir = {modeldir};
%     matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
%     matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
%     matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
%     matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
%     matlabbatch{1}.spm.stats.fmri_spec.sess.scans = ffiles;
%     matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
%     matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
%     matlabbatch{1}.spm.stats.fmri_spec.sess.regress.name = regions(i).seedName;
%     matlabbatch{1}.spm.stats.fmri_spec.sess.regress.val = seed; %  timeseries of mean seed region activity
%     matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {nfile};
%     matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
%     matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
%     matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
%     matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
%     matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
%     matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
%     matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
%     matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    % estimate  model
%     matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
%     matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
%     matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % do +1 contrast (t-test) on the seed region regressor
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(modeldir,'SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = regions(i).seedName;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
    matlabbatch{1}.spm.stats.con.delete = 1;
    
    %fname = fullfile(jobdir, sprintf('FuncConn_GLM_%s_%s_PAID',regions(i).seedName,curr_rest));
    %save(fname, 'matlabbatch');
    
    % run job
    spm_jobman('run', matlabbatch)
    
    %% find coordinates of maximum connectivity in target regions
    [targetvol, targetXYZ] = spm_read_vols(spm_vol_nifti(fullfile(unnormdir,sprintf('rw%s.nii',regions(i).targetName))));
    targetidx = find(targetvol==1);
    
    convol = spm_read_vols(spm_vol_nifti(fullfile(modeldir,'spmT_0001.nii')));
    convals = convol(targetidx);
    
    StimCoord = round(targetXYZ(:,targetidx(convals==max(convals)))');
    sprintf('coordinate for %s: x=%01d, y=%01d, z=%01d', regions(i).seedName, StimCoord(1), StimCoord(2), StimCoord(3));
    
    regions(i).StimCoord = StimCoord;
end
    % save coordinates from both GLMs
    save(fullfile(func_conn_dir, 'stimulation_coordinates.mat'), 'regions')
    
end
end
