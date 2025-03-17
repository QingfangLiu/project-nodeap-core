
% redo 
% by NOT doing 2nd time of realignment
% use 's6f*.nii' data to do func conn analysis

% but the previous f* files have been contaminated, so also redo the
% weighted combination

%processing_coordinate_pipeline_NODEAP
clear all
close all

do_combine_echos = 1;
do_select_coordinate =1;

studydir = '/Volumes/X9Pro/NODEAP/MRI';
SubIDlist = dir(fullfile(studydir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

nscans = 310; % 310 volumes for ME sequnce
nechos = 3;
TR = 1.5;
TEs = [14.60, 39.04, 63.48];


%%
for subj = 36:nSubIDlist

    sn = SubIDlist(subj).name;
    fprintf('%s\n',sn)
    
    tic

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
functdir = fullfile(niidir, 'functional_test'); % path to rsfRMI images
anatdir = fullfile(niidir, 'anat'); % path to anatomy
jobdir = fullfile(niidir, 'jobs');
maskdir = fullfile(studydir, 'ConnectivityMasks'); % path to masks for study
unnormdir = fullfile(subdir, 'UnnormedMasks');  %path to native space masks for each subejct

if do_combine_echos
    
    %%% realign and reslice all volumes within each echo
    rfilename = [];
    for j=1:nechos
        
        filename = [];
        
        n = dir(fullfile(medir, sprintf('D0_rest*e%01d.nii',j)));
        for i = 1:nscans
            filename{1}{i,1} = sprintf('%s,%01d', fullfile(medir, n(1).name),i);
            rfilename{j}{i,1} = sprintf('%s,%01d', fullfile(medir, sprintf('r%s', n(1).name)),i);
        end
    end
    
    %%% COMBINE ECHOS
    % compute tSNR maps per echo
    hdr2 = spm_vol(filename{1}{1,1}); % get one header
    hdr2.dt = [16,0];
    tSNR = [];
    tSNR_TE = [];
    for j=1:nechos
        
        tmp = [];
        % load data
        x = [];
        for i= 1:nscans
            hdr = spm_vol(rfilename{j}{i,1});
            x(:,:,:,i) = spm_read_vols(hdr);
        end
        
        % compute and save tSNR map
        tSNR{j} = mean(x,4)./std(x, [], 4);
        hdr2.fname = fullfile(niidir, sprintf('tSNR_echo-%01d.nii', j));
%         spm_write_vol(hdr2, tSNR{j});
        
        % multiply tSNR by TE for each echo
        tSNR_TE(:,:,:,j) = tSNR{j}.*TEs(j);
    end
    clear x
    
    % compute wight (w)
    w = [];
    for j=1:nechos
        
        % normalize tSNR*TE maps by sum of tSNR*TE's across echos
        w{j} = tSNR_TE(:,:,:,j)./sum(tSNR_TE,4);
        
        % save tSNR map
        hdr2.fname = fullfile(niidir, sprintf('w_tSNR_TE_echo-%01d.nii', j));
%         spm_write_vol(hdr2, w{j});
    end
    
    % combine echos using w
    if ~isdir(functdir)
        mkdir(functdir)
    end
    
    for i=1:nscans
        y = zeros(size(w{j}));
        for j=1:nechos
            hdr = spm_vol(rfilename{j}{i,1});
            x = spm_read_vols(hdr);
            y = y + (x.*w{j});
        end
        
        % save combined volume
        hdr.n = [1 1];
        hdr.fname = fullfile(functdir, sprintf('fvol_%03d.nii', i));
        spm_write_vol(hdr, y);
    end
end % end combine echos


%% Do functional image preprocessing

% Check if the data has already been preprocessed
if length(dir(fullfile(functdir,'s6f*.nii'))) ~= nscans
    
    % get functional images
    fnames = dir(fullfile(functdir,'f*.nii'));
    
    ffiles = [];
    for i = 1:length(fnames)
        ffiles = [ffiles; {fullfile(functdir,fnames(i).name)}];
    end
    
    clear matlabbatch
    
    % Smooth EPI
    matlabbatch{1}.spm.spatial.smooth.data = ffiles;
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's6';
    
    spm_jobman('run', matlabbatch)
    clear matlabbatch
    
else
    fprintf('\nData already preprocessed\n');
end

%% SET up GLM for each seed

% get resting state data
fnames = dir(fullfile(functdir,'s6f*.nii'));
tmpvol = spm_read_vols(spm_vol(fullfile(functdir, fnames(1).name)));
tmpvol = zeros([size(tmpvol), length(fnames)]);
ffiles = [];
for i = 1:length(fnames)
    ffiles = [ffiles; {fullfile(functdir, fnames(i).name)}];
    tmpvol(:,:,:,i) = spm_read_vols(spm_vol_nifti(ffiles{i}));
end

nname = dir(fullfile(medir,'rp*.txt'));
nfile = fullfile(medir, nname(1).name); % realignment parameters

for i=1:length(regions)
    matlabbatch = [];
    
    [seedvol, seedXYZ] = spm_read_vols(spm_vol_nifti(fullfile(unnormdir,sprintf('rw%s.nii', regions(i).seedName))));
    seedidx = find(seedvol==1);
    
    % make GLM folder
    modeldir = fullfile(niidir,sprintf('GLM_%s_updated2',regions(i).seedName));
    
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

toc


end


