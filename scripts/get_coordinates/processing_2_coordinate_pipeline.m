
clear all
close all

do_combine_echos = 1;
do_select_coordinate =1;

project_folder = '/Users/liuq13/project-nodeap-core';
maskdir     = fullfile(project_folder, 'atlas_masks', 'ofc_connectivity_masks');
tmpl = fullfile(project_folder, 'tissue_masks', 'rTPM_point1_NoCerebellum.nii');

dat_folder   = '/Volumes/X9Pro/NODEAP/MRI';

SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);       % keep directories only
nSub      = numel(SubIDlist);

nscans = 310; % 310 volumes for ME sequnce
nechos = 3;
TR = 1.5;
TEs = [14.60, 39.04, 63.48];

%%
for subj = 1:nSub

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

subdir = fullfile(dat_folder,sn);
niidir = fullfile(subdir, 'nii_for_TMS');
medir = fullfile(niidir, 'day0_me'); % path to multi-echo data
functdir = fullfile(niidir, 'functional'); % path to rsfRMI images
anatdir = fullfile(niidir, 'anat'); % path to anatomy
jobdir = fullfile(niidir, 'jobs');
unnormdir = fullfile(subdir, 'UnnormedMasks');  %path to native space masks for each subejct

cd(niidir)
spm fmri
if do_combine_echos
    
    %%% realign and reslice all volumes within each echo
    rfilename = [];
    for j=1:nechos
        matlabbatch = [];
        filename = [];
        
        n = dir(fullfile(medir, sprintf('D0_rest*e%01d.nii',j)));
        for i = 1:nscans
            filename{1}{i,1} = sprintf('%s,%01d', fullfile(medir, n(1).name),i);
            rfilename{j}{i,1} = sprintf('%s,%01d', fullfile(medir, sprintf('r%s', n(1).name)),i);
        end
        
        matlabbatch{1}.spm.spatial.realign.estwrite.data = filename;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % register to first !!!
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
        
        % save and run job
        if ~isfolder(jobdir)
            mkdir(jobdir);
        end
        fname = fullfile(jobdir, sprintf('me_realignment_e%01d.mat',j));
        save(fname, 'matlabbatch');
        spm_jobman('run', matlabbatch);
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
        spm_write_vol(hdr2, tSNR{j});
        
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
        spm_write_vol(hdr2, w{j});
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


%% make spherical ROIs surrounding the coordinates, constrained by the tissue probability map
if ~exist(maskdir,'dir')
    mkdir(maskdir);
end

for i = 1:length(regions)
    seedmask{i} = fullfile(maskdir, sprintf('%s.nii',regions(i).seedName));
    targetmask{i} = fullfile(maskdir, sprintf('%s.nii',regions(i).targetName));
end

%check if spherical ROIs in MNI space have already been made
if ~exist(seedmask{1}) && ~exist(targetmask{1})
    
    rget = 4;
    
    % load gray matter mask
    tpmhdr = spm_vol_nifti(tmpl);
    [tpmvol,XYZ] = spm_read_vols(tpmhdr);
    tpmidx = find(tpmvol>0);
    sz = size(tpmvol);
    [MX, MY, MZ] = ndgrid(1:sz(1),1:sz(2),1:sz(3));
    ref_vox = round(sz/2);
    radii = sqrt((MX-ref_vox(1)).^2 + (MY-ref_vox(2)).^2 + (MZ-ref_vox(3)).^2);
    radius_index = find(radii<rget) - sub2ind(sz,ref_vox(1),ref_vox(2),ref_vox(3));
    
    for i = 1:length(regions)
        
        % seed regions
        coordidx = find(ismember(XYZ',regions(i).seedMNI,'rows')>0);
        idx = intersect(radius_index + coordidx, tpmidx);
        newhdr = tpmhdr;
        newvol = zeros(sz);
        newhdr.fname =  seedmask{i};
        newvol(idx) = 1;
        spm_write_vol(newhdr,newvol);
        
        % target regions
        coordidx = find(ismember(XYZ',regions(i).targetMNI,'rows')>0);
        idx = intersect(radius_index + coordidx, tpmidx);
        newhdr = tpmhdr;
        newvol = zeros(sz);
        newhdr.fname = targetmask{i};
        newvol(idx) = 1;
        spm_write_vol(newhdr,newvol);
    end
else
    fprintf('\nSeed and Target masks already made\n');
end


%% Do functional image preprocessing

% get anatomical image
aname = dir(fullfile(anatdir,'D0_T1*.nii'));
afile = fullfile(anatdir,aname(1).name);

% Check if the data has already been preprocessed
if length(dir(fullfile(functdir,'s6rf*.nii'))) ~= nscans
    
    % get functional images
    fnames = dir(fullfile(functdir,'f*.nii'));
    
    ffiles = [];
    for i = 1:length(fnames)
        ffiles = [ffiles; {fullfile(functdir,fnames(i).name)}];
    end
    
    % Realignment
    matlabbatch{1}.spm.spatial.realign.estwrite.data = {ffiles};
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    % Coregister mean realigned EPI to T1, and bring along the realigned EPIs
    matlabbatch{2}.spm.spatial.coreg.estimate.ref = {afile};
    matlabbatch{2}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{2}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
    matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    % Normalize T1 to generate warping parameters
    matlabbatch{3}.spm.spatial.normalise.estwrite.subj.vol = {afile};
    matlabbatch{3}.spm.spatial.normalise.estwrite.subj.resample = {afile};
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.tpm =  {fullfile(spm('dir'),'tpm','TPM.nii')};
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
    matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
    
    % Smooth EPI
    matlabbatch{4}.spm.spatial.smooth.data(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{4}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{4}.spm.spatial.smooth.dtype = 0;
    matlabbatch{4}.spm.spatial.smooth.im = 0;
    matlabbatch{4}.spm.spatial.smooth.prefix = 's6';
    
    % make inverse deformation field
    matlabbatch{5}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(anatdir,sprintf('y_%s',aname(1).name))};
    matlabbatch{5}.spm.util.defs.comp{1}.inv.space = {fullfile(functdir,sprintf('s6r%s',fnames(1).name))};
    matlabbatch{5}.spm.util.defs.out{1}.savedef.ofname = 'inverse_deformation';
    matlabbatch{5}.spm.util.defs.out{1}.savedef.savedir.saveusr = {anatdir};
    
    % Run preprocessing
    fname = fullfile(jobdir, 'func_data_processing.mat');
    save(fname, 'matlabbatch');
        
    spm_jobman('run', matlabbatch)
    clear matlabbatch
    
else
    fprintf('\nData already preprocessed\n');
end

%% un-norm and reslice spherical seed and target ROIs
fnames = dir(fullfile(functdir,'s6rf*.nii'));% get functional images
unnormfile = fullfile(anatdir,'y_inverse_deformation.nii');% inverse deformation map

matlabbatch = [];
% unnormalize masks
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {unnormfile};
ii=0;
for i=1:length(regions)
    ii=ii+1;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample{ii,1} = strcat(seedmask{i},',1');
    ii=ii+1;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample{ii,1} = strcat(targetmask{i},',1');
end
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    78 112 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

% resclice unnormalized masks
matlabbatch{2}.spm.spatial.coreg.write.ref = {fullfile(functdir,fnames(1).name)};
matlabbatch{2}.spm.spatial.coreg.write.source(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';

fname = fullfile(jobdir, 'unnormalize_masks.mat');
save(fname, 'matlabbatch');
spm_jobman('run', matlabbatch)

% move un-norm files to the subject mask folder
if ~exist(unnormdir,'dir')
    mkdir(unnormdir);
end

for i=1:length(regions)
    movefile(fullfile(maskdir, sprintf('rw%s.nii', regions(i).seedName)), unnormdir);
    movefile(fullfile(maskdir, sprintf('w%s.nii', regions(i).seedName)), unnormdir);
    movefile(fullfile(maskdir, sprintf('rw%s.nii', regions(i).targetName)), unnormdir);
    movefile(fullfile(maskdir, sprintf('w%s.nii', regions(i).targetName)), unnormdir);
end


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

nname = dir(fullfile(functdir,'rp*.txt'));
nfile = fullfile(functdir, nname(1).name); % realignment parameters

for i=1:length(regions)
    matlabbatch = [];
    
    [seedvol, seedXYZ] = spm_read_vols(spm_vol_nifti(fullfile(unnormdir,sprintf('rw%s.nii', regions(i).seedName))));
    seedidx = find(seedvol==1);
    
    % make GLM folder
    modeldir = fullfile(subdir,sprintf('GLM_%s',regions(i).seedName));
    
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
    
    fname = fullfile(jobdir, sprintf('GLM_%s',regions(i).seedName));
    save(fname, 'matlabbatch');
    
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

for i=1:length(regions)
    fprintf('coordinate for %s: x=%01d, y=%01d, z=%01d \n', regions(i).seedName, regions(i).StimCoord(1), regions(i).StimCoord(2), regions(i).StimCoord(3))
end

% save coordinates
save(fullfile(subdir, 'stimulation_coordinates.mat'), 'regions')

end
