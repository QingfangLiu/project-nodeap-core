

% this one calculates map with all gm voxels for ALFF and fALFF values
% use 3mm MNI space

clear; clc; close all
studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis');

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

table_data = readtable(fullfile(HomeDir,'ProcessedData/SubConds.xlsx'));
Order = table_data.StimOrder;
StimLoc = table_data.StimLoc;

nscans = 310; % 310 volumes for ME sequnce
nechos = 3;
TR = 1.5;
TEs = [14.60, 39.04, 63.48];
rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

%%
maskpath = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis/Scripts_matlab/GlobalConn/mask');
gm_nii = fullfile(maskpath,'gm_0.1_3mm.nii'); % gray matter mask
hdr = spm_vol(gm_nii);
gm_dat = spm_read_vols(hdr);
gm_idx = find(gm_dat > 0);
n_gm_vox = length(gm_idx);

%%
for subj = 1:nSubIDlist

    sn = SubIDlist(subj).name;
    
    % prepare where to save the ALFF Results
    respath = fullfile(studydir, 'ALFF_PAID', sn);
    if ~exist(respath,'dir')
        mkdir(respath)
    end
    
    for r = 1:n_rest_names % to loop over rest_names
        
        tic
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,sn)
        
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        
        nhdr = spm_vol(fullfile(maskpath,'wm_0.9_3mm.nii'));
        ALFF_name = fullfile(respath, sprintf('%s_ALFF.nii', curr_rest));
        zALFF_name = fullfile(respath, sprintf('%s_zALFF.nii', curr_rest));
        fALFF_name = fullfile(respath, sprintf('%s_fALFF.nii', curr_rest));
        
        if exist(ALFF_name,'file') && exist(fALFF_name,'file') && exist(zALFF_name,'file') 
            continue;  % skip if both exist
        end
        
        load(fullfile(studydir, 'GlobalConn_PAID', sn, curr_rest, 'tc_filtered_3mm.mat')); 
        
        sampleInterval = TR; 
        lowFreq = 0.01; % Lower bound of frequency range in Hz
        highFreq = 0.08; % Upper bound of frequency range in Hz

        ALFF = zeros(1, n_gm_vox);          % Initialize the ALFF array
        zALFF = zeros(1, n_gm_vox);          % Initialize the zALFF array
        fALFF = zeros(1, n_gm_vox);         % Initialize the fALFF array
        pool = gcp();                       % This will start a parallel pool if one is not already started
        parfor i = 1:n_gm_vox
            [ALFF(i),fALFF(i)] = computePowerMeasures(dat(:, i), sampleInterval, lowFreq, highFreq);
            if ~mod(i,5000)
                fprintf('\t\t%d of %d\n', i, n_gm_vox);
            end
        end
        
        % save ALFF
        vol = zeros(hdr.dim);
        vol(gm_idx) = ALFF;
        nhdr.fname = ALFF_name;
        spm_write_vol(nhdr, vol); 
        
        % also calculate zALFF & save
        % note z-scoring is done for each session across all gm voxels
        vol = zeros(hdr.dim);
        vol(gm_idx) = zscore(ALFF);
        nhdr.fname = zALFF_name;
        spm_write_vol(nhdr, vol); 
        
        % save fALFF
        vol = zeros(hdr.dim);
        vol(gm_idx) = fALFF;
        nhdr.fname = fALFF_name;
        spm_write_vol(nhdr, vol); 
        
        toc
            
    end
end


%% Compute difference map & do smoothing for each subject
do_smooth = 1; % whether to do smooth
sk = 6; % smoothing kernel size (for individual contrast map before aggregation)
do_compute_diff = 1; % whether to calculate difference
agg_contrast_map = 1; % whether to aggregate contrast maps

contrast_names = {'Day1_c_vs_s','Day2_c_vs_s','Day12_c_vs_s'}; % compare cTBS vs sham, on Day1, Day2, or both
n_contrast = length(contrast_names);

Option = 'fALFF';    % to do this with ALFF or fALFF or zALFF

if do_compute_diff
for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;
    
    fprintf('computing alff maps from %s\n',SubID)
    
    % prepare where to save the connResults
    respath = fullfile(studydir, 'ALFF_PAID', SubID);
    outpath = fullfile(respath, Option);
    if ~exist(outpath,'dir')
        mkdir(outpath)
    end
    image_dim = size(spm_read_vols(spm_vol(fullfile(maskpath,'wm_0.9_3mm.nii'))));
    map = nan([image_dim,n_rest_names]);
    
    for r = 1:n_rest_names
        curr_rest = rest_names{r};
        if count_table(subj,r)<1
            continue;
        end
        map(:,:,:,r) = spm_read_vols(spm_vol(fullfile(respath,sprintf('%s_%s.nii', curr_rest, Option))));
    end
     
    if Order(subj)==123 || Order(subj)==132 
        out_Day1_c = map(:,:,:,2);
        out_Day1_s = cat(4, map(:,:,:,4), map(:,:,:,6)); 
        out_Day1_s_mean = nanmean(cat(4, map(:,:,:,4), map(:,:,:,6)),4); 
    elseif Order(subj)==213 || Order(subj)==312
        out_Day1_c = map(:,:,:,4);
        out_Day1_s = cat(4, map(:,:,:,2), map(:,:,:,6));
        out_Day1_s_mean = nanmean(cat(4, map(:,:,:,2), map(:,:,:,6)),4);
    elseif Order(subj)==231 || Order(subj)==321
        out_Day1_c = map(:,:,:,6);
        out_Day1_s = cat(4, map(:,:,:,2), map(:,:,:,4));
        out_Day1_s_mean = nanmean(cat(4, map(:,:,:,2), map(:,:,:,4)),4);
    end
    
    if Order(subj)==213 || Order(subj)==231
        out_Day2_c = map(:,:,:,3);
        out_Day2_s = cat(4, map(:,:,:,5), map(:,:,:,7)); 
        out_Day2_s_mean = nanmean(cat(4, map(:,:,:,5), map(:,:,:,7)),4); 
    elseif Order(subj)==123 || Order(subj)==321
        out_Day2_c = map(:,:,:,5);
        out_Day2_s = cat(4, map(:,:,:,3), map(:,:,:,7)); 
        out_Day2_s_mean = nanmean(cat(4, map(:,:,:,3), map(:,:,:,7)),4); 
    elseif Order(subj)==132 || Order(subj)==312
        out_Day2_c = map(:,:,:,7);
        out_Day2_s = cat(4, map(:,:,:,3), map(:,:,:,5)); 
        out_Day2_s_mean = nanmean(cat(4, map(:,:,:,3), map(:,:,:,5)),4); 
    end
    
    out_Day12_c = nanmean(cat(4, out_Day1_c, out_Day2_c), 4); 
    out_Day12_s = nanmean(cat(4, out_Day1_s, out_Day2_s), 4); 
    
    % 1: Day1_c_vs_s
    hdr.fname = fullfile(outpath,'Day1_c_vs_s.nii');
    spm_write_vol(hdr, out_Day1_c - out_Day1_s_mean);
    
    % 2: Day2_c_vs_s
    hdr.fname = fullfile(outpath,'Day2_c_vs_s.nii');
    spm_write_vol(hdr, out_Day2_c - out_Day2_s_mean);
    
    % 3: Day12_c_vs_s
    hdr.fname = fullfile(outpath,'Day12_c_vs_s.nii');
    spm_write_vol(hdr, out_Day12_c - out_Day12_s);
    
    % Smooth
    if do_smooth
        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {fullfile(outpath,'Day1_c_vs_s.nii')
                                                  fullfile(outpath,'Day2_c_vs_s.nii')
                                                  fullfile(outpath,'Day12_c_vs_s.nii')};
        matlabbatch{1}.spm.spatial.smooth.fwhm = ones(1,3)*sk;
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's6';
        spm_jobman('run', matlabbatch);
        clear matlabbatch
    end
end % subject
end

%% aggregate contrast maps across all subjects
    
if agg_contrast_map
    for c = 1:n_contrast
    % where to save aggregate contrast map results
    pathname = fullfile(studydir, 'ALFF_PAID', 'Agg', Option, contrast_names{c});
    if ~exist(pathname, 'dir')
        mkdir(pathname);
    end

    % check if t-test has been run in this folder
    if any(size(dir(fullfile(pathname,'beta*.nii')),1))
        continue;
    end
    
    ffiles = cell(nSubIDlist,1); 
    for subj = 1:nSubIDlist
        SubID = SubIDlist(subj).name;
        filename = sprintf('%s.nii', contrast_names{c});
        filename = sprintf('s6%s,1',filename); % use smoothed images
        ffiles{subj} = fullfile(studydir, 'ALFF_PAID', SubID, Option, filename);
    end
            
    clear matlabbatch
    matlabbatch{1}.spm.stats.factorial_design.dir = {pathname}; 
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = ffiles;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {fullfile(maskpath,'gm_0.1_3mm.nii')}; % use the gray matter mask
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % contrast
    matlabbatch{3}.spm.stats.con.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = sprintf('positive_%s', contrast_names{c});
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = sprintf('negative_%s', contrast_names{c});
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 1; % delete existing contrasts
    
    % save job
    jobname = fullfile(pathname, 'Job_agg.mat');
    save(jobname, 'matlabbatch')
    
    %run batch
    spm_jobman('run', matlabbatch);  
    clear matlabbatch
    
    end % end loop of contrasts
end


%% aggregate contrast maps across subgroups of subjects
% based on stim location of aOFC or pOFC
 
regions = {'aOFC','pOFC'};
if agg_contrast_map
    for c = 1:n_contrast
    for region = 1:2
    
    roi = regions{region};
    % where to save aggregate contrast map results
    pathname = fullfile(studydir, 'ALFF_PAID', 'Agg', Option, sprintf('%s_%s',contrast_names{c},roi));
    if ~exist(pathname, 'dir')
        mkdir(pathname);
    end

    % check if t-test has been run in this folder
    if any(size(dir(fullfile(pathname,'beta*.nii')),1))
        continue;
    end
    
    if strcmp(roi,'aOFC')
        tmploc = strcmp(StimLoc,'Anterior'); % find aOFC stim subjects
    else
        tmploc = strcmp(StimLoc,'Posterior'); % find pOFC stim subjects
    end
    
    ntmpsub = sum(tmploc);
    ffiles = cell(ntmpsub,1);
    subctr = 0;
    for subj = find(tmploc==1)'
        subctr = subctr + 1;
        SubID = SubIDlist(subj).name;
        filename = sprintf('%s.nii', contrast_names{c});
        filename = sprintf('s6%s,1',filename); % use smoothed images
        ffiles{subctr} = fullfile(studydir, 'ALFF_PAID', SubID, Option, filename);
    end 
      
    clear matlabbatch
    matlabbatch{1}.spm.stats.factorial_design.dir = {pathname}; 
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = ffiles;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {fullfile(maskpath,'gm_0.1_3mm.nii')}; % use the gray matter mask
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % contrast
    matlabbatch{3}.spm.stats.con.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = sprintf('positive_%s', contrast_names{c});
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = sprintf('negative_%s', contrast_names{c});
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 1; % delete existing contrasts
    
    % save job
    jobname = fullfile(pathname, 'Job_agg.mat');
    save(jobname, 'matlabbatch')
    
    %run batch
    spm_jobman('run', matlabbatch);  
    clear matlabbatch
    
    end
    end % end loop of contrasts
end

