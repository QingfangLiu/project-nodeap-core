
clc; clear;
userHome = getenv('HOME');
maskpath = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis/Scripts_matlab/GlobalConn/mask');
studydir = '/Volumes/X9Pro/NODEAP/MRI';
SubIDlist = dir(fullfile(studydir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

HomeDir = '/Volumes/X9Pro/NODEAP';
file_path = fullfile(HomeDir,'NODEAP_DataCollectionSheet.xlsx');
table_data = readtable(file_path,'NumHeaderLines',0); % need to specify header line to properly read the sheet (ignore warning)
Order = table_data.Order;
StimLoc = table_data.Anterior_Posterior;

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

Option = 'Unsigned'; % Unsigned or Signed or Squared
MRIcount = readtable('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx',"ReadRowNames",true);
do_smooth = 1; % whether to do smooth
sk = 6; % smoothing kernel size (for individual contrast map before aggregation)
agg_contrast_map = 1; % whether to aggregate contrast maps

contrast_names = {'Day1_c_vs_s','Day2_c_vs_s','Day12_c_vs_s'};
n_contrast = length(contrast_names);

%% load masks
hdr = spm_vol(fullfile(maskpath,'gm_0.1_3mm.nii'));
gm = spm_read_vols(hdr);
gm_idx = find(gm > 0);
nvox = length(gm_idx); % number of voxels

%%
for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;
    SubDir = fullfile(studydir,SubID);
    
    % prepare where to save the connResults
    respath = fullfile(studydir, 'GlobalConn', SubID, Option);
    if ~exist(respath,'dir')
        mkdir(respath)
    end
    
    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
        tic
        
        if MRIcount{SubID,curr_rest}==0
            continue;
        end
        
        if MRIcount{SubID,curr_rest}==0.5
            nscans = 205;
        else
            nscans = 310;
        end
        
        fprintf('do connectedness from %s of %s\n',curr_rest,SubID)
        
        nhdr = spm_vol(fullfile(maskpath,'wm_0.9_3mm.nii'));
        nhdr.fname = fullfile(respath, sprintf('%s.nii', curr_rest));
        
        if ~exist(nhdr.fname,'file') % check if it has been created
            % load extracted and filtered scan data
            load(fullfile(studydir, 'GlobalConn', SubID, curr_rest, 'tc_filtered_3mm.mat')); 
            dat = zscore(dat); % zcore for convenient Pearson correlation calculation
            cors = zeros(1,nvox); % initialize a vector for each voxel
        
            parfor i = 1:nvox
                tmp = (dat(:,i)' * dat)./ nscans; % Pearson correlation on z-scored values (with each voxel)

                if strcmp(Option,'Unsigned')
                    cors(i) = sum(abs(atanh(tmp(tmp<1))))/nvox; %UNSIGNED recommended (average across voxels, transformed first)
                elseif strcmp(Option,'Signed')
                    cors(i) = sum(atanh(tmp(tmp<1)))/nvox; % SIGNED
                else
                    cors(i) = sum((atanh(tmp(tmp<1))).^2)/nvox; % Squared
                end

                if ~mod(i,5000)
                    fprintf('\t\t%d of %d\n', i, nvox);
                end

            end

            vol = zeros(hdr.dim);
            vol(gm_idx) = cors;
            spm_write_vol(nhdr, vol); 
        
        end
        
        toc
        
    end % end of sessions
    
end % end of subjects


%% Compute difference map & do smoothing for each subject

for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;
    SubDir = fullfile(studydir,SubID);
    
    fprintf('computing connectedness maps from %s\n',SubID)
    
    % prepare where to save the connResults
    respath = fullfile(studydir, 'GlobalConn', SubID, Option);
    image_dim = size(spm_read_vols(spm_vol(fullfile(maskpath,'wm_0.9_3mm.nii'))));
    map = nan([image_dim,n_rest_names]);
    
    for r = 1:n_rest_names
        curr_rest = rest_names{r};
        if MRIcount{SubID,curr_rest}==0
            continue;
        end
        map(:,:,:,r) = spm_read_vols(spm_vol(fullfile(respath,sprintf('%s.nii', curr_rest))));
    end
     
    if Order(subj)==123 || Order(subj)==132
        out_Day1_c = map(:,:,:,2);
        out_Day1_s = nanmean(cat(4, map(:,:,:,4), map(:,:,:,6)), 4); 
    elseif Order(subj)==213 || Order(subj)==312
        out_Day1_c = map(:,:,:,4);
        out_Day1_s = nanmean(cat(4, map(:,:,:,2), map(:,:,:,6)), 4); 
    elseif Order(subj)==231 || Order(subj)==321
        out_Day1_c = map(:,:,:,6);
        out_Day1_s = nanmean(cat(4, map(:,:,:,2), map(:,:,:,4)), 4); 
    end
    
    if Order(subj)==213 || Order(subj)==231
        out_Day2_c = map(:,:,:,3);
        out_Day2_s = nanmean(cat(4, map(:,:,:,5), map(:,:,:,7)), 4); 
    elseif Order(subj)==123 || Order(subj)==321
        out_Day2_c = map(:,:,:,5);
        out_Day2_s = nanmean(cat(4, map(:,:,:,3), map(:,:,:,7)), 4); 
    elseif Order(subj)==132 || Order(subj)==312
        out_Day2_c = map(:,:,:,7);
        out_Day2_s = nanmean(cat(4, map(:,:,:,3), map(:,:,:,5)), 4); 
    end
    
    out_Day12_c = nanmean(cat(4, out_Day1_c, out_Day2_c), 4); 
    out_Day12_s = nanmean(cat(4, out_Day1_s, out_Day2_s), 4); 
       
    % 1: Day1_c_vs_s
    hdr.fname = fullfile(respath,'Day1_c_vs_s.nii');
    spm_write_vol(hdr, out_Day1_c - out_Day1_s);
    
    % 2: Day2_c_vs_s
    hdr.fname = fullfile(respath,'Day2_c_vs_s.nii');
    spm_write_vol(hdr, out_Day2_c - out_Day2_s);
    
    % 3: Day12_c_vs_s
    hdr.fname = fullfile(respath,'Day12_c_vs_s.nii');
    spm_write_vol(hdr, out_Day12_c - out_Day12_s);
    
    % Smooth
    if do_smooth
        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {fullfile(respath,'Day1_c_vs_s.nii')
                                                  fullfile(respath,'Day2_c_vs_s.nii')
                                                  fullfile(respath,'Day12_c_vs_s.nii')};
        matlabbatch{1}.spm.spatial.smooth.fwhm = ones(1,3)*sk;
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's6';
        spm_jobman('run', matlabbatch);
        clear matlabbatch
    end
end % subject

%% aggregate contrast maps across all subjects
    
if agg_contrast_map
    for c = 1:n_contrast
    % where to save aggregate contrast map results
    pathname = fullfile(studydir, 'GlobalConn', 'Agg', contrast_names{c});
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
        ffiles{subj} = fullfile(studydir, 'GlobalConn', SubID, Option, filename);
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
 
regions = {'aOFC','pOFC'};
if agg_contrast_map
    for c = 1:n_contrast
    for region = 1:2
    
    roi = regions{region};
    % where to save aggregate contrast map results
    pathname = fullfile(studydir, 'GlobalConn', 'Agg', sprintf('%s_%s',contrast_names{c},roi));
    if ~exist(pathname, 'dir')
        mkdir(pathname);
    end

    % check if t-test has been run in this folder
    if any(size(dir(fullfile(pathname,'beta*.nii')),1))
        continue;
    end
    
    if strcmp(roi,'aOFC')
        tmploc = strcmp(StimLoc,'Anterior');
    else
        tmploc = strcmp(StimLoc,'Posterior');
    end
    
    ntmpsub = sum(tmploc);
    ffiles = cell(ntmpsub,1);
    subctr = 0;
    for subj = find(tmploc==1)'
        subctr = subctr + 1;
        SubID = SubIDlist(subj).name;
        filename = sprintf('%s.nii', contrast_names{c});
        filename = sprintf('s6%s,1',filename); % use smoothed images
        ffiles{subctr} = fullfile(studydir, 'GlobalConn', SubID, Option, filename);
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




