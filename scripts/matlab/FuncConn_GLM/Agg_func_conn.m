
% this code tries to aggregate the func conn maps across subjects

% each of the 4 seed ROIs gives one func conn map

% want to get cTBS-sham difference, only from Day 1, Day 2, or both
% and then do one-sample t-test to aggregate
% (without considering potential order effect)

clear; clc; close all

studydir = '/Volumes/X9Pro/NODEAP';
FCdir = fullfile(studydir,'FuncConn_PAID_MNI');

userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis');
maskpath = fullfile(HomeDir,'Scripts_matlab/FuncConn/mask');

SubIDlist = dir(fullfile(FCdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

% read MRI count file
count_table = xlsread(fullfile(studydir,'MRI_func_count.xlsx'));
table_data = readtable(fullfile(HomeDir,'ProcessedData/SubConds.xlsx'));
Order = table_data.StimOrder;

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

FC_models = {'aOFC_seed','aOFC_stim','pOFC_seed','pOFC_stim'};

sk = 6; % smoothing kernel size (for individual contrast map before aggregation)
agg_contrast_map = 1;
contrast_names = {'Day1_c_vs_s','Day2_c_vs_s','Day12_c_vs_s',...
    'Day1_c','Day2_c','Day12_c',...
    'Day1_s','Day2_s','Day12_s',...
    'Day12_both','All','Day0'}; 
n_contrast = length(contrast_names);

%%
for subj = 1:nSubIDlist

    sn = SubIDlist(subj).name;
    fprintf('Calculating contrast maps in %s\n',sn)
        
    for i = 1:length(FC_models)  % for each func conn model
        
     respath = fullfile(FCdir,sn,'Contrasts',FC_models{i});
     if ~exist(respath,'dir')
         mkdir(respath)
     end

    map = nan([79,95,79,n_rest_names]); % dimension in MNI space
    
    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
        % fprintf('Loading beta values in %s of %s\n',curr_rest,sn)
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        modelname = sprintf('GLM_%s',FC_models{i});
        modeldir = fullfile(FCdir,sn,curr_rest,modelname);
        hdr = spm_vol_nifti(fullfile(modeldir,'beta_0001.nii'));  % use beta map
        map(:,:,:,r) = spm_read_vols(hdr); 
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
    
    out_Day12_both = nanmean(cat(4, out_Day1_s, out_Day2_s, out_Day1_c, out_Day2_c), 4); 
    out_all = nanmean(map,4);
    
    out_Day0 = map(:,:,:,1);
    
    % 1: Day1_c_vs_s
    if all(isnan(out_Day1_c(:))) || all(isnan(out_Day1_s_mean(:)))
        disp('Skipping saving because at least one map contains all NAN values.')
    else
        hdr.fname = fullfile(respath,'Day1_c_vs_s.nii');
         if ~exist(hdr.fname,'file')
            spm_write_vol(hdr, out_Day1_c - out_Day1_s_mean);
         end
    end
    
    % 2: Day2_c_vs_s
    if all(isnan(out_Day2_c(:))) || all(isnan(out_Day2_s_mean(:)))
        disp('Skipping saving because at least one map contains all NAN values.')
    else
        hdr.fname = fullfile(respath,'Day2_c_vs_s.nii');
         if ~exist(hdr.fname,'file')
            spm_write_vol(hdr, out_Day2_c - out_Day2_s_mean);
         end
    end
    
    % 3: Day12_c_vs_s
    if all(isnan(out_Day12_c(:))) || all(isnan(out_Day12_s(:)))
        disp('Skipping saving because at least one map contains all NAN values.')
    else
        hdr.fname = fullfile(respath,'Day12_c_vs_s.nii');
         if ~exist(hdr.fname,'file')
            spm_write_vol(hdr, out_Day12_c - out_Day12_s);
         end
    end
    
    % 4: Day1_c
    if all(isnan(out_Day1_c(:)))
        disp('Skipping saving because the map contains all NAN values.')
    else
        hdr.fname = fullfile(respath,'Day1_c.nii');
        if ~exist(hdr.fname,'file')
            spm_write_vol(hdr, out_Day1_c);
        end
    end
    
    % 5: Day2_c
    if all(isnan(out_Day2_c(:))) 
        disp('Skipping saving because the map contains all NAN values.')
    else
        hdr.fname = fullfile(respath,'Day2_c.nii');
        if ~exist(hdr.fname,'file')
            spm_write_vol(hdr, out_Day2_c);
        end
    end
    
    % 6: Day12_c
    if all(isnan(out_Day12_c(:)))
        disp('Skipping saving because the map contains all NAN values.')
    else
        hdr.fname = fullfile(respath,'Day12_c.nii');
        if ~exist(hdr.fname,'file')
            spm_write_vol(hdr, out_Day12_c);
        end
    end
    
    % 7: Day1_s
    if all(isnan(out_Day1_s_mean(:)))
        disp('Skipping saving because the map contains all NAN values.')
    else
        hdr.fname = fullfile(respath,'Day1_s.nii');
        if ~exist(hdr.fname,'file')
            spm_write_vol(hdr, out_Day1_s_mean);
        end
    end
    
    % 8: Day2_s
    if all(isnan(out_Day2_s_mean(:))) 
        disp('Skipping saving because the map contains all NAN values.')
    else
        hdr.fname = fullfile(respath,'Day2_s.nii');
        if ~exist(hdr.fname,'file')
            spm_write_vol(hdr, out_Day2_s_mean);
        end
    end
    
    % 9: Day12_s
    if all(isnan(out_Day12_s(:)))
        disp('Skipping saving because the map contains all NAN values.')
    else
        hdr.fname = fullfile(respath,'Day12_s.nii');
        if ~exist(hdr.fname,'file')
            spm_write_vol(hdr, out_Day12_s);
        end
    end
    
    % 10: Day12_both
    if all(isnan(out_Day12_both(:)))
        disp('Skipping saving because the map contains all NAN values.')
    else
        hdr.fname = fullfile(respath,'Day12_both.nii');
        if ~exist(hdr.fname,'file')
            spm_write_vol(hdr, out_Day12_both);
        end
    end
    
    % 11: All sessions
    if all(isnan(out_all(:)))
        disp('Skipping saving because the map contains all NAN values.')
    else
        hdr.fname = fullfile(respath,'All.nii');
        if ~exist(hdr.fname,'file')
            spm_write_vol(hdr, out_all);
        end
    end
    
    % 12: Day 0 only
    if all(isnan(out_Day0(:)))
        disp('Skipping saving because the map contains all NAN values.')
    else
        hdr.fname = fullfile(respath,'Day0.nii');
        if ~exist(hdr.fname,'file')
            spm_write_vol(hdr, out_Day0);
        end
    end
    
    end % end of func conn model
    
end  % end of subjects


%% aggregate contrast maps using one-sample t-tests

gm_mask = fullfile(maskpath,'gm_0.1_2mm.nii');

if agg_contrast_map
    for i = 1:length(FC_models) 
    for c = 1:n_contrast
        
    % where to save aggregate contrast map results
    pathname = fullfile(FCdir, 'Agg_t_test', FC_models{i}, contrast_names{c});
    if ~exist(pathname, 'dir')
        mkdir(pathname);
    end

    % check if t-test has been run in this folder
     if any(size(dir(fullfile(pathname,'spmT_*.nii')),1))
         continue;
     end
    
    % to account for missing data
    subctr = 0;
    clear ffiles
    for subj = 1:nSubIDlist
        SubID = SubIDlist(subj).name;
        filename = sprintf('%s.nii', contrast_names{c});
        filepath = fullfile(FCdir, SubID, 'Contrasts', FC_models{i}, filename);
        if exist(filepath,'file')
            subctr = subctr + 1;
            ffiles{subctr,1} = filepath;
        end
    end
            
    clear matlabbatch
    matlabbatch{1}.spm.stats.factorial_design.dir = {pathname}; 
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = ffiles;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {gm_mask}; % use the gray matter mask
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




%% aggregate contrast maps using multiple regression
% both day 1 and day 2 TMS sessions
% adding order as a covariate

% Initialize arrays for all scans, conditions, subjects, and session orders
all_scans = {};
all_conditions = [];
all_subjects = [];
session_order = [];

use_rest_names = {'S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_use_rest_names = length(use_rest_names);

for i = 1:length(FC_models)  % for each func conn model
modelname = sprintf('GLM_%s',FC_models{i});

for subj = 1:nSubIDlist
    sn = SubIDlist(subj).name;
    
    % coding sham and cTBS as 1 and 2, respectively
    tmp_order = Order(subj);
    switch tmp_order
        case 123
            TMS_types = [2,1,1,2,1,1];            
        case 132
            TMS_types = [2,1,1,1,1,2];       
        case 213
            TMS_types = [1,2,2,1,1,1];            
        case 231
            TMS_types = [1,2,1,1,2,1];            
        case 312
            TMS_types = [1,1,2,1,1,2];               
        case 321
            TMS_types = [1,1,1,2,2,1];
    end
    
    for r = 1:n_use_rest_names % to loop over rest_names
        curr_rest = use_rest_names{r};
        fprintf('%s of %s\n',curr_rest,sn)
        if count_table(subj,r+1)<1 % note here r skips Day0 session
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        modeldir = fullfile(studydir,'FuncConn_PAID_MNI',sn,curr_rest,modelname);
        all_scans = [all_scans; fullfile(modeldir,'beta_0001.nii')];  % use beta map
        all_subjects = [all_subjects; subj];
        all_conditions = [all_conditions; TMS_types(r)];  % Sham (1) or cTBS (2) condition
        session_order = [session_order; r];
    end % end of rest_names
end

% treat sham and cTBS as two regressors
% session order as another covariate
% (ignored subject information)

% where to save aggregate contrast map results
pathname = fullfile(FCdir, 'Agg_multi_reg', FC_models{i});
if ~exist(pathname, 'dir')
    mkdir(pathname);
end

clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = {pathname};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = all_scans;

% Define the conditions (sham and cTBS) as regressors
sham_idx = double(all_conditions == 1);
cTBS_idx = double(all_conditions == 2);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).c = sham_idx;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).cname = 'Sham';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).c = cTBS_idx;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).cname = 'cTBS';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).iCC = 1;

% Include session order as a covariate
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(3).c = session_order;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(3).cname = 'SessionOrder';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(3).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {fullfile(maskpath,'gm_0.1_2mm.nii')}; % use the gray matter mask
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Estimate the second-level model
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(pathname,'SPM.mat')};
 
% Define group-level contrasts
matlabbatch{3}.spm.stats.con.spmmat = {fullfile(pathname,'SPM.mat')};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Sham > cTBS';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1 0];  % Include zero for covariate
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'cTBS > Sham';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1 0];  % Include zero for covariate
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

% save job
jobname = fullfile(pathname, 'Job_agg.mat');
save(jobname, 'matlabbatch')

% Run the job
spm_jobman('run', matlabbatch);
clear matlabbatch

end
