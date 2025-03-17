
% This analysis compares the functional connectivity of the seed region
% with the aOFC target and the pOFC target. Theoretically we expect that
% aOFC seed region is more connected to the aOFC target region, while pOFC
% seed region is more connected to the pOFC target region. But there is no
% guarantee it is the truth. This analysis leverages the GLM used to get
% the TMS coordinates and then extract t-values in both target regions.

% This is using data processed after the PAID method and spatial smoothing
% across all TMS sessions

clear; clc; close all

studydir = '/Volumes/X9Pro/NODEAP/MRI';
SubIDlist = dir(fullfile(studydir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

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

% plot distribution for each subject?
plot_dist = 0;

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);
mean_convals = nan(length(regions),2,nSubIDlist,n_rest_names);

%%
for subj = 1:nSubIDlist

    sn = SubIDlist(subj).name;
    subdir = fullfile(studydir,sn);
    unnormdir = fullfile(subdir, 'UnnormedMasks');  %path to native space masks for each subejct
    TMSpostdir = fullfile(studydir,'TMS_post_hoc');
    if ~exist(TMSpostdir,'dir')
        mkdir(TMSpostdir)
    end
    
for r = 1:n_rest_names % to loop over rest_names
    
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,subdir)
    
    if count_table(subj,r)<1
        disp('Not enough data, skipping this iteration ...')
        continue;
    end

for i = 1:length(regions) % seed region OFC
    modeldir = fullfile(subdir,'FuncConn_PAID',curr_rest,sprintf('GLM_%s',regions(i).seedName));
    % read connectivitiy t values in each of the target regions
    for t = 1:2 % target region LPFC
        [targetvol, targetXYZ] = spm_read_vols(spm_vol_nifti(fullfile(unnormdir,sprintf('rw%s.nii',regions(t).targetName))));
        targetidx = find(targetvol==1);
        convol = spm_read_vols(spm_vol_nifti(fullfile(modeldir,'spmT_0001.nii')));
        vals = convol(targetidx);
        mean_convals(i,t,subj,r) = mean(vals);
    end
end
end
end

%%
path = '../../ProcessedData/SeedStimConn_subs_sessions.mat';
save(path, 'mean_convals'); 
