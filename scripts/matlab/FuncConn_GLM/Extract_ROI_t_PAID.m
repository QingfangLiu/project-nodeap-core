
% This analysis compares the functional connectivity of the seed region
% with the aOFC target and the pOFC target. Theoretically we expect that
% aOFC seed region is more connected to the aOFC target region, while pOFC
% seed region is more connected to the pOFC target region. But there is no
% guarantee it is the truth. This analysis leverages the GLM used to get
% the TMS coordinates and then extract t-values in both target regions.

% This is using data processed after the PAID method and spatial smoothing
% across all TMS sessions

clear; clc; close all

studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

FCdir = fullfile(studydir,'FuncConn_PAID_MNI');

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

% Right anterior OFC
regions(1).seedName = 'aOFC_seed';
regions(1).targetName = 'aOFC_target_right';

% Right posterior OFC
regions(2).seedName = 'pOFC_seed';
regions(2).targetName = 'pOFC_target_right';

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);
mean_convals = nan(length(regions),2,nSubIDlist,n_rest_names);

%%
for subj = 1:nSubIDlist

    sn = SubIDlist(subj).name;
    
for r = 1:n_rest_names % to loop over rest_names
    
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,sn)
    
    if count_table(subj,r)<1
        disp('Not enough data, skipping this iteration ...')
        continue;
    end

for i = 1:2 % seed region OFC
    modeldir = fullfile(FCdir,sn,curr_rest,sprintf('GLM_%s',regions(i).seedName));
    % read connectivitiy t values in each of the target regions
    for t = 1:2 % target region LPFC
        [targetvol, targetXYZ] = spm_read_vols(spm_vol_nifti(fullfile(studydir,'ConnectivityMasks',sprintf('%s.nii',regions(t).targetName))));
        targetidx = find(targetvol==1);
        convol = spm_read_vols(spm_vol_nifti(fullfile(modeldir,'spmT_0001.nii')));
        mean_convals(i,t,subj,r) = mean(convol(targetidx));
    end
end
end
end

%%
path = '../../ProcessedData/SeedStimConn_subs_sessions.mat';
save(path, 'mean_convals'); 
