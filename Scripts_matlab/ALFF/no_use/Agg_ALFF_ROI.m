
% aggregate the ROI-based alff and falff results

clear; clc; close all

studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

%%
nroi = 4;
MeanALFF = nan(n_rest_names,nroi,nSubIDlist);
MeanfALFF = nan(n_rest_names,nroi,nSubIDlist);

for subj = 1:nSubIDlist

    sn = SubIDlist(subj).name;
    respath = fullfile(studydir, 'ALFF_PAID', sn);
    
for r = 1:n_rest_names % to loop over rest_names
    
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,sn)
    
    if count_table(subj,r)<1
        disp('Not enough data, skipping this iteration ...')
        continue;
    end

    mat_name = fullfile(respath,sprintf('ALFF_ROI_%s.mat',curr_rest));
    load(mat_name)
    
    for k = 1:4
        tmp_idx = ALFF_ROI_vals(:,3)==k;
        MeanALFF(r,k,subj) = nanmean(ALFF_ROI_vals(tmp_idx,4));
        MeanfALFF(r,k,subj) = nanmean(ALFF_ROI_vals(tmp_idx,5));
    end
    
end
end

%%
path = '../../ProcessedData/ALFF_subs_sessions.mat';
save(path, 'MeanALFF','MeanfALFF'); 



