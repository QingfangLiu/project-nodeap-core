

% This is to count which subject and which scanning session has MRI data
% (with nscans=310 volumes)
% based on combined echoes data

clear; clc;
HomeDir = '/Volumes/X9Pro/NODEAP/MRI';

nscans = 310; % 310 volumes for ME sequnce

SubIDlist = dir(fullfile(HomeDir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);
count_table = zeros(nSubIDlist,n_rest_names);

%%

for subj = 1:nSubIDlist

SubID = SubIDlist(subj).name;
SubDir = fullfile(HomeDir,SubID);
niidir = fullfile(SubDir, 'nifti');

for r = 1:n_rest_names % to loop over rest_names
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubID)
    functdir = fullfile(niidir, 'functional', curr_rest); % path to rsfRMI images
    
    % count number of 'fvol_*.nii' volumes 
    nvol = length(dir(fullfile(functdir,'fvol_*.nii')));
    if ~exist(functdir,'dir') || nvol==0
        continue;
    end
    
    if nvol==nscans % if the functional dir has all combined scans
        count_table(subj,r) = 1;
    end 
    if nvol==205
        count_table(subj,r) = 0.5;
    end 
end
end

%%
myTable = array2table(count_table, 'VariableNames', rest_names, 'RowNames', {SubIDlist.name});
writetable(myTable, '/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx','WriteRowNames', true);

