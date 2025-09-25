%% ========================================================================
%  Script: count_functional_sessions.m
%
%  Purpose
%    Count, per subject and per resting session, whether combined
%    multi-echo functional images exist with the expected number of volumes.
%
%  Definition of counts (per cell):
%    1   = full-length run present (nvol == 310)
%    0.5 = short run present (nvol == 205)
%    0   = missing / not complete
%
%  Inputs (directory layout)
%    <dat_folder>/<SubID>/nifti/functional/<Sess>/fvol_###.nii
%
%  Output
%    /Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx  (table with RowNames=SubIDs)
% ========================================================================

%% Setup ------------------------------------------------------------------
clear; clc;

dat_folder = '/Volumes/X9Pro/NODEAP/MRI';

% Subjects
SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);
nSub      = numel(SubIDlist);

% Sessions (columns)
rest_names   = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = numel(rest_names);

% Expected volumes
FULL_NSCANS  = 310;
SHORT_NSCANS = 205;

% Results matrix
count_table = zeros(nSub, n_rest_names);

%% Iterate subjects/sessions ----------------------------------------------
for subj = 1:nSub
    SubID  = SubIDlist(subj).name;
    SubDir = fullfile(dat_folder, SubID);
    niidir = fullfile(SubDir, 'nifti');

    for r = 1:n_rest_names
        curr_rest = rest_names{r};
        functdir  = fullfile(niidir, 'functional', curr_rest);  % combined outputs

        fprintf('%s — %s\n', SubID, curr_rest);

        % Skip if session folder is missing
        if ~exist(functdir, 'dir')
            continue;
        end

        % Count combined volumes present
        nvol = numel(dir(fullfile(functdir, 'fvol_*.nii')));
        if nvol == FULL_NSCANS
            count_table(subj, r) = 1;
        elseif nvol == SHORT_NSCANS
            count_table(subj, r) = 0.5;
        else
            % leave as 0 (missing/incomplete)
        end
    end
end

%% Save table --------------------------------------------------------------
row_names = {SubIDlist.name}.';  % column cellstr
myTable   = array2table(count_table, 'VariableNames', rest_names, 'RowNames', row_names);

out_xlsx  = '/Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx';
writetable(myTable, out_xlsx, 'WriteRowNames', true);
fprintf('Saved count table to: %s\n', out_xlsx);
