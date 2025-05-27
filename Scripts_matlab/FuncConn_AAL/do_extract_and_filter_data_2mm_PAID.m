
clc; clear;
run('../start.m');

studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

maskpath = fullfile(getenv('HOME'),'NODEAP_scripts','Scripts_matlab/FuncConn_AAL/masks');
MRIcount = readtable('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx',"ReadRowNames",true);

% load masks 
gm_idx  = find(spm_read_vols(spm_vol(fullfile(maskpath, 'gm_0.1_2mm.nii')))  > 0);
wm_idx  = find(spm_read_vols(spm_vol(fullfile(maskpath, 'wm_0.9_2mm.nii'))) > 0);
csf_idx = find(spm_read_vols(spm_vol(fullfile(maskpath, 'csf_0.9_2mm.nii'))) > 0);

%%
% Build job list
taskList = {};
count = 0;
for subj = 1:length(SubIDlist)
    SubID = SubIDlist(subj).name;
    for r = 1:length(rest_names)
        curr_rest = rest_names{r};
        val = MRIcount{SubID, curr_rest};
        if val == 0, continue; end
        nscans = (val == 0.5) * 205 + (val == 1) * 310;

        count = count + 1;
        taskList{count,1} = SubID;
        taskList{count,2} = curr_rest;
        taskList{count,3} = nscans;
    end
end

% Start parallel pool if needed
if isempty(gcp('nocreate'))
    parpool('local');
end

% Parallel processing
parfor t = 1:size(taskList,1)
    SubID = taskList{t,1};
    curr_rest = taskList{t,2};
    nscans = taskList{t,3};

    try
        SubDir = fullfile(MRIdir, SubID);
        path = fullfile(SubDir, 'nifti', 'functional', curr_rest);
        save_dir = fullfile(studydir, 'FuncConn_AAL', SubID, curr_rest);
        save_name = fullfile(save_dir,'tc_filtered_2mm.mat');

        if exist(save_name, 'file')
            continue;
        end

        n = dir(fullfile(path, 's6w2*.nii'));
        if isempty(n), continue; end

        fprintf('Extracting data from %s of %s\n', curr_rest, SubID);

        data_4D = spm_read_vols(spm_vol(fullfile(path, n(1).name)));
        if size(data_4D, 4) ~= nscans
            warning('Mismatch in %s %s: expected %d scans, found %d', SubID, curr_rest, nscans, size(data_4D, 4));
            continue;
        end

        dat = reshape(data_4D, [], nscans)';  % [nscans x n_voxels]
        dat_gm = dat(:, gm_idx);
        wm_mean = mean(dat(:, wm_idx), 2);
        csf_mean = mean(dat(:, csf_idx), 2);
        gm_mean = mean(dat_gm, 2);
        drift = (1:nscans)';

        mreg_file = fullfile(studydir, 'NRegressor', SubID, sprintf('nuisance_regressors_%s.txt', curr_rest));
        if ~exist(mreg_file, 'file')
            warning('Missing regressor file for %s %s', SubID, curr_rest);
            continue;
        end
        mreg = load(mreg_file);
        mreg = [mreg, gm_mean, wm_mean, csf_mean, drift];
        mreg = zscore(mreg);
        mreg = [mreg, ones(nscans,1)];

        % Regress out nuisance
        b = mreg \ dat_gm;
        dat_filtered = dat_gm - mreg * b;

        % Save
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end

        s = struct("dat_filtered", dat_filtered);
        save(save_name, "-fromstruct", s);  

    catch ME
        fprintf('❌ Error with %s %s: %s\n', SubID, curr_rest, ME.message);
    end
end

