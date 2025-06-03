
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

% load aal atlas
aal_atlas_3D = spm_read_vols(spm_vol('./aal116/raal.nii'));  % size = [X Y Z]
aal_labels_gm = aal_atlas_3D(gm_idx);  % label of each voxel in gm

roi_labels = unique(aal_labels_gm); % unique labels (should be 0-116)
roi_labels(roi_labels == 0) = [];  % remove background 
nROIs = length(roi_labels); % 116 ROIs

nscans = 310;

%%
% Build job list
taskList = {};
count = 0;
for subj = 1:length(SubIDlist)
    SubID = SubIDlist(subj).name;
    for r = 1:length(rest_names)
        curr_rest = rest_names{r};
        val = MRIcount{SubID, curr_rest};
        if val ~= 1, continue; end  % ignore missing or half-missing sessions
        
        count = count + 1;
        taskList{count,1} = SubID;
        taskList{count,2} = curr_rest;
        taskList{count,3} = nscans;
    end
end

%%

% Start parallel pool if needed
if isempty(gcp('nocreate'))
    parpool('local');
end


%% compute filtered time series

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
        fprintf('Error with %s %s: %s\n', SubID, curr_rest, ME.message);
    end
end


%% compute functional conn using filtered data

parfor t = 1:size(taskList,1)
    SubID = taskList{t,1};
    curr_rest = taskList{t,2};
    nscans = taskList{t,3};

    try
        save_dir = fullfile(studydir, 'FuncConn_AAL', SubID, curr_rest);
        save_name = fullfile(save_dir,'conn_matrix.mat');

        if exist(save_name, 'file')
            continue;
        end
        
        fprintf('Computing FC for %s of %s\n', curr_rest, SubID);
     
        data = load(fullfile(save_dir,'tc_filtered_2mm.mat'));  
        dat_filtered = data.dat_filtered; % [nscans x n_voxels]
        ROI_signals = zeros(nscans, nROIs);

        for i = 1:nROIs
            label = roi_labels(i);
            roi_voxel_idx = find(aal_labels_gm == label);  
            if ~isempty(roi_voxel_idx)
                ROI_signals(:, i) = mean(dat_filtered(:, roi_voxel_idx), 2);
            end
        end

        FC_matrix = corr(ROI_signals);  % final FC matrix
        s = struct("FC_matrix", FC_matrix);
        save(save_name, "-fromstruct", s);  

    catch ME
        fprintf('Error with %s %s: %s\n', SubID, curr_rest, ME.message);
    end
end


%% check for NAN values

problem_sessions = {};

fprintf('Checking NaNs in filtered data and FC matrices...\n');

for subj = 1:length(SubIDlist)
    SubID = SubIDlist(subj).name;

    for r = 1:length(rest_names)
        curr_rest = rest_names{r};

        save_dir = fullfile(studydir, 'FuncConn_AAL', SubID, curr_rest);
        fpath_ts = fullfile(save_dir, 'tc_filtered_2mm.mat');
        fpath_fc = fullfile(save_dir, 'conn_matrix.mat');

        has_issue = false;

        % Check filtered time series
        if exist(fpath_ts, 'file')
            data_ts = load(fpath_ts);
            if any(isnan(data_ts.dat_filtered(:)))
                fprintf('NaN in filtered time series: %s %s\n', SubID, curr_rest);
                has_issue = true;
            end
        else
            fprintf('Missing filtered file: %s %s\n', SubID, curr_rest);
        end

        % Check FC matrix
        if exist(fpath_fc, 'file')
            data_fc = load(fpath_fc);
            if any(isnan(data_fc.FC_matrix(:)))
                fprintf('NaN in FC matrix: %s %s\n', SubID, curr_rest);
                has_issue = true;
            end
        else
            fprintf('Missing FC matrix file: %s %s\n', SubID, curr_rest);
        end

        if has_issue
            problem_sessions{end+1,1} = SubID;
            problem_sessions{end,2} = curr_rest;
        end
    end
end

fprintf('\nSummary of problematic sessions:\n');
disp(cell2table(problem_sessions, 'VariableNames', {'Subject', 'Session'}));


