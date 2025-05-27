
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

% load aal atlas
aal_atlas_3D = spm_read_vols(spm_vol('./aal116/raal.nii'));  % size = [X Y Z]
aal_labels_gm = aal_atlas_3D(gm_idx);  % label of each voxel in dat_gm

roi_labels = unique(aal_labels_gm);
roi_labels(roi_labels == 0) = [];  % remove background
nROIs = length(roi_labels); % 116 ROIs

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
        
        save_dir = fullfile(studydir, 'FuncConn_AAL', SubID, curr_rest);
        save_name = fullfile(save_dir,'conn_matrix.mat');

        if exist(save_name, 'file')
            continue;
        end
     
        data = load(fullfile(save_dir,'tc_filtered_2mm.mat'));  
        dat_filtered = data.dat_filtered;
     
        ROI_signals = zeros(nscans, nROIs);

        for i = 1:nROIs
            label = roi_labels(i);
            roi_voxel_idx = find(aal_labels_gm == label);  % column indices in dat_gm
            if ~isempty(roi_voxel_idx)
                ROI_signals(:, i) = mean(dat_filtered(:, roi_voxel_idx), 2);
            end
        end

        FC_matrix = corr(ROI_signals);  % final FC matrix

      
        % Save
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end

        s = struct("FC_matrix", FC_matrix);
        save(save_name, "-fromstruct", s);  

    catch ME
        fprintf('❌ Error with %s %s: %s\n', SubID, curr_rest, ME.message);
    end
end

