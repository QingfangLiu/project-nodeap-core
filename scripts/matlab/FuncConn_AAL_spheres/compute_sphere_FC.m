

% ========================================================================
% Compute FC between four spheres across all subjects and sessions
% ========================================================================

clear; clc; close all

studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
maskpath = fullfile(getenv('HOME'),'NODEAP_scripts','Scripts_matlab/FuncConn_AAL_spheres/masks');
MRIcount = readtable('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx',"ReadRowNames",true);

gm_idx  = find(spm_read_vols(spm_vol(fullfile(maskpath, 'gm_0.1_2mm.nii')))  > 0);
sphere_path = '/Volumes/X9Pro/NODEAP/ConnectivityMasks';
spheres = {'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'
    };  
n_spheres = length(spheres);

% ---------------------
% Get voxel indices for each sphere
% ---------------------
LOC_sphere = cell(n_spheres,1); 
for i = 1:n_spheres
    targetvol = spm_read_vols(spm_vol_nifti(fullfile(sphere_path,spheres{i})));
    tmploc = find(targetvol==1); 
    [~, LOC_sphere{i}] = ismember(tmploc,gm_idx);
end

% ---------------------
% Build task list of valid subject-session pairs
% ---------------------

taskList = {};
count = 0;
for subj = 1:length(SubIDlist)
    SubID = SubIDlist(subj).name;
    for r = 1:length(rest_names)
        curr_rest = rest_names{r};
        if MRIcount{SubID, curr_rest} ~= 1, continue; end  % ignore missing or half-missing sessions
        count = count + 1;
        taskList(count,:) = {SubID, curr_rest};
    end
end

% ---------------------
% Main loop: Compute sphere-to-sphere FC
% ---------------------

for t = 1:size(taskList,1)
    SubID = taskList{t,1};
    curr_rest = taskList{t,2};
    
    try
        save_dir = fullfile(studydir, 'FuncConn_spheres', SubID, curr_rest);
        save_path = fullfile(save_dir, 'conn_matrix_4spheres.mat');
        if exist(save_path,'file'), continue; end
        
        fprintf('Computing sphere-sphere FC for %s of %s\n', curr_rest, SubID);
        
        dat_dir = fullfile(studydir, 'FuncConn_AAL', SubID, curr_rest);
        dat_path = fullfile(dat_dir,'tc_filtered_2mm.mat');
        
        if ~exist(dat_path, 'file')
            warning('Missing data file: %s', dat_path);
            continue;
        end
        dat_struct = load(dat_path);
        dat = dat_struct.dat_filtered;
        
        % Extract mean time series for each sphere
        ts_all = zeros(size(dat,1), n_spheres);
        for i = 1:n_spheres
            idx = LOC_sphere{i}(LOC_sphere{i} ~= 0);
            ts_all(:,i) = mean(dat(:,idx), 2);
        end

        % Compute full correlation matrix between spheres
        correlation_matrix = corr(ts_all);

        if ~exist(save_dir, 'dir'), mkdir(save_dir); end
        s = struct("correlation_matrix", correlation_matrix, ...
                   "sphere_order", {spheres});
        save(save_path, "-fromstruct", s);
        
    catch ME
        fprintf('Error with %s %s: %s\n', SubID, curr_rest, ME.message);
    end
end







