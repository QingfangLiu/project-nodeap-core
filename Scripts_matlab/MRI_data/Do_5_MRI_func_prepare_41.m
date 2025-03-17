
tic
clear; clc;
HomeDir = '/Volumes/X9Pro/NODEAP/MRI';

SubIDlist = dir(fullfile(HomeDir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

spm('Defaults','fMRI')
global defaults
defaults.stats.maxmem = 2 ^ 34; % based on RAM = 32GB
defaults.stats.resmem = true; % allow temp files to be kept in memory


%%
for subj = 18

SubID = SubIDlist(subj).name;
SubDir = fullfile(HomeDir,SubID);
niidir = fullfile(SubDir, 'nifti');

past_volumes = 0;
for r = 1:6
    if count_table(subj,r)==0
        disp('Not enough data, skipping this iteration ...')
        continue;
    end
    
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubID)
    curr_func_dir = fullfile(niidir,'functional',curr_rest);
    func_4d = fullfile(curr_func_dir, 'fvol_4d.nii');
    
    if ~exist(func_4d,'file')
        disp('Working on 3D to 4D ...')
        func_files_list = dir(fullfile(curr_func_dir,'fvol*.nii'));
        func_files = arrayfun(@(x) fullfile(curr_func_dir, x.name), func_files_list, 'UniformOutput', false);
        spm_file_merge(func_files, func_4d, 0);  % Merge 3D volumes into a 4D volume
        for i = 1:numel(func_files)
            delete(func_files{i});
        end
        func_4d = fullfile(curr_func_dir, 'fvol_4d.nii');
    end
    
    V = spm_vol(func_4d);  % Use SPM to read the NIfTI file
    for i = 1:numel(V)
        mri_data(:,:,:,i + past_volumes) = spm_read_vols(V(i));
    end
    past_volumes = past_volumes + numel(V);
end

    mean_volume_file = V(1);  % Use the header information from the last session    
    mean_volume_file.fname = fullfile(niidir,'functional', 'mean_fvol_1.nii');
    spm_write_vol(mean_volume_file, mean(mri_data,4));
    
    clear mri_data
    r = 7;
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubID)
    curr_func_dir = fullfile(niidir,'functional',curr_rest);
    func_4d = fullfile(curr_func_dir, 'fvol_4d.nii');
    
    if ~exist(func_4d,'file')
        disp('Working on 3D to 4D ...')
        func_files_list = dir(fullfile(curr_func_dir,'fvol*.nii'));
        func_files = arrayfun(@(x) fullfile(curr_func_dir, x.name), func_files_list, 'UniformOutput', false);
        spm_file_merge(func_files, func_4d, 0);  % Merge 3D volumes into a 4D volume
        for i = 1:numel(func_files)
            delete(func_files{i});
        end
        func_4d = fullfile(curr_func_dir, 'fvol_4d.nii');
    end
    
    V = spm_vol(func_4d);  % Use SPM to read the NIfTI file
    for i = 1:numel(V)
        mri_data(:,:,:,i) = spm_read_vols(V(i));
    end
    
    mean_volume_file = V(1);  % Use the header information from the last session    
    mean_volume_file.fname = fullfile(niidir,'functional', 'mean_fvol_2.nii');
    spm_write_vol(mean_volume_file, mean(mri_data,4));

    
end

