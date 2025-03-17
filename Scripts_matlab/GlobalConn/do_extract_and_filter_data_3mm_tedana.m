
clc; clear;

studydir = '/Volumes/X9Pro/NODEAP/MRI';
SubIDlist = dir(fullfile(studydir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

userHome = getenv('HOME');
maskpath = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis/Scripts_matlab/GlobalConn/mask');
MRIcount = readtable('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx',"ReadRowNames",true);

%% load masks & find indexes
% mask should be resliced to match the dimension of scan data (i.e. normalized space)
% created using SPM's tissue probability map

gm_nii = fullfile(maskpath,'gm_0.1_3mm.nii'); % gray matter mask
gm_dat = spm_read_vols(spm_vol(gm_nii));
gm_idx = find(gm_dat > 0);
nvox = length(gm_idx); % number of voxels

wm_nii = fullfile(maskpath,'wm_0.9_3mm.nii'); % white matter mask
wm_dat = spm_read_vols(spm_vol(wm_nii));
wm_idx = find(wm_dat > 0);

csf_nii = fullfile(maskpath,'csf_0.9_3mm.nii'); % CSF mask
csf_dat = spm_read_vols(spm_vol(csf_nii));
csf_idx = find(csf_dat > 0);

%%

for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;
    SubDir = fullfile(studydir,SubID);
    
    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
        
        if MRIcount{SubID,curr_rest}==0
            continue;
        end
        
        if MRIcount{SubID,curr_rest}==0.5
            nscans = 205;
        else
            nscans = 310;
        end
        
        % load data
        path = fullfile (SubDir, 'tedana_output', curr_rest,'aic');
        n = dir(fullfile(path, 's6w3*.nii')); % data with 3mm voxel size
        
        WhereToSave = fullfile(studydir, 'GlobalConn_tedana', SubID, curr_rest);
        save_name = fullfile(WhereToSave,'tc_filtered_3mm.mat');
        if exist(save_name,'file')
            continue;
        end
        
        dat = zeros(nscans,nvox); % initiliaze to get data
        wm_mean = zeros(nscans,1);
        csf_mean = zeros(nscans,1);
        
        fprintf('Extraxt data from %s of %s\n',curr_rest,SubDir)
        
        data_4D = spm_read_vols(spm_vol(fullfile(path, n.name)));
               
        for i = 1:nscans
            tmp = data_4D(:,:,:,i);
            dat(i,:) = tmp(gm_idx); % gray matter
            wm_mean(i) = mean(tmp(wm_idx)); % white matter: take the mean
            csf_mean(i) = mean(tmp(csf_idx)); % CSF: take the mean
        end
                
        % filter data
        mreg = load(fullfile(studydir, 'NuisanceRegressor',SubID, sprintf('nuisance_regressors_%s.txt',curr_rest)));
        mreg = [zscore([mreg, mean(dat,2), wm_mean, csf_mean, [1:nscans]']), ones(nscans,1)]; % add mean(gm), mean(wm), mean(csf), drift, and constant
        b = inv(mreg'*mreg)*(mreg'*dat);
        dat = dat - mreg*b;
        
        % save data
        if ~exist(WhereToSave,'dir')
            mkdir(WhereToSave)
        end
        save(save_name, 'dat');       
       
    end % end of session loop
end

