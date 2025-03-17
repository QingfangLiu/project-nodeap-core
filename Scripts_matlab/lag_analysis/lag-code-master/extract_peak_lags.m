
% this script extracts values from peak lag estimates
% and plot individual results together

clc; clear; close all
studydir = '/Volumes/X9Pro/NODEAP';
userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Library/CloudStorage/Box-Box/NODEAP_data_analysis');

maskpath = fullfile(userHome,'Library/CloudStorage/Box-Box/NODEAP_data_analysis/Scripts_matlab/GlobalConn/mask');

MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

masks = {'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'
    };  % masks in 3mm space
n_masks = length(masks);

gm_nii = fullfile(maskpath,'gm_0.1_3mm.nii'); % gray matter mask
gm_dat = spm_read_vols(spm_vol(gm_nii));
gm_idx = find(gm_dat > 0);
n_gm_vox = length(gm_idx);

%% concatenante the ROIs
voxel_idx = [];
roi_idx = [];
for i = 1:n_masks
    [targetvol, ~] = spm_read_vols(spm_vol_nifti(fullfile(maskpath,masks{i})));
    tmploc = find(targetvol==1);  % loc idx in 53*63*52 space
    voxel_idx = [voxel_idx;tmploc];
    roi_idx = [roi_idx;repmat(i,length(tmploc),1)];
end

%%
array_size = size(gm_dat);  % 53*63*52
[check_exist, roi_gm_linear_indices] = ismember(voxel_idx, gm_idx);
all_idx = [voxel_idx,roi_gm_linear_indices,roi_idx];
nvox = size(all_idx,1);  % total number of voxels within the ROIs

% col1: voxel in 3mm MNI space
% col2: voxel in gray matter linear index
% col3: which roi: 1,2,3,4


%% extract peak lag values
% from 4 sections: mask2-mask1; mask2-mask3; mask4-mask1; mask4-mask3
% see result slides

sections = {'aOFC-target-aOFC-seed','aOFC-target-pOFC-seed',...
    'pOFC-target-aOFC-seed','pOFC-target-pOFC-seed'};
mask1_idx = find(all_idx(:,3) == 1);
mask2_idx = find(all_idx(:,3) == 2);
mask3_idx = find(all_idx(:,3) == 3);
mask4_idx = find(all_idx(:,3) == 4);

section_coverage = cell(4,2);
section_coverage{1,1} = mask2_idx;
section_coverage{1,2} = mask1_idx;
section_coverage{2,1} = mask2_idx;
section_coverage{2,2} = mask3_idx;
section_coverage{3,1} = mask4_idx;
section_coverage{3,2} = mask1_idx;
section_coverage{4,1} = mask4_idx;
section_coverage{4,2} = mask3_idx;

%%

Mean_peak_lags = nan(n_rest_names,length(sections),nSubIDlist);

for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name; 
    disp(['Processing ' SubID]);
    
    lag_out_dir = fullfile(studydir, 'lag_out_PAID', SubID); 
    
    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,SubID)
    
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
    
        for i = 1:length(sections) % loop over each section
            
            mat_name = fullfile(lag_out_dir,sprintf('lag_corr_%s.mat',curr_rest));
            load(mat_name)
            use_peak_lags = peak_lags(section_coverage{i,1},section_coverage{i,2});
            Mean_peak_lags(r,i,subj) = nanmean(use_peak_lags,'all');
        
        end
    end
end

%% save ALFF values for later process
path = fullfile(HomeDir,'ProcessedData','peak_lags_subs_sessions.mat');
save(path, 'Mean_peak_lags','sections'); 

