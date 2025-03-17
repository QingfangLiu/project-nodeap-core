

% use 3mm MNI space
% this one calculates map with all gm voxels for ALFF and fALFF values
% focus on 4 ROIs

clear; clc; close all
studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis');

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

table_data = readtable(fullfile(HomeDir,'ProcessedData/SubConds.xlsx'));
Order = table_data.StimOrder;

nscans = 310; % 310 volumes for ME sequnce
nechos = 3;
TR = 1.5;
TEs = [14.60, 39.04, 63.48];
rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

%%
maskpath = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis/Scripts_matlab/GlobalConn/mask');
masks = {'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'
    };  % masks in 3mm space
n_masks = length(masks);

labels = {'aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'};

gm_nii = fullfile(maskpath,'gm_0.1_3mm.nii'); % gray matter mask
hdr = spm_vol(gm_nii);
gm_dat = spm_read_vols(hdr);
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

array_size = size(gm_dat);  % 53*63*52
[check_exist, roi_gm_linear_indices] = ismember(voxel_idx, gm_idx);
all(check_exist)
% check_exist: make sure all ROI indices are within the gray matter (as it
% should)
% roi_gm_linear_indices: indices of ROIs that are in the gray matter

all_idx = [voxel_idx,roi_gm_linear_indices,roi_idx];
nvox = size(all_idx,1);

% col1: voxel in 3mm MNI space
% col2: voxel in gray matter linear index
% col3: which roi: 1,2,3,4

%%
for subj = 1:nSubIDlist

    sn = SubIDlist(subj).name;
    
    % prepare where to save the ALFF Results
    respath = fullfile(studydir, 'ALFF_PAID', sn);
    if ~exist(respath,'dir')
        mkdir(respath)
    end
    
    for r = 1:n_rest_names % to loop over rest_names
        
        tic
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,sn)
        
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        
        mat_name = fullfile(respath,sprintf('ALFF_ROI_%s.mat',curr_rest));
%         if exist(mat_name,'file') 
%             continue;  
%         end
        
        load(fullfile(studydir, 'GlobalConn_PAID', sn, curr_rest, 'tc_filtered_3mm.mat')); 
        
        sampleInterval = TR; 
        lowFreq = 0.01; % Lower bound of frequency range in Hz
        highFreq = 0.08; % Upper bound of frequency range in Hz

        ALFF = zeros(nvox,1);          % Initialize the ALFF array
        fALFF = zeros(nvox,1);         % Initialize the fALFF array
        pool = gcp();                       % This will start a parallel pool if one is not already started
        parfor i = 1:nvox
            use_idx = all_idx(i,2);        % index in the gray matter
            [ALFF(i),fALFF(i)] = computePowerMeasures(dat(:, use_idx), sampleInterval, lowFreq, highFreq);
            if ~mod(i,50)
                fprintf('\t\t%d of %d\n', i, nvox);
            end
        end
        
        ALFF_ROI_vals = [all_idx,ALFF,fALFF];
        save(mat_name,'ALFF_ROI_vals')
        toc
            
    end
end




