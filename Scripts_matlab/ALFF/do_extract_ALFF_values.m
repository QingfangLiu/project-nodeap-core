
%% This script extracts ALFF values in each seed and target regions

clc; clear;
userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis');

studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

file_path = fullfile(HomeDir,'ProcessedData','SubConds.xlsx');
table_data = readtable(file_path); 

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

MRIcount = readtable(fullfile(studydir,'MRI_func_count.xlsx'),"ReadRowNames",true);

%% read ROI idx 
roiname = {'aOFC_seed_right','aOFC_target_right','pOFC_seed_right','pOFC_target_right'};
clear roiidx
maskpath = fullfile(HomeDir,'Scripts_matlab/GlobalConn/mask');

for i = 1:length(roiname)
    roi_vol = spm_read_vols(spm_vol_nifti(fullfile(maskpath,sprintf('%s.nii',roiname{i}))));
    roiidx{i} = find(roi_vol);
end

Option = 'zALFF';    % to do this with ALFF or fALFF or zALFF

%% extract ALFF values
MeanALFF = nan(n_rest_names,length(roiidx),nSubIDlist);

for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;   
    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
        
        fprintf('%s of %s\n',curr_rest,SubID)
    
        if MRIcount{SubID,curr_rest}<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
    
        for i = 1:length(roiidx) % loop over each ROI
            ALFFvalues = spm_read_vols(spm_vol_nifti(fullfile(studydir, 'ALFF_PAID', SubID, sprintf('%s_%s.nii', curr_rest, Option)))); % in 3mm space, no smoothing
            MeanALFF(r,i,subj) = nanmean(ALFFvalues(roiidx{i}));    
        end
    end
end

%% save ALFF values for later process
path = fullfile(HomeDir,'ProcessedData',sprintf('ALFF_subs_sessions_%s.mat',Option));
save(path, 'MeanALFF','roiname'); 

                