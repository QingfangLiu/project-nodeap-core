
%% This script tries to decode rs-fMRI sessions based on patterns of ALFF values

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
Order = table_data.StimOrder;

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

MRIcount = readtable(fullfile(studydir,'MRI_func_count.xlsx'),"ReadRowNames",true);
maskpath = fullfile(HomeDir,'Scripts_matlab/GlobalConn/mask');

addpath('/Users/liuq13/libsvm-3.25/matlab') % add package for svm

Option = 'fALFF';    % to do this with ALFF or fALFF or zALFF

%% read ROI idx 
roiname = {'aOFC_seed_right','aOFC_target_right','pOFC_seed_right','pOFC_target_right'};

%% extract ALFF values
% do this for each roi
AllALFF = cell(length(roiname),1);

for i = 1:length(roiname) % loop over each ROI
    
    roi_vol = spm_read_vols(spm_vol_nifti(fullfile(maskpath,sprintf('%s.nii',roiname{i}))));
    roi_idx = find(roi_vol);
    roi_size = length(roi_idx);
    ALFF = nan(roi_size,n_rest_names,nSubIDlist);
    mat_TMS_types = nan(n_rest_names,nSubIDlist);
    
for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;   
    fprintf('%s of %s\n',SubID,roiname{i})
    
    % coding sham and cTBS as 1 and 2, respectively
    % null as 3
    tmp_order = Order(subj);
    switch tmp_order
        case 123
            TMS_types = [3,2,1,1,2,1,1];            
        case 132
            TMS_types = [3,2,1,1,1,1,2];       
        case 213
            TMS_types = [3,1,2,2,1,1,1];            
        case 231
            TMS_types = [3,1,2,1,1,2,1];            
        case 312
            TMS_types = [3,1,1,2,1,1,2];               
        case 321
            TMS_types = [3,1,1,1,2,2,1];
    end
    mat_TMS_types(:,subj) = TMS_types;
    
    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
       
        if MRIcount{SubID,curr_rest}<1  % only with complete time-series
            disp('Not enough data, skipping this iteration ...')
            mat_TMS_types(r,subj) = nan; % also update the mat containing TMS types if missing data
            continue;
        end
        ALFFvalues = spm_read_vols(spm_vol_nifti(fullfile(studydir, 'ALFF_PAID', SubID, sprintf('%s_%s.nii', curr_rest, Option)))); % in 3mm space, no smoothing
        ALFF(:,r,subj) = ALFFvalues(roi_idx);    
    end
end

    % put data from ROIs together
    AllALFF{i,1} = ALFF;

end

%%
% two things to decode: sham vs cTBS, aOFC vs pOFC
% decoding is done via leave-one-subject-out CV
% this way it keeps individual sessions
accvals = zeros(nSubIDlist,length(roiname)); % pre-define before saving acc vals in each LOOCV 
SubVec = 1:nSubIDlist;

for i = 1:length(roiname) % loop over each ROI
    
    for subj = 1:nSubIDlist
        use_data = AllALFF{i};
        use_data = squeeze(use_data(:,:,subj));
        use_data = use_data';

        use_labels = mat_TMS_types(:,subj);
        to_include = use_labels == 1 | use_labels ==2;

        Labels{subj,1} = use_labels(to_include);
        data{subj,1} = use_data(to_include);

    end
   
    for subj = 1:nSubIDlist

        subj_train = SubVec(SubVec~=subj); % leave one subj cross-validation
        subj_test = subj;

        vectors_train = [];
        labels_train = [];

        for t = subj_train
            vectors_train = [vectors_train; data{t,1}];
            labels_train = [labels_train; Labels{t,1}];
        end

        vectors_test = data{subj_test,1};
        labels_test = Labels{subj_test,1};

        model = svmtrain(labels_train,vectors_train,'-s 0 -t 0 -q -c 0.0001'); % -s svm_type; -t kernel_type; -q quite mode
        [~, accuracy, ~] = svmpredict(labels_test,vectors_test,model,'-q'); % pred_label: labels of testing data predicted from the trained model
        accvals(subj,i) = accuracy(1);
    end   
end

%% average across sessions before decoding
% so sham and cTBS are fully balanced
accvals = zeros(nSubIDlist,length(roiname)); % pre-define before saving acc vals in each LOOCV 
SubVec = 1:nSubIDlist;

for i = 1:length(roiname) % loop over each ROI
    
    % organize data and Labels first
    for subj = 1:nSubIDlist
        use_data = AllALFF{i};
        use_data = squeeze(use_data(:,:,subj));
        use_data = use_data';

        use_labels = mat_TMS_types(:,subj);
        sham_idx = use_labels == 1;
        cTBS_idx = use_labels == 2;
        data_sham = mean(use_data(sham_idx,:),1);
        data_cTBS = mean(use_data(cTBS_idx,:),1);

        Labels{subj,1} = [1;2];
        data{subj,1} = [data_sham;data_cTBS];
    end
   
    for subj = 1:nSubIDlist
        subj_train = SubVec(SubVec~=subj); % leave one subj cross-validation
        subj_test = subj;
        vectors_train = [];
        labels_train = [];
        for t = subj_train
            vectors_train = [vectors_train; data{t,1}];
            labels_train = [labels_train; Labels{t,1}];
        end
        vectors_test = data{subj_test,1};
        labels_test = Labels{subj_test,1};
        model = svmtrain(labels_train,vectors_train,'-s 0 -t 0 -q -c 0.0001'); % -s svm_type; -t kernel_type; -q quite mode
        [~, accuracy, ~] = svmpredict(labels_test,vectors_test,model,'-q'); % pred_label: labels of testing data predicted from the trained model
        accvals(subj,i) = accuracy(1);
    end
    
end



                