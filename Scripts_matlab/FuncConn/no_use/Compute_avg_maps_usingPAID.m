
% this code tries to calculate the average func conn map across TMS
% sessions and subjects, to have a rough idea of how the conn map changes
% with TMS types


clear; clc; close all

studydir = '/Volumes/X9Pro/NODEAP/MRI';
SubIDlist = dir(fullfile(studydir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis');
SubInfo = readtable(fullfile(HomeDir,'ProcessedData/SubConds.xlsx'));

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

clear regions
% Right anterior OFC
regions(1).seedMNI = [34, 54, - 14];
regions(1).seedName = 'aOFC_seed_right';
regions(1).targetMNI = [44, 28, 38];
regions(1).targetName = 'aOFC_target_right';

% Right posterior OFC
regions(2).seedMNI = [28, 38, -16];
regions(2).seedName = 'pOFC_seed_right';
regions(2).targetMNI = [46, 38, 14];
regions(2).targetName = 'pOFC_target_right';

convol = cell(nSubIDlist,n_rest_names);
% Initialize cell array TMS_types with the dimensions 7x48
TMS_types = cell(nSubIDlist,n_rest_names);

%%
for subj = 1:nSubIDlist

    sn = SubIDlist(subj).name;
    subdir = fullfile(studydir,sn);
    unnormdir = fullfile(subdir, 'UnnormedMasks');  %path to native space masks for each subejct
   
for r = 1:n_rest_names % to loop over rest_names
    
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,subdir)
    
    if count_table(subj,r)<1
        disp('Not enough data, skipping this iteration ...')
        continue;
    end
    
    if subj==18 && r==7 % not include this one because it has diff dim 104*104*78
        continue;
    end
    
for i = 1:length(regions) % seed region OFC
    modeldir = fullfile(subdir,'FuncConn_PAID',curr_rest,sprintf('GLM_%s',regions(i).seedName));
    hdr = spm_vol_nifti(fullfile(modeldir,'spmT_0001.nii'));
    convol{subj,r} = spm_read_vols(hdr);   
    
end
end
    tmp_order = SubInfo.StimOrder(subj);
    switch tmp_order
        case 123
            TMS_types(subj,:) = {'S', 'C', 'S', 'S', 'C', 'S', 'S'};
        case 132
            TMS_types(subj,:) = {'S', 'C', 'S', 'S', 'S', 'S', 'C'};
        case 213
            TMS_types(subj,:) = {'S', 'S', 'C', 'C', 'S', 'S', 'S'};
        case 231
            TMS_types(subj,:) = {'S', 'S', 'C', 'S', 'S', 'C', 'S'};
        case 312
            TMS_types(subj,:) = {'S', 'S', 'S', 'C', 'S', 'S', 'C'};
        case 321
            TMS_types(subj,:) = {'S', 'S', 'S', 'S', 'C', 'C', 'S'};
    end
end

%%
% Convert the TMS_types cell array to a logical array for 'C' and 'S'
isC = cellfun(@(x) strcmp(x, 'C'), TMS_types);
isS = cellfun(@(x) strcmp(x, 'S'), TMS_types);

% Extract cell arrays of 'C' and 'S' matrices
convolC = convol(isC);
convolS = convol(isS);

% Remove empty cells
convolC = convolC(~cellfun('isempty', convolC));
convolS = convolS(~cellfun('isempty', convolS));

% Convert cell arrays to 4D arrays
convolC_4D = cat(4, convolC{:});
convolS_4D = cat(4, convolS{:});

% Calculate the sum of each 4D array along the 4th dimension
sumC = sum(convolC_4D, 4);
sumS = sum(convolS_4D, 4);

% Calculate the number of 'C' and 'S' elements
countC = numel(convolC);
countS = numel(convolS);

% Calculate the average matrices
avgC = sumC / countC;
avgS = sumS / countS;

hdr.fname = fullfile(HomeDir,'ProcessedData/AvgFuncConnMap_C.nii');
spm_write_vol(hdr, avgC);
        
hdr.fname = fullfile(HomeDir,'ProcessedData/AvgFuncConnMap_S.nii');
spm_write_vol(hdr, avgS);

hdr.fname = fullfile(HomeDir,'ProcessedData/AvgFuncConnMap_CminusS.nii');
spm_write_vol(hdr, avgC - avgS);

        




