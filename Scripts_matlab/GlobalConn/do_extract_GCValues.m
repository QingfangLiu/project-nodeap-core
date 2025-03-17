
%% This script extracts GC values in each seed and target regions

clc; clear;
userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis');
maskpath = fullfile(HomeDir,'Scripts_matlab/GlobalConn/mask');

studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

file_path = fullfile(HomeDir,'ProcessedData','SubConds.xlsx');
table_data = readtable(file_path); 

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

Option = 'Unsigned'; % Unsigned or Signed
MRIcount = readtable(fullfile(studydir,'MRI_func_count.xlsx'),"ReadRowNames",true);


%% make ROIs on 3mm space

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


%% make spherical ROIs surrounding the coordinates, constrained by the tissue probability map
% in 3mm space

for i = 1:length(regions)
    seedmask{i} = fullfile(maskpath, sprintf('%s.nii',regions(i).seedName));
    targetmask{i} = fullfile(maskpath, sprintf('%s.nii',regions(i).targetName));
end

%check if spherical ROIs in MNI space have already been made
if ~exist(seedmask{2},'file') && ~exist(targetmask{2},'file')
    
    rget = 3; % updated from 4 (*2mm) to 3 (*3mm)
    
    % load gray matter mask in 3mm space
    % ROIs are constrained by the gray matter mask
    tmpl = fullfile(maskpath,'gm_0.1_3mm.nii');
    tpmhdr = spm_vol_nifti(tmpl);
    [tpmvol,XYZ] = spm_read_vols(tpmhdr);
    tpmidx = find(tpmvol>0);
    sz = size(tpmvol);
    [MX, MY, MZ] = ndgrid(1:sz(1),1:sz(2),1:sz(3));
    ref_vox = round(sz/2);
    radii = sqrt((MX-ref_vox(1)).^2 + (MY-ref_vox(2)).^2 + (MZ-ref_vox(3)).^2);
    radius_index = find(radii<rget) - sub2ind(sz,ref_vox(1),ref_vox(2),ref_vox(3));
    
    for i = 1:length(regions)
        
        % seed regions
        % coordidx = find(ismember(XYZ',regions(i).seedMNI,'rows')>0);
        [~, coordidx] = min(mean(abs(XYZ' - regions(i).seedMNI),2)); % find the idx closest to the seed MNI
        idx = intersect(radius_index + coordidx, tpmidx);
        newhdr = tpmhdr;
        newvol = zeros(sz);
        newhdr.fname =  seedmask{i};
        newvol(idx) = 1;
        spm_write_vol(newhdr,newvol);
        
        % target regions
        % coordidx = find(ismember(XYZ',regions(i).targetMNI,'rows')>0);
        [~, coordidx] = min(mean(abs(XYZ' - regions(i).targetMNI),2)); % find the idx closest to the target MNI
        idx = intersect(radius_index + coordidx, tpmidx);
        newhdr = tpmhdr;
        newvol = zeros(sz);
        newhdr.fname = targetmask{i};
        newvol(idx) = 1;
        spm_write_vol(newhdr,newvol);
    end
else
    fprintf('\nSeed and Target masks already made\n');
end

%% read ROI idx 
roiname = {'aOFC_seed_right','aOFC_target_right','pOFC_seed_right','pOFC_target_right'};
clear roiidx
for i = 1:length(roiname)
    roi_vol = spm_read_vols(spm_vol_nifti(fullfile(maskpath,sprintf('%s.nii',roiname{i}))));
    roiidx{i} = find(roi_vol);
end

%% extract global conn values

MeanGC = nan(n_rest_names,length(roiidx),nSubIDlist);

for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;   
    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
        
        fprintf('%s of %s\n',curr_rest,SubID)
    
        if MRIcount{SubID,curr_rest}==0
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
    
        for i = 1:length(roiidx) % loop over each ROI
            GCvalues = spm_read_vols(spm_vol_nifti(fullfile(studydir, 'GlobalConn_PAID', SubID, Option, sprintf('%s.nii', curr_rest)))); % in 3mm space, no smoothing
            MeanGC(r,i,subj) = nanmean(GCvalues(roiidx{i}));    
        end
    end
end

% save this for later process
path = fullfile(HomeDir,'ProcessedData',sprintf('GC_subs_sessions_%s.mat',Option));
save(path, 'MeanGC','roiname'); 

                