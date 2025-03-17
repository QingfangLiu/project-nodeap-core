clear all
root = 'C:\Data\POETMS\';
maskpath = 'C:\Users\kahntt2\Box\scripts\POETMS\masks';
addpath('C:\spm12\')

stim = [2, 4, 5, 6, 7, 12, 13, 16, 21, 24, 27, 33, 38, 40, 41, 44, 46, 61, 63, 65, 66, 72, 74, 77, 78, 79];
sham = [11, 15, 18, 23, 25, 29, 30, 35, 43, 45, 47, 48, 49, 50, 51, 60, 62, 64, 67, 68, 70, 71, 73, 75, 81]; %39 excluded
subjects = [stim, sham];


nruns = 2;
nscans = 250;
sk = 's6';

%% load masks
hdr = spm_vol(fullfile(maskpath,'r3mm_SeedRegion_lOFC.nii'));
nhdr = hdr;
lOFC = round(spm_read_vols(hdr));
hdr = spm_vol(fullfile(maskpath,'r3mm_SeedRegion_rOFC.nii'));
rOFC = round(spm_read_vols(hdr));
OFC = lOFC+rOFC;

hdr = spm_vol(fullfile(maskpath,'r3mm_TargetRegion_lLPFC.nii'));
lLPFC = round(spm_read_vols(hdr));
hdr = spm_vol(fullfile(maskpath,'r3mm_TargetRegion_rLPFC.nii'));
rLPFC = round(spm_read_vols(hdr));
LPFC = lLPFC+rLPFC;

hdr =spm_vol(fullfile(maskpath,'r3mm_rTPM_point1_NoCerebellum.nii'));
yy = round(spm_read_vols(hdr));

hdr =spm_vol(fullfile(maskpath,'r3mm_aal.nii'));
y = round(spm_read_vols(hdr));

hdr = spm_vol(fullfile(maskpath,'r3mm_wm_0.9.nii'));
wm = spm_read_vols(hdr);

hdr = spm_vol(fullfile(maskpath,'r3mm_csf_0.9.nii'));
csf = spm_read_vols(hdr);


%% find voxel index
gmindx{1} = find(lOFC>0);
gmindx{2} = find(rOFC>0);
gmindx{3} = find(OFC>0);
gmindx{4} = find(lLPFC>0);
gmindx{5} = find(rLPFC>0);
gmindx{6} = find(LPFC>0);
nclust = length(gmindx);

%% find wb ROI index
nrois = length(unique(y))-1;
for j=1:nrois
    wbindx{j} = find(y==j);
end

wmindx = find(wm>0);
csfindx = find(csf>0);
glindx = find(yy>0);

for s = subjects
    sn =  sprintf('POETMS_%02d',s);
    datapath = fullfile(root, sn);
    
    
    for r = 1:nruns
        fprintf('Working on %s DAY %01d\n', sn, r);
        
        path = fullfile (datapath, 'nii', sprintf('d%01d', r));
        n = dir(fullfile(path, sprintf('%sw3f*.nii', sk)));
        
        
        % load data;
        try
            load(fullfile(root, sn, sprintf('tc_%s_SEED_STIM_AAL_filtered_3mm_d%01d.mat', sk, r)));

        catch 
        x = zeros(nscans,nclust);
        y = zeros(nscans,nrois);
        global_mean = zeros(nscans,1);
        wm_mean = zeros(nscans,1);
        csf_mean = zeros(nscans,1);

        for i=1:nscans
            tmp = spm_read_vols(spm_vol(fullfile(path, n(i).name)));
            global_mean(i) = mean(tmp(glindx)); % global mean
            wm_mean(i) = mean(tmp(wmindx)); % white matter
            csf_mean(i) = mean(tmp(csfindx)); % CSF

            for j=1:nclust
                x(i,j) = mean(tmp(gmindx{j})); % mean seed/stim cluster activity
            end

            for j=1:nrois
                y(i,j) = mean(tmp(wbindx{j})); % mean wb roi activity
            end

        end

            % filter data
            mreg = load(fullfile(path, sprintf('nusiance_regresssors_%s_run_%01d.txt', sn,  r)));
            
            % make time lagged global mean regressors
            startlag = -5; % startlag in TR
            stoplag = 5;
            global_mean = zscore(global_mean);
            global_mean_lag = [];
            for lag=startlag:1:stoplag
                if sign(lag)==-1
                    global_mean_lag = [global_mean_lag, [global_mean(abs(lag)+1:end); zeros(abs(lag),1)] ];
                elseif lag == 0
                    global_mean_lag = [global_mean_lag, global_mean];
                elseif sign(lag)==+1
                    global_mean_lag = [global_mean_lag, [zeros(abs(lag),1); global_mean(1:end-abs(lag))] ];
                end
            end

            mreg = [zscore([mreg, global_mean_lag,  wm_mean, csf_mean,  [1:nscans]', ]), ones(nscans,1)]; 

            for j=1:nclust
                b = inv(mreg'*mreg)*(mreg'*x(:,j));
                x(:,j) = x(:,j) - mreg*b; 
            end

           for j=1:nrois
                b = inv(mreg'*mreg)*(mreg'*y(:,j));
                y(:,j) = y(:,j) - mreg*b; 
            end
            
            % save data
           save(fullfile(root, sn, sprintf('tc_%s_SEED_STIM_AAL_filtered_3mm_d%01d.mat', sk, r)), 'x', 'y', 'wbindx', 'gmindx');
        end

        cm{r}(:,:,s) = atanh(corrcoef([x, y]));
    end
    TMS_diff(:,:,s) = cm{2}(:,:,s) - cm{1}(:,:,s); 
end

ii = [ mean(TMS_diff(:,:,sham),3), mean(TMS_diff(:,:,stim),3), [mean(TMS_diff(:,:,sham),3)-mean(TMS_diff(:,:,stim),3)] ];
ii = ii(:);

% figure, 
% subplot(1,3,1)
% imagesc(mean(TMS_diff(:,:,sham),3), [min(ii), max(ii)]),title('sham (D2-D1)')
% subplot(1,3,2)
% imagesc(mean(TMS_diff(:,:,stim),3), [min(ii), max(ii)]), title('stim (D2-D1)')
% subplot(1,3,3)
% imagesc(mean(TMS_diff(:,:,sham),3)-mean(TMS_diff(:,:,stim),3), [min(ii), max(ii)]), title('sham(D2-D1)-stim(D2-D1)')
% colorbar


figure, 
subplot(1,3,1)
imagesc(mean(TMS_diff(1:6,1:end,sham),3)', [min(ii), max(ii)]),title('sham (D2-D1)')
subplot(1,3,2)
imagesc(mean(TMS_diff(1:6,1:end,stim),3)', [min(ii), max(ii)]), title('stim (D2-D1)')
subplot(1,3,3)
imagesc(mean(TMS_diff(1:6,1:end,sham),3)'-mean(TMS_diff(1:6,1:end,stim),3)', [min(ii), max(ii)]), title('sham(D2-D1)-stim(D2-D1)')
colorbar



load(fullfile(maskpath, 'aal.mat'))
ROIs = {'left OFC', 'right OFC', 'OFC', 'left LPFC', 'right LPFC', 'LPFC'};
for i=1:6
    for j=1:nrois+6
        [h,p, ci, stats] = ttest2(TMS_diff(i,j,sham), TMS_diff(i,j,stim));
        if h==1 && p<(0.05/116)
            if j<7
                disp(sprintf('p=%0.5f, t=%1.3f, %s, %s', p, stats.tstat, ROIs{i}, ROIs{j}))
            else
                disp(sprintf('p=%0.5f, t=%1.3f, %s, %s, %1.0f', p, stats.tstat,  ROIs{i}, char(name{j-6}), j-6))
            end
        end
    end
end

