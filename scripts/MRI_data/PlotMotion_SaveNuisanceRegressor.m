
% This code
% 1. plot realignment motion parameters across volumes, for each session
% and each subject
% 2. do within-volumn calculation 
% for multi-echo data, motion estimated from e1
% 3. saves everything as nuisance regressor

%% setup
studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

MRIcount = readtable('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx',"ReadRowNames",true);
if_plot_motion = true;

%%
for subj = 1:nSubIDlist

SubID = SubIDlist(subj).name;
SubDir = fullfile(MRIdir,SubID);
niidir = fullfile(SubDir, 'nifti');

% initiate some variables for curr subject
slicediff = cell(1,n_rest_names);
slicevar = cell(1,n_rest_names);
adj_slicediff = cell(1,n_rest_names);
adj_slicevar = cell(1,n_rest_names);
NR = cell(1,n_rest_names);
badvol = cell(1,n_rest_names);      % bad volumes per session
ctr_badvol = zeros(n_rest_names,1); % count # of bad volumes
ctr_NR = zeros(n_rest_names,1);     % count # of nuisance regressors

% where to save nuisance regressors
NROutdir = fullfile(studydir,'NRegressor',SubID);
if ~exist(NROutdir,'dir')
    mkdir(NROutdir)
end

disp('Working on motion parameters ...')

for r = 1:n_rest_names % to loop over rest_names
    curr_rest = rest_names{r};
    if MRIcount{SubID,curr_rest}==0
        continue;
    end
    fprintf('%s of %s\n',curr_rest,SubID)
    
    if MRIcount{SubID,curr_rest}==0.5
        nscans = 205;
    else
        nscans = 310;
    end
    
    % find the motion parameters
    funcRunDir = fullfile(niidir,sprintf('%s_me',curr_rest));
    mpfilename = dir(fullfile(funcRunDir, 'rp_*.txt'));    %rp file contains 6 columns represnting positions  
    mp = load(fullfile(funcRunDir, mpfilename(1).name));          % mp motion parameters
        
    % add diff, squared mp, and squared diff
    NR{r} = [mp, [zeros(1,6); diff(mp)], mp.^2, [zeros(1,6); diff(mp).^2]];

    if if_plot_motion
        figure;
        subplot(2,4,1); plot(NR{r}(:,1:3)); xlabel('scans'); ylabel('mm'); 
        subplot(2,4,2); plot(NR{r}(:,4:6)); xlabel('scans'); ylabel('radian');
        subplot(2,4,3); plot(NR{r}(:,7:9)); xlabel('scans'); ylabel('diff(mm)'); 
        subplot(2,4,4); plot(NR{r}(:,10:12)); xlabel('scans'); ylabel('diff(radian)');
        subplot(2,4,5); plot(NR{r}(:,13:15)); xlabel('scans'); ylabel('mm^2');
        subplot(2,4,6); plot(NR{r}(:,16:18)); xlabel('scans'); ylabel('radian^2');
        subplot(2,4,7); plot(NR{r}(:,19:21)); xlabel('scans'); ylabel('diff(mm)^2'); 
        subplot(2,4,8); plot(NR{r}(:,22:24)); xlabel('scans'); ylabel('diff(radian)^2');
        set(gcf, 'PaperPosition', [2 1 18 8]);
    
        % save the plot
        tmpfigname = fullfile(funcRunDir,sprintf('MotionParameters_%s_%s.bmp',SubID,curr_rest));
        saveas(gcf,tmpfigname)
        movefile(tmpfigname,NROutdir)
    end

%% analyze slice differences and variances
% this is done using e1 before any processing (as we did realignment using e1)
    idx1 = 1:2:60;   % odd slices
    idx2 = 2:2:60;   % even slices

    funcfilenames = dir(fullfile(funcRunDir,sprintf('%s*_e1.nii',curr_rest))); 
    if isempty(funcfilenames)
        gzfilenames = dir(fullfile(funcRunDir,sprintf('%s*_e1.nii.gz',curr_rest))); 
        gunzip(fullfile(funcRunDir,gzfilenames(1).name))
        funcfilenames = dir(fullfile(funcRunDir,sprintf('%s*_e1.nii',curr_rest))); 
    end
    
    y = spm_vol(fullfile(funcRunDir,funcfilenames(1).name)); % 310*1 struct
    
    for vol = 1:nscans
        v = spm_read_vols(y(vol));   %dimensions of v show that there are 60 slices
        odd_slice = nanmean(nanmean(nanmean(v(:,:,idx1),3)));
        even_slice = nanmean(nanmean(nanmean(v(:,:,idx2),3)));
        slicediff{r}(vol) = odd_slice - even_slice; % diff b/t odd and even slices
        slicevar{r}(vol) = var(nanmean(squeeze(nanmean(v)))); % var across all 60 slices
    end
    
end % 

save(fullfile(NROutdir,'slice_diff_var.mat'), 'slicediff', 'slicevar');   % save the two variables for later access

all_slicediff = [slicediff{:}]; % aggregate all slice differences across sessions and runs
all_slicevar = [slicevar{:}]; % aggregate all slice variances across sessions and runs
aztmp = [max(abs(zscore(all_slicediff))),max(abs(zscore(all_slicevar)))];
crit = [4,4];        % criterion for slicediff and slicevar, eyeball it looking at the peaks

disp('More plotting ...')

close all
figure

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
    
    fprintf('%s of %s\n',curr_rest,SubID)
    
        adj_slicediff{r} = abs((slicediff{r} - mean(all_slicediff))./std(all_slicediff));
        adj_slicevar{r} = abs((slicevar{r} - mean(all_slicevar))./std(all_slicevar));

        sh = subplot(n_rest_names, 2, r*2-1);    
        plot(adj_slicediff{r})
        set(sh, 'YLim', [0,max(aztmp)+1]);
        title(sprintf('slice diff, %s',curr_rest));
        hold on
        plot(find(adj_slicediff{r}> crit(1)), adj_slicediff{r}(adj_slicediff{r}> crit(1)), '*r');
        
        sh = subplot(n_rest_names, 2, r*2);
        plot(adj_slicevar{r})
        set(sh, 'YLim', [0,max(aztmp)+1]);
        title(sprintf('slice var, %s',curr_rest));
        hold on
        plot(find(adj_slicevar{r}> crit(2)), adj_slicevar{r}(adj_slicevar{r}> crit(2)), '*r');
        
        % include regressors for bad volumes
        tmpbadvol = unique([find(adj_slicediff{r}> crit(1)), find(adj_slicevar{r}> crit(2))]);
        badvol{r} = tmpbadvol;
        ctr_badvol(r) = length(tmpbadvol);         % count bad volumes
        bv = zeros(nscans,ctr_badvol(r));
        for i = 1:ctr_badvol(r)
            bv(tmpbadvol(i),i) = 1;
        end
            
        % add slicediff and diff of slicediff and bad volume regressors
        NR{r} = [NR{r},adj_slicediff{r}', adj_slicevar{r}',...
            [0; diff(adj_slicediff{r})'], [0; diff(adj_slicevar{r})'],...
            [0; diff(adj_slicediff{r})'].^2, [0; diff(adj_slicevar{r})'].^2,...
            bv];   
        ctr_NR(r) = size(NR{r},2);
        
        % save z-scored nuisance regressors as txt
        NRfilename = fullfile(NROutdir,sprintf('nuisance_regressors_%s.txt', curr_rest));
        dlmwrite(NRfilename, zscore(NR{r}), 'delimiter', ' ', 'precision', '%.24f');   %z scoring
    
end % end of rest session

        % Save Slice diff and var figures 
        set(gcf, 'PaperPosition', [2 1 12 16]);
        saveas(gcf,fullfile(NROutdir,'SliceDiffVar.bmp'))
        
        % save variables for current subject: badvol, ctr_badvol, ctr_NR
        save(fullfile(NROutdir,'BadVolumes.mat'), 'badvol', 'ctr_badvol', 'ctr_NR');  
        
        close all

end % subject

