
%% ========================================================================
%  Script: make_nuisance_regressors_motion_slice.m
%
%  Purpose
%    1) For each subject/session, load SPM realignment params (from e1) and
%       build motion-derived nuisance regressors:
%         - 6 motion params (X,Y,Z, pitch, roll, yaw)
%         - temporal diffs (Δ), squared terms (²), squared diffs (Δ²)
%    2) Compute slice-based artifacts from echo-1 (preprocessing space):
%         - odd–even slice mean difference (per volume)
%         - across-slice variance (per volume)
%       Then z-score across ALL sessions (subject-level) and flag outliers.
%    3) Save per-session regressors (z-scored) and diagnostic plots.
%
%  Inputs (per subject)
%    <dat_folder>/<SubID>/nifti/<Sess>_me/  (contains rp_*.txt, <Sess>*_e1.nii(.gz))
%    Availability matrix: /Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx
%      - rows: SubID order from directory listing
%      - cols: {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'}
%      - values: 1 (310 vols), 0.5 (205 vols), 0 (missing)
%
%  Outputs
%    NRegressor/<SubID>/nuisance_regressors_<Sess>.txt   % z-scored regressors
%    NRegressor/<SubID>/SliceDiffVar.bmp                 % diagnostics figure
%    NRegressor/<SubID>/slice_diff_var.mat               % raw slice metrics
%    NRegressor/<SubID>/BadVolumes.mat                   % indices and counts
%
%  Notes
%    - Motion parameters are taken from e1 (as realignment was run on e1).
%    - Slice metrics adapt to the actual #slices (e.g., 60 vs 78).
%    - Thresholds for outlier detection are fixed at 4 SD (modifiable).
% ========================================================================

%% 0) Setup ----------------------------------------------------------------
clear; clc;

dat_folder = '/Volumes/X9Pro/NODEAP/MRI';
NR_dir = '/Volumes/X9Pro/NODEAP/NRegressor';
if ~exist(NR_dir, 'dir'); mkdir(NR_dir); end

% Subject list
SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);
nSub      = numel(SubIDlist);

% Availability table (rows=subjects matching SubIDlist order; cols=rest_names)
count_xlsx  = '/Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx';
count_table = readmatrix(count_xlsx);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

if_plot_motion = true;

% Outlier thresholds (in SD after subject-level z-scoring)
crit = [4,4];        % criterion for slicediff and slicevar, eyeball it looking at the peaks

%% 1) Loop: subjects -------------------------------------------------------
for subj = 1:nSub

SubID = SubIDlist(subj).name;
SubDir = fullfile(dat_folder,SubID);
niidir = fullfile(SubDir, 'nifti');

% Per-subject containers
slicediff = cell(1,n_rest_names);
slicevar = cell(1,n_rest_names);
adj_slicediff = cell(1,n_rest_names);
adj_slicevar = cell(1,n_rest_names);
NR = cell(1,n_rest_names);
badvol = cell(1,n_rest_names);      % bad volumes per session
ctr_badvol = zeros(n_rest_names,1); % count # of bad volumes
ctr_NR = zeros(n_rest_names,1);     % count # of nuisance regressors

% Output folder for this subject
NR_sub_dir = fullfile(NR_dir, SubID);
if ~exist(NR_sub_dir, 'dir'); mkdir(NR_sub_dir); end

fprintf('\n================ %s ================\n', SubID);
disp('Working on motion parameters and slice metrics ...');

%% 2) Loop: sessions ----------------------------------------------------
for r = 1:n_rest_names 
    curr_rest = rest_names{r};
    
    % Skip if marked missing in count_table
    val = count_table(subj, r+1);
    if val == 0
        continue;
    elseif val == 0.5
        nscans = 205;
    else
        nscans = 310;
    end
        
    fprintf('%s — %s\n', SubID, curr_rest);
    
    %% 2a) Motion regressors (from rp_*.txt)
    funcRunDir = fullfile(niidir,sprintf('%s_me',curr_rest));
    mpfilename = dir(fullfile(funcRunDir, 'rp_*.txt'));    %rp file contains 6 columns represnting positions  
    mp = load(fullfile(funcRunDir, mpfilename(1).name));          % mp motion parameters
        
    % Build regressors: 6, Δ6, 6², (Δ6)²  => 24 columns
    NR{r} = [mp, [zeros(1,6); diff(mp)], mp.^2, [zeros(1,6); diff(mp).^2]];

    % Optional motion plots
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
        movefile(tmpfigname,NR_sub_dir)
    end

    %% 2b) Slice metrics from echo-1 (preprocessing space)
    % Find e1 NIfTI (decompress if needed)
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

%% 3) Aggregate slice metrics across sessions (subject-level z-scoring)
save(fullfile(NR_sub_dir,'slice_diff_var.mat'), 'slicediff', 'slicevar');   % save the two variables for later access

all_slicediff = [slicediff{:}]; % aggregate all slice differences across sessions and runs
all_slicevar = [slicevar{:}]; % aggregate all slice variances across sessions and runs
aztmp = [max(abs(zscore(all_slicediff))),max(abs(zscore(all_slicevar)))];


%% 4) Per-session outliers, add regressors, save -----------------------

    close all
    figure
    
    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
        
        % Skip if session never processed
        if isempty(NR{r}) || isempty(slicediff{r}) || isempty(slicevar{r})
            continue;
        end
    
        adj_slicediff{r} = abs((slicediff{r} - mean(all_slicediff))./std(all_slicediff));
        adj_slicevar{r} = abs((slicevar{r} - mean(all_slicevar))./std(all_slicevar));

        % Panels
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
            
        % Append slice metrics and their diffs (and squared diffs) + bad vols
        NR{r} = [NR{r},adj_slicediff{r}', adj_slicevar{r}',...
            [0; diff(adj_slicediff{r})'], [0; diff(adj_slicevar{r})'],...
            [0; diff(adj_slicediff{r})'].^2, [0; diff(adj_slicevar{r})'].^2,...
            bv];   
        ctr_NR(r) = size(NR{r},2);
        
        % save z-scored nuisance regressors as txt
        NRfilename = fullfile(NR_sub_dir,sprintf('nuisance_regressors_%s.txt', curr_rest));
        dlmwrite(NRfilename, zscore(NR{r}), 'delimiter', ' ', 'precision', '%.24f');   %z scoring
    
end % end of rest session

        % Save Slice diff and var figures 
        set(gcf, 'PaperPosition', [2 1 12 16]);
        saveas(gcf,fullfile(NR_sub_dir,'SliceDiffVar.bmp'))
        
        % save variables for current subject: badvol, ctr_badvol, ctr_NR
        save(fullfile(NR_sub_dir,'BadVolumes.mat'), 'badvol', 'ctr_badvol', 'ctr_NR');  
        
        close all

end % subject

