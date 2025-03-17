%% Main lag script used to compute time delay matrix and lag projection
%
%% General information for a few variables and their expected file types
%
%  subjects:  various compatible file formats; read as string array
%  BOLD:    .mat file; Time x Space (where the spatial dimension corresponds to a vector of spatial units...in other words, BOLD is a 2D matrix, where spatial units are vectorized)
%  good:    .mat file; Vector with length = number of spatial units; used to mask out bad spatial units
%  format:  various compatible file formats (e.g., text file); vector with length = number of time points; used to mask out bad time points
%%%%%%%%
%
clc; clear; close all
studydir = '/Volumes/X9Pro/NODEAP';
userHome = getenv('HOME');
maskpath = fullfile(userHome,'Library/CloudStorage/Box-Box/NODEAP_data_analysis/Scripts_matlab/GlobalConn/mask');

MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);
masks = {'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'
    };  % masks in 3mm space
n_masks = length(masks);
nscans = 310;

gm_nii = fullfile(maskpath,'gm_0.1_3mm.nii'); % gray matter mask
gm_dat = spm_read_vols(spm_vol(gm_nii));
gm_idx = find(gm_dat > 0);
n_gm_vox = length(gm_idx);

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

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

num_nodes = nvox;    % number of time series (i.e., spatial dimension)

% this analysis is done in the multiple of TRs
lag_lim = 4;    % lag limit (in seconds)
lags = -4:4;    % range of TR shifts; max(lags) = round(lag_lim/tr + 1)
tr = 1.5;       % sampling interval in seconds (changed to 1.5s)

min_block_durn = (max(lags)+1)*tr;   % min. block duration (in seconds)

%% Loop over subjects

for subj = 1:nSubIDlist
    
    tic
    SubID = SubIDlist(subj).name;
    disp(['Processing ' SubID]);

    SubDir = fullfile(MRIdir,SubID);
    
    % load info of bad volumes
    NROutdir = fullfile(studydir,'NRegressor',SubID);
    load(fullfile(NROutdir,'BadVolumes.mat'))
    
    lag_out_dir = fullfile(studydir, 'lag_out_PAID', SubID); 
    if ~exist(lag_out_dir,'dir')
        mkdir(lag_out_dir)
    end

    for r = 1:n_rest_names
    
    if count_table(subj,r)<1
        disp('Not enough data, skipping this iteration ...')
        continue;
    end
    
    curr_rest = rest_names{r};
    path = fullfile (SubDir, 'nifti', 'functional', curr_rest);
    n = dir(fullfile(path, 's6w3*.nii')); % data with 3mm voxel size, smoothed
    dat = zeros(nscans,n_gm_vox); % initiliaze to get data
    data_4D = spm_read_vols(spm_vol(fullfile(path, n.name)));              
    for i = 1:nscans
        tmp = data_4D(:,:,:,i);
        dat(i,:) = tmp(gm_idx); % gray matter
    end
        
    % initialize subject matrices (num_nodes * num_nodes)
    peak_lags = single(nan(num_nodes)); % peak lags
    r_ZL = peak_lags;   % zero-lag correlation
    r_peak = peak_lags; % peak correlation (correlation at optimal lag)

    BOLD = dat(:,all_idx(:,2));   % only use those voxels from the ROIs
    
    good = true(1, num_nodes);    % a logical vector of spatial dimension (no use here)
    
    format = true(1, nscans);     % a logical vector of temporal dim
    format(badvol{1,r}) = false;  % mark bad volumes as false

    % separating the whole run into blocks using location of bad volumes
    FORMAT = create_blocks(format,min_block_durn,tr);
    nblocks = numel(FORMAT);

    %% Do the lagged correlation/covariance computation of TD matrices
    Cov = single(zeros([sum(good) sum(good) numel(lags)]));
    nframes = 0;  % count number of frames
    
    % De-mean time series
    run_mean = nanmean(BOLD(format,:),1);
    BOLD = BOLD - run_mean;
        
    % Loop over blocks of contiguous frames
    % adds up covariance across frames
    for j = 1:numel(FORMAT)
        nframes = nframes + numel(FORMAT{j});
        FHCR = false(1,nscans);
        FHCR(FORMAT{j}) = true; % mark which volumes to include
        Cov = Cov + lagged_cov(BOLD(FHCR,good),BOLD(FHCR,good),max(lags));
    end
    
    % Normalize pairwise cross-covariance functions based on entire run
    for k = 1:numel(lags)
        Cov(:,:,k) = Cov(:,:,k)/(nframes - abs(lags(k))*nblocks);
    end
    
    % Parabolic interpolation to get peak lag/correlation
    [pl,pc] = parabolic_interp(Cov,tr);  % pl on second unit
    pl(abs(pl) > lag_lim) = nan; % Exclude long lags (generally occur when CCF is flat)
    
    % Get zero-lag correlation
    temp = Cov(:,:,lags==0);  % zero-lag correlation
    d = zeros(size(temp));
    d(logical(eye(length(temp)))) = sqrt(diag(temp));
    temp = d^(-1)*temp/d;
    temp = atanh(temp); % Fisher z transform
    temp(isnan(pl)) = nan;
    
    % now we have the lags, ZL corr, and peak corr
    % save for the current session of rs-fMRI data
    peak_lags(good,good) = pl;
    % peak_lags: diagonal=0, anti-symmetrical
    r_ZL(good,good) = temp;
    r_peak(good,good) = pc;
    
    % plot r_peak and r_ZL against each other
    figure
    subplot(1,2,1)
    imagesc(r_peak);
    colorbar; % Add a colorbar to show the scale
    subplot(1,2,2)
    imagesc(r_ZL);
    colorbar; % Add a colorbar to show the scale
    
    
    fig = figure('Visible', 'off'); % plot without showing
    imagesc(peak_lags);
    colorbar; % Add a colorbar to show the scale
    caxis([-4, 4]); % Adjust these limits based on your data
    title(sprintf('Voxel-voxel connectivity peak lags of %s',curr_rest));
    xlabel('');
    ylabel('');
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
    colormap(jet); % 
    
    % Determine the boundaries for each group
    group_boundaries = find(diff(all_idx(:,3))) + 1;

    % Add lines to separate the groups
    hold on;
    for i = 1:length(group_boundaries)
        x = group_boundaries(i) - 0.5;
        y = group_boundaries(i) - 0.5;
        plot([x, x], [0.5, 349.5], 'k', 'LineWidth', 2);
        plot([0.5, 349.5], [y, y], 'k', 'LineWidth', 2);
    end
    
    tick_positions = [mean([1, group_boundaries(1)]), ...
                  mean([group_boundaries(1), group_boundaries(2)]), ...
                  mean([group_boundaries(2), group_boundaries(3)]), ...
                  mean([group_boundaries(3), 349])];

    labels = {'aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'};
    set(gca, 'XTick', tick_positions, 'XTickLabel', labels, ...
         'YTick', tick_positions, 'YTickLabel', labels, ...
         'TickLength', [0 0], 'FontSize', 12, 'FontWeight', 'bold', 'XColor', 'k', 'YColor', 'k');

    % Adjust tick label positions for better visualization
    ax = gca;
    ax.XAxis.TickLabelInterpreter = 'none';
    ax.YAxis.TickLabelInterpreter = 'none';
    ax.XAxisLocation = 'top'; % Move X axis labels to the top

    hold off;
    
    img_name = fullfile(lag_out_dir,sprintf('peak_lags_%s.bmp',curr_rest));
    saveas(gcf, img_name);
    
    % save variables 
    mat_name = fullfile(lag_out_dir,sprintf('lag_corr_%s.mat',curr_rest));
    save(mat_name,'peak_lags', 'r_ZL', 'r_peak')
    
    end
    
    toc

end

%%
% average lag matrix across sessions & across subjects

rctr = 0;
all_all_peak_lags = single(nan([num_nodes num_nodes n_rest_names*nSubIDlist]));

for subj = 1:nSubIDlist
    
    SubID = SubIDlist(subj).name;
    disp(['Processing ' SubID]);
    
    lag_out_dir = fullfile(studydir, 'lag_out_PAID', SubID); 
    all_peak_lags = single(nan([num_nodes num_nodes n_rest_names]));
    
    for r = 1:n_rest_names
        rctr = rctr + 1;
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        curr_rest = rest_names{r};
        mat_name = fullfile(lag_out_dir,sprintf('lag_corr_%s.mat',curr_rest));
        load(mat_name)
        all_peak_lags(:,:,r) = peak_lags;
        all_all_peak_lags(:,:,rctr) = peak_lags;
    end
    
    avg_peak_lags = nanmean(all_peak_lags,3);
    
    %% plot
    
    fig = figure('Visible', 'off'); % plot without showing
    imagesc(avg_peak_lags);
    colorbar; % Add a colorbar to show the scale
    %caxis([-4, 4]); % Adjust these limits based on your data
    title('Voxel-voxel connectivity peak lags of all sessions');
    xlabel('');
    ylabel('');
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
    colormap(jet); % 
    
    % Determine the boundaries for each group
    group_boundaries = find(diff(all_idx(:,3))) + 1;

    % Add lines to separate the groups
    hold on;
    for i = 1:length(group_boundaries)
        x = group_boundaries(i) - 0.5;
        y = group_boundaries(i) - 0.5;
        plot([x, x], [0.5, 349.5], 'k', 'LineWidth', 2);
        plot([0.5, 349.5], [y, y], 'k', 'LineWidth', 2);
    end
    
    tick_positions = [mean([1, group_boundaries(1)]), ...
                  mean([group_boundaries(1), group_boundaries(2)]), ...
                  mean([group_boundaries(2), group_boundaries(3)]), ...
                  mean([group_boundaries(3), 349])];

    labels = {'aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'};
    set(gca, 'XTick', tick_positions, 'XTickLabel', labels, ...
         'YTick', tick_positions, 'YTickLabel', labels, ...
         'TickLength', [0 0], 'FontSize', 12, 'FontWeight', 'bold', 'XColor', 'k', 'YColor', 'k');

    % Adjust tick label positions for better visualization
    ax = gca;
    ax.XAxis.TickLabelInterpreter = 'none';
    ax.YAxis.TickLabelInterpreter = 'none';
    ax.XAxisLocation = 'top'; % Move X axis labels to the top

    hold off;
    
    img_name = fullfile(lag_out_dir,'peak_lags_all.bmp');
    saveas(gcf, img_name);
    
end

%% plot across all subs and all sessions
avg_peak_lags = nanmean(all_all_peak_lags,3);
lag_out_dir = fullfile(studydir, 'lag_out_PAID', 'agg'); 
if ~exist(lag_out_dir,'dir')
    mkdir(lag_out_dir)
end
    
fig = figure('Visible', 'off'); % plot without showing
imagesc(avg_peak_lags);
colorbar; % Add a colorbar to show the scale
%caxis([-4, 4]); % Adjust these limits based on your data
title('Voxel-voxel connectivity peak lags of all sessions & all subjects');
xlabel('');
ylabel('');
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
colormap(jet); % 

% Determine the boundaries for each group
group_boundaries = find(diff(all_idx(:,3))) + 1;

% Add lines to separate the groups
hold on;
for i = 1:length(group_boundaries)
    x = group_boundaries(i) - 0.5;
    y = group_boundaries(i) - 0.5;
    plot([x, x], [0.5, 349.5], 'k', 'LineWidth', 2);
    plot([0.5, 349.5], [y, y], 'k', 'LineWidth', 2);
end

tick_positions = [mean([1, group_boundaries(1)]), ...
              mean([group_boundaries(1), group_boundaries(2)]), ...
              mean([group_boundaries(2), group_boundaries(3)]), ...
              mean([group_boundaries(3), 349])];

labels = {'aOFC-seed', 'aOFC-stim', 'pOFC-seed', 'pOFC-stim'};
set(gca, 'XTick', tick_positions, 'XTickLabel', labels, ...
     'YTick', tick_positions, 'YTickLabel', labels, ...
     'TickLength', [0 0], 'FontSize', 12, 'FontWeight', 'bold', 'XColor', 'k', 'YColor', 'k');

% Adjust tick label positions for better visualization
ax = gca;
ax.XAxis.TickLabelInterpreter = 'none';
ax.YAxis.TickLabelInterpreter = 'none';
ax.XAxisLocation = 'top'; % Move X axis labels to the top

hold off;

img_name = fullfile(lag_out_dir,'peak_lags_all.bmp');
saveas(gcf, img_name);
    
    


