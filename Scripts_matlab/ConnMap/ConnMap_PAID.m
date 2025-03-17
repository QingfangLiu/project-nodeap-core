
% compute conn maps within the OFC and LPFC clusters, using each voxel as a
% seed
% use 3mm MNI space
% use filtered dat in global conn analysis

clear; clc; close all
studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

nscans = 310; % 310 volumes for ME sequnce
nechos = 3;
TR = 1.5;
TEs = [14.60, 39.04, 63.48];
rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

have_smooth = false; % to use smoothed data or non-smoothed data

%%
userHome = getenv('HOME');
maskpath = fullfile(userHome,'Documents/Projects/NODEAP_data_analysis/Scripts_matlab/GlobalConn/mask');
masks = {'aOFC_seed_right.nii', ...
    'aOFC_target_right.nii', ...
    'pOFC_seed_right.nii', ...
    'pOFC_target_right.nii'
    };  % masks in 3mm space
n_masks = length(masks);

%%
gm_nii = fullfile(maskpath,'gm_0.1_3mm.nii'); % gray matter mask
gm_dat = spm_read_vols(spm_vol(gm_nii));
gm_idx = find(gm_dat > 0);

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
    subdir = fullfile(MRIdir,sn);
    niidir = fullfile(subdir, 'nifti');
    
for r = 1:n_rest_names % to loop over rest_names
        
        curr_rest = rest_names{r};
        fprintf('%s of %s\n',curr_rest,sn)
        
        func_conn_dir = fullfile(studydir, 'ConnMap_PAID', sn); % path to save the ConnMap from PAID result
        if ~exist(func_conn_dir,'dir')
            mkdir(func_conn_dir)
        end
        
        if count_table(subj,r)<1
            disp('Not enough data, skipping this iteration ...')
            continue;
        end
        
        if have_smooth
            load(fullfile(studydir, 'GlobalConn_PAID', sn, curr_rest, 'tc_filtered_3mm.mat')); 
        else
            load(fullfile(studydir, 'GlobalConn_PAID', sn, curr_rest, 'tc_filtered_3mm_no_sm.mat')); 
        end
        
        % reduce the 'dat' to voxels we care about
        dat = dat(:,all_idx(:,2));
        dat = zscore(dat); % zcore for convenient Pearson correlation calculation
                
%%
    corr_map = nan(nvox,nvox);
    for i = 1:nvox
        corr_map(i,:) = (dat(:,i)' * dat)./ nscans; % Pearson correlation on z-scored values (with each voxel)
        if ~mod(i,50)
            fprintf('\t\t%d of %d\n', i, nvox);
        end
    end
    
    fig = figure('Visible', 'off'); % plot without showing
    imagesc(corr_map);
    colorbar; % Add a colorbar to show the scale
    caxis([-1, 1]); % Adjust these limits based on your data
    title(sprintf('Voxel-voxel connectivity map of %s',curr_rest));
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
    
    if have_smooth
        note = '';
    else
        note = '_no_sm';
    end
    
    img_name = fullfile(func_conn_dir,sprintf('conn_map_%s_%s.bmp',curr_rest,note));
    saveas(gcf, img_name);
    
    % save this corr mat 
    mat_name = fullfile(func_conn_dir,sprintf('conn_map_%s_%s.mat',curr_rest,note));
    save(mat_name,'corr_map')
    
    close all    
    
end
end







