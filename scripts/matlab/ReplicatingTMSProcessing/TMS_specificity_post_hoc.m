
% This analysis compares the functional connectivity of the seed region
% with the aOFC target and the pOFC target. Theoretically we expect that
% aOFC seed region is more connected to the aOFC target region, while pOFC
% seed region is more connected to the pOFC target region. But there is no
% guarantee it is the truth. This analysis leverages the GLM used to get
% the TMS coordinates and then extract t-values in both target regions.
% This is only using the Day0 resting-state data.

clear
close all

studydir = '/Volumes/X9Pro/NODEAP/MRI';
SubIDlist = dir(fullfile(studydir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

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

% plot distribution for each subject?
plot_dist = 0;

%%
for subj = 1:nSubIDlist

    sn = SubIDlist(subj).name;
    subdir = fullfile(studydir,sn);
    unnormdir = fullfile(subdir, 'UnnormedMasks');  %path to native space masks for each subejct
    TMSpostdir = fullfile(studydir,'TMS_post_hoc');
    if ~exist(TMSpostdir,'dir')
        mkdir(TMSpostdir)
    end

for i=1:length(regions)
    modeldir = fullfile(subdir,'nii_for_TMS',sprintf('GLM_%s',regions(i).seedName));
    % read connectivitiy t values in each of the target regions
    for t = 1:2
        [targetvol, targetXYZ] = spm_read_vols(spm_vol_nifti(fullfile(unnormdir,sprintf('rw%s.nii',regions(t).targetName))));
        targetidx = find(targetvol==1);
        convol = spm_read_vols(spm_vol_nifti(fullfile(modeldir,'spmT_0001.nii')));
        vals = convol(targetidx);
        convals{i,t} = vals;
        mean_convals(i,t,subj) = mean(convol(targetidx));
        
        % coord with max convals
        loc_max = targetidx(vals==max(vals));
        StimCoord = round(targetXYZ(:,loc_max)');
        
    end
end

%% plot the distribution for each subject
if plot_dist
figure('Position', [100, 100, 700, 600]);
for i = 1:size(convals,1)
    subplot(size(convals,1),1,i);
    hold on;
    for t = 1:size(convals,2)
        histogram(convals{i,t}, 'DisplayName', regions(t).targetName);
    end
    hold off;
    title(regions(i).seedName,'Interpreter', 'none');
    xlabel('Connectivity values');
    l=legend('show');
    set(l, 'Interpreter', 'none')
end
filename = fullfile(TMSpostdir,'Comp_convals_subs',sprintf('Comp_convals_%s.png',sn));
saveas(gcf, filename);
close all
end

        
end


%%
% Extract the four entries from the 3D array
entry1 = reshape(mean_convals(1,1,:),[],1); % aOFC seed - aOFC stim site
entry2 = reshape(mean_convals(1,2,:),[],1); % aOFC seed - pOFC stim site
entry3 = reshape(mean_convals(2,1,:),[],1); % pOFC seed - aOFC stim site
entry4 = reshape(mean_convals(2,2,:),[],1); % pOFC seed - pOFC stim site

% make this to be a matrix
vec_sub = repmat({SubIDlist.name}',4,1);
vec_conn = [entry1;entry2;entry3;entry4];
vec_seed = [repmat({'aOFC seed'},nSubIDlist*2,1);repmat({'pOFC seed'},nSubIDlist*2,1)];
vec_stim = [repmat({'aOFC stim'},nSubIDlist,1);repmat({'pOFC stim'},nSubIDlist,1)];
vec_stim = repmat(vec_stim,2,1);
conn_mat_sub = table(vec_sub,vec_seed,vec_stim,vec_conn,...
    'VariableNames', {'SubID', 'seed', 'stim', 'conn'});
writetable(conn_mat_sub, '../ProcessedData/SeedStimConn_subs.xlsx')

