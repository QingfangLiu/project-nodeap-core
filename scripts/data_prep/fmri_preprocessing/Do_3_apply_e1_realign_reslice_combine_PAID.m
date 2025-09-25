%% ========================================================================
%  Script: apply_e1_realign_reslice_combine_PAID.m
%
%  Purpose
%    (1) Apply realignment (voxel-to-world) from echo-1 to echoes 2–3
%    (2) Reslice volumes for all echoes (writes 'r*.nii')
%    (3) Combine echoes using PAID weights (w ∝ tSNR * TE) -> fvol_###.nii
%
%  Assumed layout (per subject)
%    <dat_folder>/<SubID>/
%       ├─ nifti/
%       │   ├─ <Sess>_me/                % contains <Sess>_rest*_e[1-3].nii
%       │   ├─ functional/<Sess>/        % combined output fvol_###.nii
%       │   ├─ tSNRmaps/                 % per-echo tSNR maps (written here)
%       │   └─ w_tSNR/                   % w ∝ tSNR*TE maps (written here)
%       └─ MRIjobs/
%
%  Notes
%    - Default nscans = 310; NODEAP_44/S1D1 = 205
%    - Realignment parameters are assumed to be estimated already on e1
%    - This script writes resliced images (prefix 'r') for ALL echoes
%    - PAID combine: weights w_j = (tSNR_j * TE_j) / sum_k (tSNR_k * TE_k)
% ========================================================================

%% 0) Setup -----------------------------------------------------------------
clear; clc;

% Echo metadata
TEs    = [14.60, 39.04, 63.48];
nechos = 3;

% Data root (fallback to prior location if needed)
dat_folder = '/Volumes/X9Pro/NODEAP/MRI';

% Subjects
SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);
nSub      = numel(SubIDlist);

% Sessions (order matters)
rest_names   = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = numel(rest_names);

% SPM defaults (keep your memory settings)
spm('Defaults','fMRI');
global defaults
defaults.stats.maxmem = 2 ^ 34;   % ~16GiB max for stats ops (based on 32GB RAM)
defaults.stats.resmem = true;     % allow temp arrays in memory
spm_jobman('initcfg');

%% 1) Loop over subjects ----------------------------------------------------
for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;
    SubDir = fullfile(dat_folder,SubID);
    niidir = fullfile(SubDir, 'nifti');
    jobdir = fullfile(SubDir, 'MRIjobs');
    
    fprintf('\n================ %s ================\n', SubID);

%% 2) Loop over rest sessions ------------------------------------------
for r = 1:n_rest_names % to loop over rest_names
    
    curr_rest = rest_names{r};
    medir = fullfile(niidir, sprintf('%s_me',curr_rest)); % path to me data
    n_e1 = dir(fullfile(medir, sprintf('%s_rest*_e1.nii',curr_rest)));
    functdir = fullfile(niidir, 'functional', curr_rest); % path to rsfRMI images
    
    fprintf('%s of %s\n',curr_rest,SubID)
    
    % this session has fewer volumes due to scanning cutoff
    if strcmp(SubID,'NODEAP_44') && strcmp(curr_rest,'S1D1')
        nscans = 205;
    else
        nscans = 310;
    end
    
    % skip current session if the functional dir has all combined scans
    if numel(dir(fullfile(functdir, 'fvol_*.nii')))==nscans 
        fprintf('functional images have been converted for %s\n',functdir)
        continue;
    end
    
    if ~exist(medir,'dir')
        continue;
    end
    
    %% (1) Apply e1 realignment transforms to e2/e3 headers -------------
    for i = 1:nscans
        V = spm_get_space(sprintf('%s,%01d', fullfile(medir, n_e1(1).name),i));
        for iecho = 2:nechos
            n = dir(fullfile(medir, sprintf('%s_rest*%d.nii',curr_rest,iecho))); 
            spm_get_space(sprintf('%s,%01d', fullfile(medir, n(1).name),i),V);
        end
    end

%% (2) Reslice all echoes (writes 'r*.nii') -------------------------
fimages = [];
func_ctr = 0;
rfilename = [];

for iecho = 1:nechos
    n = dir(fullfile(medir, sprintf('%s_rest*_e%d.nii',curr_rest,iecho)));
    for i = 1:nscans
        func_ctr = func_ctr + 1;
        fimages{func_ctr,1} = sprintf('%s,%01d', fullfile(medir, n(1).name),i);
        rfilename{iecho}{i,1} = sprintf('%s,%01d', fullfile(medir, sprintf('r%s', n(1).name)),i);
    end
end

matlabbatch = [];
matlabbatch{1}.spm.spatial.realign.write.data = fimages;
matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 0]; % reslice all images, but no need to create the mean
matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';

% save and run job
fname = fullfile(jobdir, sprintf('me_reslice_all_echoes_%s.mat',curr_rest));
save(fname, 'matlabbatch');
spm_jobman('run', matlabbatch);

    
%% (3) PAID combine: tSNR, weights, fused fvol_###.nii --------------
fprintf('  ↳ Computing tSNR & PAID weights; combining echoes ...\n');

tSNPmapsdir = fullfile(niidir,'tSNPmaps');
if ~exist(tSNPmapsdir,'dir')
    mkdir(tSNPmapsdir)
end

hdr2 = spm_vol(fimages{1,1}); % get one header
hdr2.dt = [16,0];
tSNR = [];
tSNR_TE = [];
for j = 1:nechos

    tmp = [];
    % load data
    x = [];
    for i = 1:nscans
        hdr = spm_vol(rfilename{j}{i,1});
        x(:,:,:,i) = spm_read_vols(hdr);
    end

    % compute and save tSNR map
    tSNR{j} = mean(x,4)./std(x, [], 4);
    hdr2.fname = fullfile(tSNPmapsdir, sprintf('tSNR_echo_%d_%s.nii', j,curr_rest));
    spm_write_vol(hdr2, tSNR{j});
    
    % multiply tSNR by TE for each echo
    tSNR_TE(:,:,:,j) = tSNR{j}.*TEs(j);
end
clear x

w_tSNRdir = fullfile(niidir,'w_tSNR');
if ~exist(w_tSNRdir,'dir')
    mkdir(w_tSNRdir)
end
    
% compute weight (w)
w = [];
for j = 1:nechos

    % normalize tSNR*TE maps by sum of tSNR*TE's across echos
    w{j} = tSNR_TE(:,:,:,j)./sum(tSNR_TE,4);

    % save w_tSNR map
    hdr2.fname = fullfile(w_tSNRdir, sprintf('w_tSNR_TE_echo_%d_%s.nii', j,curr_rest));
    spm_write_vol(hdr2, w{j});
end

if ~exist(functdir,'dir')
    mkdir(functdir)
end

% combine echos using w
for i = 1:nscans
    y = zeros(size(w{j}));
    for j=1:nechos
        hdr = spm_vol(rfilename{j}{i,1});
        x = spm_read_vols(hdr);
        y = y + (x.*w{j});
    end

    % save combined volume
    hdr.n = [1 1];
    hdr.fname = fullfile(functdir, sprintf('fvol_%03d.nii', i));
    spm_write_vol(hdr, y);
end

end % end of sessions
end % end of subjects



