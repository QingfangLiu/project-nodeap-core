
%% this was done at the early stage of MRI data processing
% to convert dicom data to nifti (3D)
% note that later dicom data was deleted to save space
% also the 3D nifti files were converted to 4D in another analysis code

tic
% this runs on my mac mini
clear; clc;
HomeDir = '/Volumes/X9Pro/NODEAP/MRI';
dcm2niipath = '/Applications/MRIcroGL.app/Contents/Resources/dcm2niix';
nscans = 310; % 310 volumes for ME sequnce
TEs = [14.60, 39.04, 63.48];
nechos = 3;
SubIDlist = dir(fullfile(HomeDir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

OriDir = '/Users/liuq13/NODEAP/MRI';

%%
for subj = 1 %:nSubIDlist

SubID = SubIDlist(subj).name;
SubDir = fullfile(HomeDir,SubID);
niidir = fullfile(SubDir, 'nifti');
jobdir = fullfile(SubDir, 'MRIjobs');

if ~exist(jobdir,'dir') 
   mkdir(jobdir)
end

%% convert dicom to nifti files
% nifti files are saved within a new folder 'nii'
% convert T1 file
anatpath = fullfile(OriDir, SubID, 'dicoms', 'D0_T1');
outdir = fullfile(niidir, 'anat');
if ~exist(outdir,'dir')
    mkdir(outdir);
    system(sprintf('%s -o %s %s', dcm2niipath, outdir, anatpath))
end

%%
% convert resting-state data
rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);
rctr = 0; % count rest sessions
filename = [];

for r = 1:n_rest_names % to loop over rest_names
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubDir)
    
    rctr = rctr + 1; % to account for potential missing sessions
    
    % this session has fewer volumes due to scanning cutoff
    if strcmp(SubID,'NODEAP_44') && strcmp(curr_rest,'S1D1')
        nscans = 205;
    else
        nscans = 310;
    end
    
    func_dicoms_path = fullfile(OriDir, SubID, 'dicoms',sprintf('%s_rest',curr_rest));
    medir = fullfile(niidir, sprintf('%s_me',curr_rest)); % path to me data
    
    % skip current session if no dicom data
    if ~exist(func_dicoms_path,'dir') || numel(dir(func_dicoms_path)) <=2
        fprintf('No such dicom data %s\n',func_dicoms_path)
        continue;
    end
    
    if ~exist(medir,'dir')
        mkdir(medir);
        system(sprintf('%s -o %s %s', dcm2niipath, medir, func_dicoms_path))
    end

%% realign and reslice all volumes using estimates from 1st echo

% first realign estimate using echo 1 images
% realign uses all session's data

n_e1 = dir(fullfile(medir, sprintf('%s_rest*_e1.nii',curr_rest)));
for i = 1:nscans
    filename{rctr}{i,1} = sprintf('%s,%01d', fullfile(medir, n_e1(1).name),i);
end

end % end of all sessions


%% realign uses all session's data
spm fmri
matlabbatch = [];
matlabbatch{1}.spm.spatial.realign.estimate.data = filename;
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 0; % register to first
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';

% save and run job
fname = fullfile(jobdir, 'me_realignment_e1_all_sessions.mat');
save(fname, 'matlabbatch');
spm_jobman('run', matlabbatch);

% move the spm ps file to current nii directory
spm_ps_file = dir(fullfile(pwd,'spm_*.ps'));
if ~isempty(spm_ps_file)
    movefile(fullfile(pwd,spm_ps_file.name),niidir)
end

%%
for r = 1:n_rest_names % to loop over rest_names
    curr_rest = rest_names{r};
    medir = fullfile(niidir, sprintf('%s_me',curr_rest)); % path to me data
    n_e1 = dir(fullfile(medir, sprintf('%s_rest*_e1.nii',curr_rest)));
    
    fprintf('%s of %s\n',curr_rest,SubID)
    
    % this session has fewer volumes due to scanning cutoff
    if strcmp(SubID,'NODEAP_44') && strcmp(curr_rest,'S1D1')
        nscans = 205;
    else
        nscans = 310;
    end
    
    functdir = fullfile(niidir, 'functional', curr_rest); % path to rsfRMI images
    
    if numel(dir(fullfile(functdir, 'fvol_*.nii')))==nscans % skip current session if the functional dir has all combined scans
        fprintf('functional images have been converted for %s\n',functdir)
        continue;
    end
    
% apply the voxel-to-world mapping from the realigned echo to other echoes by updating header files
for i = 1:nscans
    V = spm_get_space(sprintf('%s,%01d', fullfile(medir, n_e1(1).name),i));
    for iecho = 2:nechos
        n = dir(fullfile(medir, sprintf('%s_rest*%d.nii',curr_rest,iecho))); 
        spm_get_space(sprintf('%s,%01d', fullfile(medir, n(1).name),i),V);
    end
end

%%
% reslice volumes in all echoes
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
matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';

% save and run job
fname = fullfile(jobdir, sprintf('me_reslice_all_echoes_%s.mat',curr_rest));
save(fname, 'matlabbatch');
spm_jobman('run', matlabbatch);

    
%% COMBINE ECHOS
% compute tSNR maps per echo
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

toc