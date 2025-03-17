
% This is copied from TK's code to recreate the TMS folders

clear all
close all

% what do you want to do
do_dicom_import = 1;
do_bias_correction = 1;     % save bias corrected anatomical image for TMS naviation

% add path to dcm2niix
dcm2niipath = '/Applications/MRIcroGL.app/Contents/Resources/dcm2niix';

studydir = '/Volumes/X9Pro/NODEAP/MRI';
SubIDlist = dir(fullfile(studydir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

for subj = 2:nSubIDlist

    sn = SubIDlist(subj).name;

% get data path
subdir = fullfile(studydir,sn);

%%% DICOM IMPORT 
if do_dicom_import
    matlabbatch = [];
    
    % all dicoms must be in their subfolders within 'dicoms', e.g.,.../NODEAP/NODEAP_pilot01/dicoms
    anatpath = fullfile(subdir, 'dicoms', 'D0_T1');
    functpath = fullfile(subdir, 'dicoms', 'D0_rest');

    % convert T1 file
    outdir = fullfile(subdir,'nii_for_TMS', 'anat');
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    system(sprintf('%s -o %s %s', dcm2niipath, outdir, anatpath))
    
    
    % convert resting-state data
    outdir = fullfile(subdir,'nii_for_TMS', 'day0_me');
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    system(sprintf('%s -o %s %s', dcm2niipath, outdir, functpath))
    
end

if do_bias_correction
    matlabbatch = [];
    
    % set filenames and paths
    subfolder = fullfile(subdir, 'nii_for_TMS', 'anat'); % path with anatomical data
    anat = dir(fullfile(subfolder,'D0_T1*.nii'));
    
    % change folder to anatomy folder
    cd(subfolder)
    
    fname = fullfile(subfolder, anat(1).name);
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {fname};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
    ent = [1, 1, 2, 3, 4, 2];
    for i = 1:6
        matlabbatch{1}.spm.spatial.preproc.tissue(i).tpm = {fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%01d', i))};
        matlabbatch{1}.spm.spatial.preproc.tissue(i).ngaus = ent(i);
        matlabbatch{1}.spm.spatial.preproc.tissue(i).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(i).warped = [0 0];
    end
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
    
    spm_jobman('run', matlabbatch);
    
    for i=1:6
        delete(fullfile(subfolder, sprintf('c%01d%s', i, anat(1).name)))
    end
    
end

end


