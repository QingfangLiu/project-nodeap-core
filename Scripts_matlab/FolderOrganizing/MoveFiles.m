

% This script moves TMS-preprocessing related results to the nii_for_TMS
% folder

clc
clear
close all

External = '/Volumes/X9Pro/NODEAP/MRI';
subdirs = dir(fullfile(External,'NODEAP*'));

%%
tic
parfor j = 1:length(subdirs)
    
    fprintf('Working on %s\n',subdirs(j).name)

suborigindir = fullfile(External,subdirs(j).name,'GLM_aOFC_seed_right');
subtargetdir = fullfile(External,subdirs(j).name,'nii_for_TMS');
if exist(suborigindir,'dir')
     movefile(suborigindir,subtargetdir)
end

suborigindir = fullfile(External,subdirs(j).name,'GLM_pOFC_seed_right');
subtargetdir = fullfile(External,subdirs(j).name,'nii_for_TMS');
if exist(suborigindir,'dir')
     movefile(suborigindir,subtargetdir)
end

suborigindir = fullfile(External,subdirs(j).name,'stimulation_coordinates.mat');
subtargetdir = fullfile(External,subdirs(j).name,'nii_for_TMS');
if exist(suborigindir,'file')
     movefile(suborigindir,subtargetdir)
end

end
toc
