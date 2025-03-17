

% This script copies all dicomes folder to external drive

clc
clear
close all

studydir = '/Users/liuq13/NODEAP/MRI';
External = '/Volumes/QL1T/NODEAP/MRI';
subdirs = dir(fullfile(studydir,'NODEAP*'));

%%
tic
parfor j = 1:length(subdirs)
    
    fprintf('Working on %s\n',subdirs(j).name)

suborigindir = fullfile(studydir,subdirs(j).name,'dicoms');
subtargetdir = fullfile(External,subdirs(j).name,'dicoms');
if ~exist(subtargetdir,'dir')
     mkdir(subtargetdir)
     copyfile(suborigindir,subtargetdir)
end
end
toc
