

% This script moves TMS-preprocessing related results to the nii_for_TMS
% folder

clc
clear
close all

External = '/Volumes/X9Pro/NODEAP';
subdirs = dir(fullfile(External,'MRI','NODEAP*'));

%%
tic
parfor j = 1:length(subdirs)
    
    fprintf('Working on %s\n',subdirs(j).name)

    suborigindir = fullfile(External,'MRI',subdirs(j).name,'FuncConn_tedana1');
    
    
    
    subtargetdir = fullfile(External,'FuncConn_tedana1',subdirs(j).name);
    if exist(suborigindir,'dir')
         movefile(suborigindir,subtargetdir)
    end

end
toc
