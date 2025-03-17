
% delete dicoms.zip file under each subfolder

clear; clc;
HomeDir = '/Volumes/X9Pro/NODEAP/MRI';
SubIDlist = dir(fullfile(HomeDir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

for subj = 1:nSubIDlist

SubID = SubIDlist(subj).name;
SubDir = fullfile(HomeDir,SubID);

fprintf('%s\n',SubDir)

zipfile = fullfile(SubDir,'dicoms.zip');
if exist(zipfile,'file')
    delete(zipfile)
end

niidir = fullfile(SubDir,'nifti');
if exist(niidir,'dir')
    rmdir(niidir,'s')
end

tedanadir = fullfile(SubDir,'tedana_output');
if exist(tedanadir,'dir')
    rmdir(tedanadir,'s')
end

end
