
% delete some wrong files

clear; clc;
HomeDir = '/Volumes/X9Pro/NODEAP/FuncConn_PAID';
SubIDlist = dir(fullfile(HomeDir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

%%

for subj = 1:nSubIDlist

SubID = SubIDlist(subj).name;
SubDir = fullfile(HomeDir,SubID,'Contrasts');

if exist(SubDir,'dir')
    rmdir(SubDir,'s')
end

% delete all files with the same names
filePattern = fullfile(HomeDir,SubID, '*', '*', 'w2beta_0001.nii');
fileList = dir(filePattern);

disp('Files to be deleted:');
for i = 1:length(fileList)
    disp(fullfile(fileList(i).folder, fileList(i).name));
end

% Delete each target file
for i = 1:length(fileList)
    delete(fullfile(fileList(i).folder, fileList(i).name));
    fprintf('Deleted file: %s\n', fullfile(fileList(i).folder, fileList(i).name));
end

end




