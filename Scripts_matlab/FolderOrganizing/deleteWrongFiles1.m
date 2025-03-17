
% delete some wrong files

clear; clc;
HomeDir = '/Volumes/X9Pro/NODEAP/MRI';
SubIDlist = dir(fullfile(HomeDir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

for subj = 1:nSubIDlist

SubID = SubIDlist(subj).name;
SubDir = fullfile(HomeDir,SubID);

for r = 1:n_rest_names % to loop over rest_names
    
    curr_rest = rest_names{r};
        
    tmp1 = fullfile(SubDir,'nii','functional',curr_rest,'w3fvol_4d.nii');
    tmp2 = fullfile(SubDir,'nii','functional',curr_rest,'s6w3fvol_4d.nii');
    tmp3 = fullfile(SubDir,'nii','functional',curr_rest,'w2fvol_4d.nii');
    tmp4 = fullfile(SubDir,'nii','functional',curr_rest,'s6w2fvol_4d.nii');

if exist(tmp1,'file')
    delete(tmp1)
end

if exist(tmp2,'file')
    delete(tmp2)
end

if exist(tmp3,'file')
    delete(tmp3)
end

if exist(tmp4,'file')
    delete(tmp4)
end


end

end
