
%% to convert dicom data to nifti (4D)

clear; clc;
HomeDir = '/Volumes/X9Pro/NODEAP/MRI';
OriDir = '/Users/liuq13/NODEAP/MRI'; % the original MRI data (dicom)

[~, scriptName, ~] = fileparts(mfilename('fullpath'));

dcm2niipath = '/Applications/MRIcroGL.app/Contents/Resources/dcm2niix';
TEs = [14.60, 39.04, 63.48];
nechos = 3;
SubIDlist = dir(fullfile(HomeDir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

spm fmri

%%
for subj = 1:nSubIDlist

    tic
    output_file = sprintf('output_%s_%s.txt', scriptName, string(datetime));
    diary(output_file); % keep notes of outputs
    
    SubID = SubIDlist(subj).name;
    SubDir = fullfile(HomeDir,SubID);
    niidir = fullfile(SubDir, 'nifti');
    jobdir = fullfile(SubDir, 'MRIjobs');
    if ~exist(jobdir,'dir') 
       mkdir(jobdir)
    end
    
    %% convert dicom to nifti files
    % convert T1 file
    anatpath = fullfile(OriDir, SubID, 'dicoms', 'D0_T1');
    outdir = fullfile(niidir, 'anat');
    if ~exist(outdir,'dir')
        mkdir(outdir);
        system(sprintf('%s -o %s %s', dcm2niipath, outdir, anatpath))
    end

%%
% convert resting-state data
rctr = 0; % count rest sessions (in case of missing session)
filename = []; % all e1 scans across all sessions

for r = 1:n_rest_names % to loop over rest_names
    curr_rest = rest_names{r};
    fprintf('%s of %s\n',curr_rest,SubID)
    disp('convert dicom to nifti files')
    
    func_dicoms_path = fullfile(OriDir, SubID, 'dicoms',sprintf('%s_rest',curr_rest));
    medir = fullfile(niidir, sprintf('%s_me',curr_rest)); % path to me data
    
    % this session has fewer volumes due to scanning cutoff
    if strcmp(SubID,'NODEAP_44') && strcmp(curr_rest,'S1D1')
        nscans = 205;
    else
        nscans = 310;
    end
    
    % skip current session if no dicom data
    if ~exist(func_dicoms_path,'dir') || numel(dir(func_dicoms_path)) <=2
        fprintf('No such dicom data %s\n',func_dicoms_path)
        continue;
    end
    
    rctr = rctr + 1; % to account for potential missing sessions
    
    if ~exist(medir,'dir')
        mkdir(medir);
        system(sprintf('%s -o %s %s', dcm2niipath, medir, func_dicoms_path))
    end
    
    % organize echo 1 images across sessions
    n_e1 = dir(fullfile(medir, sprintf('%s_rest*_e1.nii',curr_rest)));
    for i = 1:nscans
        filename{rctr}{i,1} = sprintf('%s,%01d', fullfile(medir, n_e1(1).name),i);
    end

end % end of all sessions

diary off;
movefile(fullfile(pwd,output_file),jobdir)
toc

end


