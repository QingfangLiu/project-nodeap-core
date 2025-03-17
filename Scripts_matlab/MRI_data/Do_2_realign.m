

clear; clc;
HomeDir = '/Volumes/X9Pro/NODEAP/MRI';
[~, scriptName, ~] = fileparts(mfilename('fullpath'));

SubIDlist = dir(fullfile(HomeDir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

spm fmri

%%
for subj = 1:nSubIDlist

    SubID = SubIDlist(subj).name;
    SubDir = fullfile(HomeDir,SubID);
    niidir = fullfile(SubDir, 'nifti');
    jobdir = fullfile(SubDir, 'MRIjobs');
    
    output_file = sprintf('output_%s_%s.txt', scriptName, string(datetime));
    diary(output_file); % keep notes of outputs
    
    rctr = 0; % count rest sessions (in case of missing session)
    filename = []; % all e1 scans across all sessions

    for r = 1:n_rest_names % to loop over rest_names
        curr_rest = rest_names{r};
        medir = fullfile(niidir, sprintf('%s_me',curr_rest)); % path to me data

        % this session has fewer volumes due to scanning cutoff
        if strcmp(SubID,'NODEAP_44') && strcmp(curr_rest,'S1D1')
            nscans = 205;
        else
            nscans = 310;
        end
    
        if ~exist(medir,'dir')
            continue;
        end

        rctr = rctr + 1; % to account for potential missing sessions
        
        % organize echo 1 images across sessions
        n_e1 = dir(fullfile(medir, sprintf('%s_rest*_e1.nii',curr_rest)));
        for i = 1:nscans
            filename{rctr}{i,1} = sprintf('%s,%01d', fullfile(medir, n_e1(1).name),i);
        end

    end
    
%% realign using e1 from all session's data

disp('Realign using e1 data from all sessions.')
matlabbatch = [];

if strcmp(SubID,'NODEAP_41') % for NODEAP_41: S3D2 session is realigned separately
    matlabbatch{1}.spm.spatial.realign.estimate.data = filename(1:6);
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 0; % register to first to be faster
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';

    matlabbatch{2}.spm.spatial.realign.estimate.data = filename(7);
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.sep = 4;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.rtm = 0; % register to first to be faster
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.interp = 2;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.weight = '';
else
    matlabbatch{1}.spm.spatial.realign.estimate.data = filename;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 0; % register to first to be faster
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';
end

% save and run job
fname = fullfile(jobdir, 'me_realignment_e1_all_sessions.mat');
save(fname, 'matlabbatch');
spm_jobman('run', matlabbatch);

% move the spm ps file to current nii directory & convert to pdf
spm_ps_file = dir(fullfile(pwd,'spm_*.ps'));
if ~isempty(spm_ps_file)
    inputPS = spm_ps_file.name;
    outputPDF = 'spm_realign.pdf';
    ghostscriptPath = '/opt/homebrew/bin/gs'; 
    command = sprintf('%s -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -sOutputFile="%s" "%s"', ghostscriptPath, outputPDF, inputPS);
    system(command);
    movefile(fullfile(pwd,inputPS),niidir)
    movefile(fullfile(pwd,outputPDF),niidir)
end

diary off;
movefile(fullfile(pwd,output_file),jobdir)

end

