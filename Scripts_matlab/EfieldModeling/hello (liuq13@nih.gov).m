

clear; clc; close all

userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Library/CloudStorage/Box-Box/NODEAP_data_analysis');

studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');

SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);

ROIdir = fullfile(studydir,'ConnectivityMasks');


% add SimNIBS toolbox
addpath('/Users/liuq13/Applications/SimNIBS-4.1/matlab_tools')

EfieldDir = fullfile(studydir,'EfieldModeling');

% more examples can be found 
% /Users/<USER_NAME>/Applications/SimNIBS.app/examples


%% create a head model
% need T1 image as the input
% Segmentation and meshing
% this takes long time (~1h) to run to find a good solution, segmenting,
% meshing, etc
% can use parallel computing to speed it up
% or use biowulf to run it
% https://hpc.nih.gov/apps/simnibs.html

for subj = 2 %:nSubIDlist

%subj = 1;
SubID = SubIDlist(subj).name;
SubDir = fullfile(MRIdir,SubID);
niidir = fullfile(SubDir, 'nifti');
OutDir = fullfile(EfieldDir,SubID);
% if exist(OutDir,'dir')
%     continue;
% end
mkdir(OutDir)
cd(OutDir) % change working directory so new folder will be under this

aname = dir(fullfile(niidir,'anat','D0_T1*.nii'));
afile = fullfile(niidir,'anat',aname(1).name);
charm_path = '/Users/liuq13/Applications/SimNIBS-4.1/bin/charm';  
command = sprintf('%s %s %s', charm_path, SubID, afile);
[status, result] = system(command, '-echo');

end


%% Starting a SESSION and Selecting a Head Mesh

s = sim_struct('SESSION');  % Initialize a session
s.subpath = sprintf('m2m_%s',SubID);  % Name of head mesh
option = 'aOFC';
s.pathfem = sprintf('SimuOut_%s/',option); % put stim target as part of the output folder name

% also to update Output and Post-processing options

%% Setting up a TMS simulation
% add a TMS list to the SESSION structure and select a coil model

s.poslist{1} = sim_struct('TMSLIST'); % Initialize a list of TMS simulations
s.poslist{1}.fnamecoil = fullfile('legacy_and_other','Magstim_70mm_Fig8.ccd');  % Select coil

% set a position for the coil
% should be coordinates in subject space
s.poslist{1}.pos(1).centre = mni2subject_coords([42, 52, 35], s.subpath);
% Point the coil handle posteriorly, we just add 10 mm to the original M1 "y" coordinate
s.poslist{1}.pos(1).pos_ydir = mni2subject_coords([42+10, 52, 35], s.subpath);
% how to define the direction here?
s.poslist{1}.pos(1).distance = 4; % 4 mm distance from coil surface to head surface
% s.poslist{1}.pos(1).didt = 20e6; % define dI/dt=20e6 A/s

% Do I use target coordinates or entry coordinates (in localite language)?
% here we can also set stimulator intensity, how?

%% Run Simulations
tic
run_simnibs(s)
toc
% 
% %% Visualize Simulations 
% msh_name = dir(fullfile(s.pathfem,'*.msh'));
% msh_file = fullfile(s.pathfem,msh_name(1).name);
% m = mesh_load_gmsh4(msh_file);
% mesh_show_surface(m);


