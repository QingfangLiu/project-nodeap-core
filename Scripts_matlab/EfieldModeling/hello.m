

clear; clc; close all

userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Library/CloudStorage/Box-Box/NODEAP_data_analysis');
studydir = '/Volumes/X9Pro/NODEAP';
MRIdir = fullfile(studydir,'MRI');
EfieldDir = fullfile(studydir,'EfieldModeling');

% read MRI count file
count_table = xlsread('/Volumes/X9Pro/NODEAP/MRI_func_count.xlsx');

SubIDlist = dir(fullfile(MRIdir, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]); % only keep directories
nSubIDlist = length(SubIDlist);

rest_names = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = length(rest_names);
ROIdir = fullfile(studydir,'ConnectivityMasks');

% add SimNIBS toolbox
addpath('/Users/liuq13/Applications/SimNIBS-4.1/matlab_tools')

% more examples can be found 
% /Users/<USER_NAME>/Applications/SimNIBS.app/examples

StimCoordFile = fullfile(HomeDir,'TMSvisualization','StimCoord_Daria.xlsx');
StimCoord = readtable(StimCoordFile);

%% create a head model
% need T1 image as the input
% Segmentation and meshing
% this takes long time (~1h) to run to find a good solution, segmenting,
% meshing, etc
% can use parallel computing to speed it up (but with weird behaviors)
% or use biowulf to run it
% https://hpc.nih.gov/apps/simnibs.html

charm_path = '/Users/liuq13/Applications/SimNIBS-4.1/bin/charm';  

for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;
    fprintf('%s\n',SubID)
    
    SubDir = fullfile(MRIdir,SubID);
    niidir = fullfile(SubDir, 'nifti');
    OutDir = fullfile(EfieldDir,SubID);
    if ~exist(OutDir,'dir')
        mkdir(OutDir)
        cd(OutDir) % change working directory so new folder will be under this
        aname = dir(fullfile(niidir,'anat','D0_T1*.nii'));
        afile = fullfile(niidir,'anat',aname(1).name);
        command = sprintf('%s %s %s', charm_path, SubID, afile);
        [status, result] = system(command, '-echo');
    end
    
    cd(OutDir) % change working directory so new folder will be under this
    
    %% Starting a SESSION and Selecting a Head Mesh

    s = sim_struct('SESSION');  % Initialize a session
    s.subpath = sprintf('m2m_%s',SubID);  % Name of head mesh
    option = string(StimCoord{subj, 'StimLoc'});
    s.pathfem = sprintf('SimuOut_%s/',option); % put stim target as part of the output folder name
    
    if exist(fullfile(OutDir,s.pathfem),'dir')
        %continue;                            % skip this
        rmdir(fullfile(OutDir,s.pathfem),'s') % delete and redo
    end

    % also to update Output and Post-processing options

    %% Setting up a TMS simulation
    % add a TMS list to the SESSION structure and select a coil model

    s.poslist{1} = sim_struct('TMSLIST'); % Initialize a list of TMS simulations
    s.poslist{1}.fnamecoil = fullfile('Drakaki_BrainStim_2022','MagVenture_Cool-B65.ccd');  % updated from default

    % set a position for the coil
    % should be coordinates in subject space
    use_coords = table2array(StimCoord(subj,[3,4,5]));
    use_coords_direction = use_coords;
    use_coords_direction(2) = use_coords_direction(2) + 10;
    % how to define the direction here?

    s.poslist{1}.pos(1).centre = use_coords;
    s.poslist{1}.pos(1).pos_ydir = use_coords_direction;
    s.poslist{1}.pos(1).distance = 4; % 4 mm distance from coil surface to head surface
    % s.poslist{1}.pos(1).didt = 20e6; % define dI/dt=20e6 A/s
    % here we can also set stimulator intensity, how?

    tic
    try
        run_simnibs(s)
    catch
        warning('Problem using function.');    
    end
    toc

end

%% Visualize Simulations 
magnE_all = nan(nSubIDlist,2);

for subj = 1:nSubIDlist
    SubID = SubIDlist(subj).name;
    fprintf('%s\n',SubID)
    OutDir = fullfile(EfieldDir,SubID);
    option = string(StimCoord{subj, 'StimLoc'});
        
    % load mesh
    msh_name = dir(fullfile(OutDir,sprintf('SimuOut_%s/',option),'*.msh'));
    
    if isempty(msh_name)
        continue;
    end
    
    msh_file = fullfile(msh_name(1).folder,msh_name(1).name);
    m = mesh_load_gmsh4(msh_file);
    
%     % display some key results for whole cortex
%     disp(' ')
%     disp('whole cortex:')
%     summary=mesh_get_fieldpeaks_and_focality(m,'field_idx',2);
% 
%     % show field on the GM surface
%     mesh_show_surface(m);
    
    
    % -------------------------------------------------------
% EXAMPLE 2
% extract a spherical ROI with 10 mm radius around the peak position
%
% the "summary" structure contains the peak values
% together with their positions:
%   "summary.percentiles" lists the tested percentile cutoffs - 
%   the 99.9 percentile is the 3rd entry
%
%   "summary.perc_values" lists the corresponding values, as they
%   are also displayed by mesh_get_fieldpeaks_and_focality
%
%   "summary.XYZ_perc" lists the corresponding center positions -
%   the center position for the 99.9 percentile is the 3rd row:
% peak_pos=summary.XYZ_perc(3,:);
% 
% % distance to peak position
% dist=sqrt(sum(bsxfun(@minus,m.nodes,peak_pos).^2,2));
% 
% % extract nodes closer than 10 mm, and belonging to tetrahedra with region
% % number 2 (gray matter has region number 2)
% node_idx=dist<10;
% m_ROI=mesh_extract_regions(m, 'elemtype','tet','region_idx',2,'node_idx', node_idx);
% 
% % get some key results for the spherical ROI
% disp(' ')
% disp('10 mm spherical ROI around peak position:')
% mesh_get_fieldpeaks_and_focality(m_ROI,'field_idx',2);
% 
% % show the extracted ROI:
% % show the whole cortex semi-transparent
% mesh_show_surface(m,'showSurface',true,'facealpha',0.4);
% % create new mesh containing the enclosing surface of the ROI tetrahedra
% TR = triangulation(m_ROI.tetrahedra,m_ROI.nodes(:,1),m_ROI.nodes(:,2),m_ROI.nodes(:,3));
% m_vis=mesh_empty;
% m_vis.nodes=m_ROI.nodes;
% m_vis.triangles=freeBoundary(TR);
% m_vis.triangle_regions=ones(size(m_vis.triangles,1),1);
% % add the he enclosing surface to the plot
% mesh_show_surface(m_vis,'showSurface',true,'surfaceColor',[1 0 0],'haxis',gca); 
% title('10 mm spherical ROI around peak position');


%% ROI analysis

% Crop the mesh so we only have gray matter volume elements (tag 2 in the mesh)
gray_matter = mesh_extract_regions(m, 'region_idx', 2);

% define ROI
use_dir = fullfile(OutDir,sprintf('m2m_%s',SubID));
use_coords_aOFC = mni2subject_coords([34, 54, -14], use_dir);
use_coords_pOFC = mni2subject_coords([28, 38, -16], use_dir);
use_coords_both = [use_coords_aOFC;use_coords_pOFC];

for k = 1:2

% we will use a sphere of radius 10 mm
r = 10;

% Electric fields are defined in the center of the elements
% get element centers
elm_centers = mesh_get_tetrahedron_centers(gray_matter);
% determine the elements in the ROI
roi = sqrt(sum(bsxfun(@minus, elm_centers, use_coords_both(k,:)).^2, 2)) < r;
% get element volumes, we will use those for averaging
elm_vols = mesh_get_tetrahedron_sizes(gray_matter);

%% Get field and calculate the mean
% Get the field of interest
field_name = 'magnE';
field_idx = get_field_idx(gray_matter, field_name, 'elements');
field = gray_matter.element_data{field_idx}.tetdata;

% Calculate the mean
avg_field_roi = sum(field(roi) .* elm_vols(roi))/sum(elm_vols(roi));
fprintf('mean %s in ROI: %f\n', field_name, avg_field_roi)

magnE_all(subj,k) = avg_field_roi;

end


end
