
%% ========================================================================
%  Script: motion_posthoc_ctbs_vs_sham.m
%
%  Purpose
%    Post-hoc analysis of motion (FD) to compare cTBS vs sham sessions.
%    - Loads per-session realignment params (rp_*.txt)
%    - Computes framewise displacement (FD) per volume and sums across vols
%    - Reads bad-volume counts saved earlier (BadVolumes.mat)
%    - Builds a tidy table with TMS labels derived from StimOrder
%    - Runs mixed-effects model (motion ~ TMS_types + StimLoc + (1|SubID))
%    - Reports means ± SD for cTBS/sham × Anterior/Posterior
%
%  Assumptions
%    • rp_*.txt in:  <SubID>/nifti/<Sess>_me/
%    • BadVolumes.mat in:  NRegressor/<SubID>/    (may be absent)
%    • SubInfo in: beh_data_found/SubConds.xlsx (SubID, StimLoc, StimOrder)
%    • count table: experiment_metadata/MRI_func_count.xlsx
%       - columns correspond to 7 sessions (may include an ID column)
%
%  Output variables (in workspace)
%    all_fd             [nSub × 7] summed FD per session
%    all_ctr_badvol     [nSub × 7] bad volume counts
%    df_motion          tidy table of subject×session rows with TMS labels
%    stats / lme / ANOVA for reporting
% ========================================================================

clc; clear;

%% ----------------------------- Setup ------------------------------------
dat_folder = '/Volumes/X9Pro/NODEAP/MRI';
NR_dir     = '/Volumes/X9Pro/NODEAP/NRegressor';
subs_meta  = '/Users/liuq13/project-nodeap-core/beh_data_found/SubConds.xlsx';
count_xlsx = '/Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx';

rest_names   = {'D0','S1D1','S1D2','S2D1','S2D2','S3D1','S3D2'};
n_rest_names = numel(rest_names);

% Subject list (disk order)
SubIDlist = dir(fullfile(dat_folder, 'NODEAP*'));
SubIDlist = SubIDlist([SubIDlist.isdir]);
SubIDs    = string({SubIDlist.name})';
nSub      = numel(SubIDs);

% Availability table (rows must align with SubID order above or contain ID col)
count_table_raw = readmatrix(count_xlsx);
% Heuristic: if table has >=8 cols, assume col 1 is an ID and sessions are cols 2..8
if size(count_table_raw,2) >= 8
    count_table = count_table_raw(:, 2:8);
else
    count_table = count_table_raw(:, 1:7);
end

% Subject meta (SubID, StimLoc, StimOrder)
SubInfo = readtable(subs_meta, 'TextType','string');
% Map SubIDlist order to SubInfo order
idxMap  = zeros(nSub,1);
for i = 1:nSub
    hit = find(SubInfo.SubID==SubIDs(i), 1);
    if isempty(hit), error('SubID %s not found in SubInfo.', SubIDs(i)); end
    idxMap(i) = hit;
end
StimLocSub    = SubInfo.StimLoc(idxMap);     % aligned to SubIDlist order
StimOrderSub  = SubInfo.StimOrder(idxMap);   % numeric codes (123, 231, ...)

%% ------- Collect FD (summed) and bad-volume counts per session ----------
all_ctr_badvol = nan(nSub, n_rest_names);   % # bad vols
all_fd         = nan(nSub, n_rest_names);   % summed FD

head_radius = 50;                            % mm; for radians→mm conversion

for s = 1:nSub
    SubID  = SubIDs(s);
    niidir = fullfile(dat_folder, SubID, 'nifti');

    % Load saved bad-volume info if present
    bv_file = fullfile(NR_dir, SubID, 'BadVolumes.mat');
    if exist(bv_file,'file')
        S = load(bv_file);                   % expects: ctr_badvol
        if isfield(S,'ctr_badvol')
            all_ctr_badvol(s, :) = S.ctr_badvol(:)';
        end
    end

    % Loop sessions
    for r = 1:n_rest_names
        % Skip if count table marks session missing
        val = count_table(s, r);
        if val == 0,   continue; end
        % nscans not strictly needed for FD sum, but kept for sanity
        nscans = (val==0.5) * 205 + (val==1) * 310;

        funcRunDir = fullfile(niidir, sprintf('%s_me', rest_names{r}));
        mpfilename = dir(fullfile(funcRunDir, 'rp_*.txt'));  % SPM realignment params
        if isempty(mpfilename)
            warning('No rp_*.txt for %s %s; skipping.', SubID, rest_names{r});
            continue;
        end

        mp = load(fullfile(funcRunDir, mpfilename(1).name)); % [vols × 6]
        if size(mp,2) ~= 6
            warning('rp_*.txt in %s has %d cols (expected 6).', funcRunDir, size(mp,2));
            continue;
        end

        % FD per Power et al. style: sum(|Δx|+|Δy|+|Δz|) + R*sum(|Δα|+|Δβ|+|Δγ|)
        dmp = [zeros(1,6); diff(mp)];                         % Δ parameters
        translational = sum(abs(dmp(:, 1:3)), 2);
        rotational    = head_radius * sum(abs(dmp(:, 4:6)), 2); % radians→mm
        fd_series     = translational + rotational;

        all_fd(s, r)  = nansum(fd_series);                    % summed FD
        if numel(fd_series) ~= nscans
            warning('Scan count mismatch for %s %s (FD=%d, expected=%d).', ...
                SubID, rest_names{r}, numel(fd_series), nscans);
        end
    end
end

%% ----------------- Build tidy table (subject × session) ------------------
% session-major stacking (D0..S3D2 for subj1, then subj2, ...)
N               = nSub * n_rest_names;
Sess_col        = repmat(categorical(rest_names,rest_names,rest_names).', nSub, 1);
SubID_col       = repmat(SubIDs, n_rest_names, 1);
StimLoc_col     = repmat(StimLocSub, n_rest_names, 1);
StimOrder_col   = repmat(StimOrderSub, n_rest_names, 1);
motion_col      = reshape(all_fd.', [], 1);                   % summed FD

% TMS_types by StimOrder (per subject → expand to 7 sessions)
TMS_types_col = strings(N,1);
row = 0;
for s = 1:nSub
    patt = tms_types_from_order(StimOrderSub(s));             % 1×7 cellstr
    TMS_types_col(row+(1:n_rest_names)) = string(patt(:));
    row = row + n_rest_names;
end
TMS_types_col = categorical(TMS_types_col, {'cTBS','null','sham'});

% time 1..7
time_col = repmat((1:n_rest_names).', nSub, 1);

df_motion = table(SubID_col, Sess_col, motion_col, StimLoc_col, StimOrder_col, ...
                  TMS_types_col, time_col, ...
    'VariableNames', {'SubID','Sess','motion','StimLoc','StimOrder','TMS_types','time'});

%% --------- Filter to cTBS & sham; fit mixed-effects; summarize ----------
use_df = df_motion(df_motion.TMS_types=='cTBS' | df_motion.TMS_types=='sham', :);
use_df.StimLoc   = categorical(use_df.StimLoc, {'Anterior','Posterior'});
use_df.TMS_types = removecats(use_df.TMS_types);
use_df.SubID     = categorical(use_df.SubID);

% Mixed-effects model: random intercept for subject
lme  = fitlme(use_df, 'motion ~ TMS_types + StimLoc + (1|SubID)');
alme = anova(lme);               % ANOVA table for fixed effects
disp(alme);

% Means ± SD (cTBS/sham × Anterior/Posterior)
stats = groupsummary(use_df, {'TMS_types','StimLoc'}, {'mean','std'}, 'motion');
getter = @(t,tms,loc,field) t.(field)(t.TMS_types==tms & t.StimLoc==loc);

c_a_mean = getter(stats,'cTBS','Anterior','mean_motion');
c_a_sd   = getter(stats,'cTBS','Anterior','std_motion');
c_p_mean = getter(stats,'cTBS','Posterior','mean_motion');
c_p_sd   = getter(stats,'cTBS','Posterior','std_motion');

s_a_mean = getter(stats,'sham','Anterior','mean_motion');
s_a_sd   = getter(stats,'sham','Anterior','std_motion');
s_p_mean = getter(stats,'sham','Posterior','mean_motion');
s_p_sd   = getter(stats,'sham','Posterior','std_motion');

fprintf('\nSummary (FD = summed per session):\n');
fprintf('cTBS  — Anterior: %.1f (±%.1f), Posterior: %.1f (±%.1f)\n', ...
    c_a_mean, c_a_sd, c_p_mean, c_p_sd);
fprintf('Sham  — Anterior: %.1f (±%.1f), Posterior: %.1f (±%.1f)\n\n', ...
    s_a_mean, s_a_sd, s_p_mean, s_p_sd);

%% -------------------- Helper: StimOrder → TMS labels --------------------
function types = tms_types_from_order(code)
% Returns a 1×7 cell array (D0..S3D2) using {'null','cTBS','sham'}
    switch code
        case 123, types = {'null','cTBS','sham','sham','cTBS','sham','sham'};
        case 132, types = {'null','cTBS','sham','sham','sham','sham','cTBS'};
        case 213, types = {'null','sham','cTBS','cTBS','sham','sham','sham'};
        case 231, types = {'null','sham','cTBS','sham','sham','cTBS','sham'};
        case 312, types = {'null','sham','sham','cTBS','sham','sham','cTBS'};
        case 321, types = {'null','sham','sham','sham','cTBS','cTBS','sham'};
        otherwise, types = repmat({'null'},1,7);
    end
end
