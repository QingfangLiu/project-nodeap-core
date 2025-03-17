
clear

rootDir = '/Users/liuq13/NODEAP';

userHome = getenv('HOME');
HomeDir = fullfile(userHome,'Library/CloudStorage/Box-Box/NODEAP_data_analysis');

file_path = fullfile(HomeDir,'ProcessedData','SubConds.xlsx');
table_data = readtable(file_path); 

Subs = table_data.SubID;
Sweet_Savory = table_data.StartOdor;
Order = table_data.StimOrder;
StimLoc = table_data.StimLoc;
nSubs = length(Subs);

AllData = [];
AllChoiceData = [];
AllConditioningData = [];

if_process_pleasant_data = 0;
if_process_choice_data = 1;
if_process_conditioning_data = 0;

% this code only deals with NODEAP_17 because its 2nd session did not
% complete conditioning task and the data were saved badly

% find which subject is NODEAP_17
j = find(strcmp(Subs, 'NODEAP_17'));

%%
%for j = 1:nSubs
SubID = Subs{j};
SubDir = fullfile(rootDir,'NODEAP_Behavioral_Data_Reupload',SubID);
filename = fullfile(SubDir,sprintf('complete_DAY2_POST_%s_SESS_03.mat',SubID));
load(filename)

% disp(res.odors)
% confirmed that 1st odor is always sweet (1,2,3,4) & 2nd odor is savory
% (5,6,7,8)

if if_process_pleasant_data

%% pleasantness rating
% organize this into a df
% col1: session number in time order (1,2,3)
% col2: pre (1) or post (2)
% col3: odor index 1 or 2
% col4: pleasantness rating (averaged from 3 times)
% col5: specific odor index (1-8)

% pre pleasantness rating
all_df = [];
for i = [1,3]
    filename = fullfile(SubDir,sprintf('pleasantness_ratings_PRE_%s_SESS_%02d.mat',SubID,i));
    load(filename)
    x = pleasantness_ratings_PRE;
    idx1 = find(x(:,1)==1);
    idx2 = find(x(:,1)==2);
    mean_ratings = [[1;2],[mean(x(idx1,5));mean(x(idx2,5))]];
    df = [repmat(i,2,1),repmat(1,2,1),mean_ratings];
    all_df = [all_df;df];
end

% post pleasantness rating
for i = [1,3]
    filename = fullfile(SubDir,sprintf('pleasantness_ratings_POST_%s_SESS_%02d.mat',SubID,i));
    load(filename)
    x = pleasantness_ratings_POST;
    idx1 = find(x(:,1)==1);
    idx2 = find(x(:,1)==2);
    mean_ratings = [[1;2],[mean(x(idx1,5));mean(x(idx2,5))]];
    df = [repmat(i,2,1),repmat(2,2,1),mean_ratings];
    all_df = [all_df;df];
end

%% 
varNames = {'Sess', 'PrePost', 'Odor', 'Pleasant'};
Data = array2table(all_df, 'VariableNames', varNames);

Data.OdorType(Data.Odor == 1) = 'W';
Data.OdorType(Data.Odor == 2) = 'A';

% add subject ID
Data.SubID = repmat(SubID, size(all_df,1), 1);
% add StimLoc
Data.StimLoc = repmat(StimLoc(j), size(all_df,1), 1);
% add stim order
Data.StimOrder = repmat(Order(j), size(all_df,1), 1);

% add a Devalued column: which odor was devalued at current session
if strcmp(Sweet_Savory(j),'Sweet')
    Data.Devalued(Data.Sess==1) = 'W';
    Data.Devalued(Data.Sess==2) = 'A';
    Data.Devalued(Data.Sess==3) = 'W';
else
    Data.Devalued(Data.Sess==1) = 'A';
    Data.Devalued(Data.Sess==2) = 'W';
    Data.Devalued(Data.Sess==3) = 'A';
end

% add a Cond column: X (cTBS-sham), Y (sham-cTBS), z (sham-sham)
if Order(j)==123
    Data.Cond(Data.Sess==1) = 'X';
    Data.Cond(Data.Sess==2) = 'Y';
    Data.Cond(Data.Sess==3) = 'Z';
    
    elseif Order(j)==132
    Data.Cond(Data.Sess==1) = 'X';
    Data.Cond(Data.Sess==2) = 'Z';
    Data.Cond(Data.Sess==3) = 'Y';
    
    elseif Order(j)==213
    Data.Cond(Data.Sess==1) = 'Y';
    Data.Cond(Data.Sess==2) = 'X';
    Data.Cond(Data.Sess==3) = 'Z';
    
    elseif Order(j)==231
    Data.Cond(Data.Sess==1) = 'Y';
    Data.Cond(Data.Sess==2) = 'Z';
    Data.Cond(Data.Sess==3) = 'X';
    
    elseif Order(j)==312
    Data.Cond(Data.Sess==1) = 'Z';
    Data.Cond(Data.Sess==2) = 'X';
    Data.Cond(Data.Sess==3) = 'Y';
    
    elseif Order(j)==321
    Data.Cond(Data.Sess==1) = 'Z';
    Data.Cond(Data.Sess==2) = 'Y';
    Data.Cond(Data.Sess==3) = 'X';
end

% if this odor was devalued at this session
Data.IfDevalue(Data.Devalued == Data.OdorType) = 1;

AllData = [AllData;Data];
clear Data

end

if if_process_choice_data

%% choice data
% cols: 
    % sess (in time order): 1,2,3
    % pre/post: 1 pre, 2 post
    % trial: manually added for each choice task (1-30 or 1-60)
    
    % 1     cue pair (1:16; 1:3 = SW/0 set A, 4:6 = SA/0 set A, 7:15 = SW/SA set A)
    % 2     choice type: (1=sweet; 2=savory; 3=both)
    % 3     odor set: (1=A, 2=B (B is not shown on DAY2 at PRE))
    % 4     left/right (1=odor cue on left, i.e. not flipped, 2=flipped l/r)
    % 6     odor chosen: (1=yes, 0=no for S+/-; 1=sweet,2=savory for S+/+)
    % 8     rt: from onset of choice
    % 11    cue # on the left
    % 12    cue # on the right
    % 13    selected cue #

all_df = [];
% 30 trials per session in pre-trials
for i = [1,3]
    tmp_mat_name = sprintf('choice_PRE_%s_SESS_%02d.mat',SubID,i);
    tmp_mat_file = fullfile(SubDir,tmp_mat_name);
    load(tmp_mat_file)
    x = choice_PRE;
    df = [repmat(i,30,1),repmat(1,30,1),(1:30)',x(:,[1,2,3,4,6,8,11,12,13])];
    all_df = [all_df;df];
end

% 60 trials per session in post-trials
for i = [1,3]
    tmp_mat_name = sprintf('choice_POST_%s_SESS_%02d.mat',SubID,i);
    tmp_mat_file = fullfile(SubDir,tmp_mat_name);
    load(tmp_mat_file)
    x = choice_POST;
    df = [repmat(i,60,1),repmat(2,60,1),(1:60)',x(:,[1,2,3,4,6,8,11,12,13])];
    all_df = [all_df;df];
end

varNames = {'Sess', 'PrePost', 'Trial', 'CuePair', 'ChoiceType', 'Set', 'OdorLR', 'ChosenOdor','rt','CueLeft','CueRight','SelectedCue'};
Data = array2table(all_df, 'VariableNames', varNames);

% add subject ID
Data.SubID = repmat(SubID, size(all_df,1), 1);
% add StimLoc
Data.StimLoc = repmat(StimLoc(j), size(all_df,1), 1);
% add stim order
Data.StimOrder = repmat(Order(j), size(all_df,1), 1);

if strcmp(Sweet_Savory(j),'Sweet')
    Data.Devalued(Data.Sess==1) = 'W';
    Data.Devalued(Data.Sess==2) = 'A';
    Data.Devalued(Data.Sess==3) = 'W';
else
    Data.Devalued(Data.Sess==1) = 'A';
    Data.Devalued(Data.Sess==2) = 'W';
    Data.Devalued(Data.Sess==3) = 'A';
end

% add a Cond column
if Order(j)==123
    Data.Cond(Data.Sess==1) = 'X';
    Data.Cond(Data.Sess==2) = 'Y';
    Data.Cond(Data.Sess==3) = 'Z';
    
    elseif Order(j)==132
    Data.Cond(Data.Sess==1) = 'X';
    Data.Cond(Data.Sess==2) = 'Z';
    Data.Cond(Data.Sess==3) = 'Y';
    
    elseif Order(j)==213
    Data.Cond(Data.Sess==1) = 'Y';
    Data.Cond(Data.Sess==2) = 'X';
    Data.Cond(Data.Sess==3) = 'Z';
    
    elseif Order(j)==231
    Data.Cond(Data.Sess==1) = 'Y';
    Data.Cond(Data.Sess==2) = 'Z';
    Data.Cond(Data.Sess==3) = 'X';
    
    elseif Order(j)==312
    Data.Cond(Data.Sess==1) = 'Z';
    Data.Cond(Data.Sess==2) = 'X';
    Data.Cond(Data.Sess==3) = 'Y';
    
    elseif Order(j)==321
    Data.Cond(Data.Sess==1) = 'Z';
    Data.Cond(Data.Sess==2) = 'Y';
    Data.Cond(Data.Sess==3) = 'X';
end

AllChoiceData = [AllChoiceData;Data];


end


if if_process_conditioning_data


%% conditioning data
% from Day 1
% cols: 
    % sess (in time order): 1,2,3
    % run: 1-5
    % trial: manually added for each choice task (1-24)
    % 1     cue: pair (1:12; 1:3 = SW/0 A, 4:6 = SW/0 set B, 7:9 = SA/0 set A, 10:12 = SA/0 set B)
    % 2     choice type: (1=sweet; 2=savory)
    % 3     odor set: (1=A, 2=B)
    % 4     left/right (1=odor cue on left, i.e. not flipped, 2=flipped l/r)';...
    % 5     choice (1=left, 2=right)';...
    % 6     odor chosen: (1=yes, 0=no)
    % 8     rt from onset of choice (3s max)
    % 11    cue # on the left
    % 12    cue # on the right
    % 13    selected cue #
    
% here I first tried including partial data from 2nd run, but apparently
% conditioning was far from complete so updated this to only include 1st
% and 3rd sessions
    
all_df = [];
for i = [1,3]
    for k = 1:5 % 5 runs in conditioning
        % conditioning data are also saved in a 'conditioning' variable
        tmp_mat_name = sprintf('conditioning_%s_SESS_%02d_run_%02d.mat',SubID,i,k);
        tmp_mat_file = fullfile(SubDir,tmp_mat_name);
        if exist(tmp_mat_file,'file')
            load(tmp_mat_file)
            x = conditioning;
            df = [repmat(i,24,1),repmat(k,24,1),(1:24)',x(:,[1,2,3,4,5,6,8,11,12,13])];
            all_df = [all_df;df];
        end
    end
end

varNames = {'Sess', 'Run', 'Trial','CuePair','ChoiceType', 'Set','OdorLR','ChoiceLR','OdorChosen','rt','CueLeft','CueRight','SelectedCue'};
Data = array2table(all_df, 'VariableNames', varNames);

% add subject ID
Data.SubID = repmat(SubID, size(all_df,1), 1);
% add StimLoc
Data.StimLoc = repmat(StimLoc(j), size(all_df,1), 1);
% add stim order
Data.StimOrder = repmat(Order(j), size(all_df,1), 1);

if strcmp(Sweet_Savory(j),'Sweet')
    Data.Devalued(Data.Sess==1) = 'W';
    Data.Devalued(Data.Sess==2) = 'A';
    Data.Devalued(Data.Sess==3) = 'W';
else
    Data.Devalued(Data.Sess==1) = 'A';
    Data.Devalued(Data.Sess==2) = 'W';
    Data.Devalued(Data.Sess==3) = 'A';
end

% add a Cond column
if Order(j)==123
    Data.Cond(Data.Sess==1) = 'X';
    Data.Cond(Data.Sess==2) = 'Y';
    Data.Cond(Data.Sess==3) = 'Z';
    
    elseif Order(j)==132
    Data.Cond(Data.Sess==1) = 'X';
    Data.Cond(Data.Sess==2) = 'Z';
    Data.Cond(Data.Sess==3) = 'Y';
    
    elseif Order(j)==213
    Data.Cond(Data.Sess==1) = 'Y';
    Data.Cond(Data.Sess==2) = 'X';
    Data.Cond(Data.Sess==3) = 'Z';
    
    elseif Order(j)==231
    Data.Cond(Data.Sess==1) = 'Y';
    Data.Cond(Data.Sess==2) = 'Z';
    Data.Cond(Data.Sess==3) = 'X';
    
    elseif Order(j)==312
    Data.Cond(Data.Sess==1) = 'Z';
    Data.Cond(Data.Sess==2) = 'X';
    Data.Cond(Data.Sess==3) = 'Y';
    
    elseif Order(j)==321
    Data.Cond(Data.Sess==1) = 'Z';
    Data.Cond(Data.Sess==2) = 'Y';
    Data.Cond(Data.Sess==3) = 'X';
end

AllConditioningData = [AllConditioningData;Data];


end

%end

%%
% save table into xlsx
% reorder columns before saving

if if_process_pleasant_data
    AllData = AllData(:, {'SubID','StimLoc','StimOrder','Sess','Cond','Devalued','PrePost','OdorType','IfDevalue','Pleasant'});
    writetable(AllData, '../ProcessedData/Pleasant_NODEAP_17.xlsx')
end

if if_process_choice_data
    AllChoiceData = AllChoiceData(:, {'SubID','StimLoc','StimOrder','Sess','Cond','Devalued','PrePost','Trial','Set','ChoiceType','OdorLR','CuePair','CueLeft','CueRight','SelectedCue','ChosenOdor','rt'});
    writetable(AllChoiceData, '../ProcessedData/Choices_NODEAP_17.xlsx')
end

if if_process_conditioning_data
    AllConditioningData = AllConditioningData(:, {'SubID','StimLoc','StimOrder','Sess','Cond','Devalued','Run','Trial','Set','ChoiceType','OdorLR','CuePair','CueLeft','CueRight','ChoiceLR','SelectedCue','OdorChosen','rt'});
    writetable(AllConditioningData, '../ProcessedData/Conditioning_NODEAP_17.xlsx')
end



