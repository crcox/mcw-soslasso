%% MCW EXPERIMENT SCRIPT
help mcw_help>how_to_exp

%% Prepare workspace
clear; clc; close all

% This is where things like dprime() reside.
addpath(genpath('~/matlab')) 

% This is where mcw_metadata.mat and original data resides.
addpath(genpath('~/mcw-data'))

% This is where all my specific code resides.
addpath(genpath('~/mcw-soslasso')) 

%% Set parameters
params.lamset       = 10:20:190;
params.classify     = true;
params.numcvs       = uint32(4);
params.GroupSize    = uint32(21);
params.Task         = 'ani-art';
params.BLUR         = false;
params.whatmethod   = 3; % 1 = lasso, 2 = glasso, 3 = soslasso, 4 = oglasso
params.MeanCenter   = true;
params.NormVariance = false;
params.Save         = true;

%% Load Metadata
metadata = load('mcw_metadata.mat',...
    'XYZ_tlrc',...
    'subjects',...
    'ReduxObj',...
    'Blocks', ...
    'Ambiguous',...
    'TrueAnimals');

%% Omit subjects
metadata.subjects(metadata.subjects==1225) = [];
metadata.subjects(metadata.subjects==2347) = [];

%% Identify Subset
metadata.Words = [true(240,1); false(960,1)] & ~metadata.Ambiguous;

%% Build CV blocks
metadata.Blocks = metadata.Blocks(metadata.Words,:);
metadata.CVBlocks = false(sum(metadata.Words),5);
for i = 1:5
    metadata.CVBlocks(:,i) = any([metadata.Blocks(:,i),metadata.Blocks(:,i+5)],2);
end
clear i;

%% Specify the targets for each subject
Y = cell(length(metadata.subjects),1);
for i = 1:length(metadata.subjects);
    temp = metadata.TrueAnimals(metadata.Words);
    Y{i} = -ones(length(temp),1);
    Y{i}(temp) = 1;
end
clear temp i;

%% Set-Up Data
[X,GroupInfo] = mcw_prep_data(metadata,params);

%% Run Experiment
[final_test,final_train,cv_test,cv_train] = classify_CRC(Y,X,GroupInfo,metadata,params);
final_test
