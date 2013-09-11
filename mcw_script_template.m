%% MCW EXPERIMENT SCRIPT
help mcw_help>how_to_exp

%% Prepare workspace
clear; clc; close all

% This is where things like dprime() reside.
addpath(genpath('~/matlab')) 

% This is where mcw_metadata.mat and original data resides.
addpath(genpath('~/mcw-dummydata'))

% This is where all my specific code resides.
addpath(genpath('~/mcw-soslasso')) 

%% Set parameters
% temp = linspace(2^-4,2^-1,20);
params.lamset       = exp(-4); % matches Nikhil.
% params.lamset       = 1; % for demonstration
params.numcvs       = 4;
params.classify     = true;
params.GroupSize    = uint32(3);
params.overlap      = true;
params.AccommodateGroups = true;
% specify a function of GroupSize or make the function yield a constant.
params.offset       = @(GS) idivide(GS,2,'floor');
params.Task         = 'ani-art';
params.BLUR         = false;
params.whatmethod   = 3; % soslasso
params.MeanCenter   = true;
params.NormVariance = false;
params.Save         = false;

params.RecoveryMode = 1;
% Recovery Modes:
% 0: Run all optimization routines from scratch.
% 1: load in prior Betahats and Cs from recovery.mat; skip CV module.
% 2: load in prior Betahats and Cs from recovery.mat adn recovery_FINAL.mat; skip all optimization, just compute final errors.

params.CVMode = 1;
% CV Modes:
% 0: local, cvs in serial
% 1: local, cvs in parallel
% 2: distributed manually, cvs in parallel. !!! Will stop running after CV
% module. Then compile all CV results to one computer, and proceed using
% recovery mode 1.

%% Load Metadata
metadata = load('metadata_dummy.mat',...
    'XYZ_tlrc',...
    'subjects',...
    'ReduxObj_aseg',...
    'Blocks', ...
    'Ambiguous',...
    'TrueAnimals');

%% Rename ReduxObj for compatibility throughout.
metadata.ReduxObj = metadata.ReduxObj_aseg;
metadata = rmfield(metadata,'ReduxObj_aseg');

%% Omit subjects
metadata.subjects(metadata.subjects==1225) = [];
metadata.subjects(metadata.subjects==2347) = [];
metadata.XYZ_tlrc = rmfield(metadata.XYZ_tlrc,'s1225');
metadata.XYZ_tlrc = rmfield(metadata.XYZ_tlrc,'s2347');
metadata.ReduxObj = rmfield(metadata.ReduxObj,'s1225');
metadata.ReduxObj = rmfield(metadata.ReduxObj,'s2347');

%% Identify Rows to include
good_rows = false(1200,length(metadata.subjects));
for s=1:length(metadata.subjects)
    S = sprintf('s%d',metadata.subjects(s));
    good_rows(:,s) = metadata.ReduxObj.(S).words;
end
good_rows = all(good_rows,2);
metadata.Words = [true(240,1); false(960,1)] & ~metadata.Ambiguous & good_rows;


%% Build CV blocks
metadata.Blocks = metadata.Blocks(metadata.Words,:);
metadata.CVBlocks = false(sum(metadata.Words),5);
b = 0;
for i = 1:2:9
    b = 1 + b;
    metadata.CVBlocks(:,b) = any([metadata.Blocks(:,i),metadata.Blocks(:,i+1)],2);
end
metadata.FinalTestBlock = metadata.CVBlocks(:,5);
metadata.CVBlocks(:,5) = [];
metadata = rmfield(metadata,'Blocks');
clear i b;

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
% params.RecoveryMode = ?;
if params.CVMode > 0
    [final_test,final_train,cv_test,cv_train] = mcw_classify_par(Y,X,GroupInfo,metadata,params);
else
    [final_test,final_train,cv_test,cv_train] = mcw_classify(Y,X,GroupInfo,metadata,params);
end
