%% Step 0: Create a new directory.
% It's probably best to create a directory for each simulation effort you
% make. Make one directory that is a test bed for experimenting with
% different parameters. Nothing in this directory "matters". It's safe to
% assume any results in this directory could be deleted, and no one would
% cry themselves to sleep. Once you have things set up the way you want
% them, create a new directory, copy over the sim script that has
% everything just the way you want it, and then execute it all over again.
% Ideally, this will keep things neat, and replicable.  

%% Prepare the workspace
% This will completely clear the workspace---either skip this line or save
% anything valuable.
clear;     % clear workspace
clc;       % clear the console
close all; % close all plots
restoredefaultpath; % ensure nothing funky is on the search path

addpath(genpath('~/matlab'));       % local non-specific code
addpath(genpath('~/mcw-soslasso')); % project specific code

%% Set Parameters
simparams = soslasso_sim_setup();
% By default:
% simparams = 
% 
%       nsubjects: 16
%         nvoxels: 1024
%       groupsize: 64
%      groupshift: 32
%         ntrials: 64
%      nactgroups: 4
%      nactvoxels: 8
%     DataTypeInd: 1
%         Modular: 0
%
% Modify these parameters to change the simulation. For example:
%     simparams.nsubjects = 12;
% would change the number of datasets created. Run soslasso_sim_help for
% documentation. All of the functions have associated "help". Typing 
% "help <function>" without the quotes will pull it up.
%
% soslasso_sim_help

%% Generate X_truth and G (for grouping)
[simparams, X_truth, Y] = soslasso_sim_setup(simparams);
save('simsetup.mat','simparams','X_truth','Y');

%% Setup data for SOSLasso
GroupSize  = 64;
GroupShift = 32;
a = 1:GroupShift:(simparams.nvoxels-GroupSize+1); % group start ind
b = GroupSize:GroupShift:simparams.nvoxels; % group end ind
G = cell(length(a),1);
for i=1:length(a)
    G{i} = a(i):b(i);
end
clear a b;

[RepIndex, groups, group_arr] = define_rep_space(G);
sosdata.G = G;
sosdata.RepIndex = RepIndex;
sosdata.groups = groups;
sosdata.group_arr = group_arr;
save('sosdata.mat','sosdata');

%% Bury signal in gaussian noise
sigma = 0.5;
X = cellfun(@(x) x + randn(size(x))*sigma, X_truth, 'Unif', 0);

%% Define lambda
lambda = 0.1;

%% Identify Active Voxels by Subjects
ActiveVoxels = cell2mat(cellfun(@any,X_truth, 'Unif', 0));

%% Recover Signal
[soslasso,lasso,univariate] = soslasso_sim_recoversignal(X,Y,ActiveVoxels,lambda,sosdata);
