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
[simparams, X_truth, G] = soslasso_sim_setup(simparams);

%% Generate X, and data elements
simdata = soslasso_sim_makesimdata(X_truth,G,simparams);

%%



mkdir('~/'