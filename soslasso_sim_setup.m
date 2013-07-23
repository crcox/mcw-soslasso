function [simparams,X,G] = soslasso_sim_setup(varargin)
%  SOSLASSO_SIM_SETUP Setup a simulation experiment and create a data
%  template.
%
%  USAGE:
%    simparams = SOSLASSO_SIM_SETUP() Returns default simparams structure, 
%    which can be modified and passed back into the function.
%
%    [simparams, X, G] = SOSLASSO_SIM_SETUP(simparams) Using the information 
%    in simparams, generate the data template and group assignments, and log 
%    additional information in simparams.
%    
%    [...] = SOSLASSO_SIM_SETUP(simparams,'verbose') Same as above, but 
%    print plots and other useful information to the screen.
%
%    Key to Data Types:
%    1 Same Sparse Groups
%    2 Shifted Sparse Groups
%    3 Different Sparse Groups
%    4 Identical No Groups
%    5 Different No Groups

if nargin == 0
    simparams.nsubjects = uint16(16);	% number of subjects
    simparams.nvoxels = uint16(1024);   % number of voxels
    simparams.groupsize = uint16(64);	% group size
    simparams.groupspace = uint16(4);   % group space
    simparams.ntrials = uint16(64);     % number of trials
    simparams.nactgroups = uint16(4);   % number of active groups
    simparams.nactvoxels = uint16(8);   % number of active voxels per trial
    simparams.DataTypeInd = uint16(1);  % See help.
    return
else
    simparams = varargin{1};
end
if nargin > 1
    VERBOSE = true;
end

if simparams.DataTypeInd == 1 
	simparams.DataType = 'Same Sparse Groups';
elseif simparams.DataTypeInd == 2 
	simparams.DataType = 'Shifted Sparse Groups';
elseif simparams.DataTypeInd == 3 
	simparams.DataType = 'Different Sparse Groups';
elseif simparams.DataTypeInd == 4 
	simparams.DataType = 'Identical No Groups';
elseif simparams.DataTypeInd == 5 
	simparams.DataType = 'Different No Groups';
end

simparams.StartDate = date;

%% Define Groups
a = 1:simparams.groupspace:(simparams.nvoxels-simparams.groupsize+1);    % group start ind
b = (simparams.groupsize):simparams.groupspace:simparams.nvoxels;        % group end ind
ngroups = uint16(((simparams.nvoxels-simparams.groupsize)/simparams.groupspace)+1);
G = cell(1,ngroups);
for i=1:ngroups
    G{i} = a(i):b(i);
end
clear a b;

%% Define Data
% Key to variables:
% simparams.nsubjects: number of subjects           
% simparams.nvoxels: number of voxels
% simparams.groupsize: group size (uniform)       
% simparams.groupspace: space between start of each group (offset for overlap).
% simparams.ntrials: Number of trials  
% simparams.nactgroups: number of ``active'' groups.
% simparams.nactvoxels: The number of voxels that will be active on each trial.         
% V: Number of (unique) voxels across all selected groups.
% X: The cell array holding the pattern for each subject.
% X_noise: X{i} + randn(size(X{i}))*sigma
% R: X_noise, replicated (see makeA_multitask_efficient.m for routine).

% Explanation of data types:
% Same Sparse Groups: The group selection is the same for each subject, but
% the patterns vary within and across subjects. Sparsity within and among
% groups, but the selection is applied to all subjects.
%
% Shifted Sparse Groups: One set of active groups is chosen before entering
% the subject loop, but the active groups for each subject are randomly
% offset (shifted positively or negatively by some amount).  This means
% that the selected groups will be anatomically similar, but not identical
% across subjects. The voxels active on each trial are randomly sampled
% from the active groups.

X = cell(simparams.nsubjects,1);
switch simparams.DataType 
    case 'Same Sparse Groups'
        active_groups = uint16(randperm(ngroups,simparams.nactgroups));
        g = cell2mat(G(active_groups));
        g_set = unique(g);
        V = uint16(length(g_set));
        for i=1:simparams.nsubjects
            X{i} = zeros(simparams.ntrials,simparams.nvoxels);
            for j=1:idivide(simparams.ntrials,2,'floor')
                temp = uint16(randperm(V,simparams.nactvoxels));
                active_voxels = g_set(temp);
                X{i}(j,active_voxels) = 1;
            end
        end
        if VERBOSE
            disp(active_groups);
        end
        
    case 'Shifted Sparse Groups'
        active_groups = uint16(randperm(ngroups,simparams.nactgroups));
        shift_log = zeros(1,simparams.nsubjects);
        for i=1:simparams.nsubjects
            ss = 2; % step size (minimum misalignment will be ss*simparams.groupspace)
            shift = randi(idivide(simparams.groupsize,simparams.groupspace*ss,'floor'))*ss*sign(randn(1));
            shift_log(i) = shift;
            active_groups = mod(active_groups + shift,ngroups) + 1;
            g = cell2mat(G(active_groups));
            g_set = unique(g);
            V = uint16(length(g_set));
            X{i} = zeros(simparams.ntrials,simparams.nvoxels);
            for j=1:idivide(simparams.ntrials,2,'floor')
                temp = uint16(randperm(V,simparams.nactvoxels));
                active_voxels = g_set(temp);
                X{i}(j,active_voxels) = 1;
            end
            if VERBOSE
                disp(active_groups)
            end
        end
 
    case 'Different Sparse Groups'
        for i=1:simparams.nsubjects
            X{i} = zeros(simparams.ntrials,simparams.nvoxels);
            active_groups = uint16(randperm(ngroups,simparams.nactgroups));
            g = cell2mat(G(active_groups));
            g_set = unique(g);
            V = uint16(length(g_set));
            for j=1:idivide(simparams.ntrials,2,'floor')
                temp = uint16(randperm(V,simparams.nactvoxels));
                active_voxels = g_set(temp);
                X{i}(j,active_voxels) = 1;
            end
            if VERBOSE
                disp(active_groups)
            end
        end
        
    case 'Identical No Groups'
        % Groups are selected only to derive V, number of voxels.
        active_groups = uint16(randperm(ngroups,simparams.nactgroups));
        g = cell2mat(G(active_groups));
        g_set = unique(g);
        V = uint16(length(g_set));
        % A random set of voxels is selected; patterns for each trial will
        % derive from v_set.
        v_set = uint16(randperm(simparams.nvoxels,V));
        for i=1:simparams.nsubjects
            X{i} = zeros(simparams.ntrials,simparams.nvoxels);
            for j=1:idivide(simparams.ntrials,2,'floor')
                temp = uint16(randperm(V,simparams.nactvoxels));
                active_voxels = v_set(temp);
                X{i}(j,active_voxels) = 1;
            end
        end

    
    case 'Different No Groups'
        % Groups are selected only to derive V, number of voxels.
        active_groups = uint16(randperm(ngroups,simparams.nactgroups));
        g = cell2mat(G(active_groups));
        g_set = unique(g);
        V = uint16(length(g_set));
        for i=1:simparams.nsubjects
            % A random set of voxels is selected; patterns for each trial will
            % derive from v_set.
            v_set = uint16(randperm(simparams.nvoxels,V));
            X{i} = zeros(simparams.ntrials,simparams.nvoxels);
            for j=1:idivide(simparams.ntrials,2,'floor')
                temp = uint16(randperm(V,simparams.nactvoxels));
                active_voxels = v_set(temp);
                X{i}(j,active_voxels) = 1;
            end
        end
end
if VERBOSE
    for i=1:6
        subplot(2,3,i);
        imagesc(X{i})
    end
end
