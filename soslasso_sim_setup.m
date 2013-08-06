function [simparams,X_truth,G] = soslasso_sim_setup(varargin)
%  SOSLASSO_SIM_SETUP Setup a simulation experiment and create a data
%  template.
%
%  USAGE:
%    simparams = SOSLASSO_SIM_SETUP() Returns default simparams structure, 
%    which can be modified and passed back into the function.
%
%    [simparams, X_truth, G] = SOSLASSO_SIM_SETUP(simparams) Using the information 
%    in simparams, generate the data template and group assignments, and log 
%    additional information in simparams.
%    
%    [...] = SOSLASSO_SIM_SETUP(simparams,'verbose') Same as above, but 
%    print plots and other useful information to the screen.
%
%    Key to Data Types: (N.B. Groups can be ``Modular'' by logical flag)
%    1 Same Sparse Groups
%    2 Shifted Sparse Groups
%    3 Different Sparse Groups
%    4 Identical No Groups
%    5 Different No Groups
%
%    If Modular==true, then active members within groups will all be 
%    responsive to the same category.
%
%    See also:
%    SOSLASSO_SIM_SETUP>KEY_VARIABLES

if nargin == 0
    simparams.nsubjects = uint32(16);	% number of subjects
    simparams.nvoxels = uint32(1024);   % number of voxels
    simparams.groupsize = uint32(64);	% group size
    simparams.groupshift = uint32(32);   % group shift (i.e. distance from G1(1)
                                        %     to G2(1) in voxel space.) 
    simparams.ntrials = uint32(64);     % number of trials
    simparams.nactgroups = uint32(4);   % number of active groups
    simparams.nactvoxels = uint32(8);   % number of active voxels per trial
    simparams.DataTypeInd = uint32(1);  % See help.
    simparams.Modular = logical(0);     % See help.
    return
end
% N.B. the group size and group shift define the ground truth---how the
% underlying signal is actually going to be created. The group size and shift
% considered by SOSlasso can be defined separately. These are separate because
% in practice the true groups are probably unknown, and we would like to know
% how robust SOSlasso is to searching over groups that don't exactly correspond
% to the truth.
if nargin > 1
    VERBOSE = true;
else
	VERBOSE = false;
end

simparams = varargin{1};
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
a = 1:simparams.groupshift:(simparams.nvoxels-simparams.groupsize+1);    % group start ind
b = (simparams.groupsize):simparams.groupshift:simparams.nvoxels;        % group end ind
ngroups = uint32(((simparams.nvoxels-simparams.groupsize)/simparams.groupshift)+1);
G = cell(ngroups,1);
for i=1:ngroups
    G{i} = a(i):b(i);
end
clear a b;

%% Define Data
% SOSLASSO_SIM_SETUP>KEY_VARIABLES
X_truth = cell(simparams.nsubjects,1);
switch simparams.DataType 
    case 'Same Sparse Groups'
        active_groups = uint32(randperm(ngroups,simparams.nactgroups));
        
        if simparams.Modular
            ani_groups = active_groups(1:idivide(simparams.nactgroups,2,'floor'));
            art_groups = active_groups(idivide(simparams.nactgroups,2,'ceil'):end);
            g_ani = unique(cell2mat(G(ani_groups))); % N.B. it is possible for ani and art
            g_art = unique(cell2mat(G(art_groups))); % voxels to overlap, whenever groups  
                                                     % overlap.
        else
            g = unique(cell2mat(G(active_groups)));
            V = uint32(length(g));
            g_randperm = g(randperm(V));
            g_ani = g_randperm(1:idivide(V,2,'floor')); % N.B. ani and art
            g_art = g_randperm(idivide(V,2,'ceil'):V);  % voxels will never overlap 
        end
        
        V_ani = uint32(length(g_ani));
        V_art = uint32(length(g_art));
        for i = 1:simparams.nsubjects
            X_truth{i} = zeros(simparams.ntrials,simparams.nvoxels);

            % Animal Loop (first half of trials)
            for j=1:idivide(simparams.ntrials,2,'floor');
                temp = uint32(randperm(V_ani,simparams.nactvoxels));
                active_voxels = g_ani(temp);
                X_truth{i}(j,active_voxels) = 1;
            end
            % Artifact Loop (second half of trials)
            for j=idivide(simparams.ntrials,2,'floor'):simparams.ntrials;
                temp = uint32(randperm(V_art,simparams.nactvoxels));
                active_voxels = g_art(temp);
                X_truth{i}(j,active_voxels) = 1;
            end
        end
        if VERBOSE
            disp(active_groups);
        end
    %%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  %%% 
    %%%% METHODS 2-5 NEED TO IMPLEMENT 2 CATEGORIES AND MODULARITY %%%%
    %%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  %%% 
    case 'Shifted Sparse Groups' 
        active_groups = uint32(randperm(ngroups,simparams.nactgroups));
        shift_log = zeros(1,simparams.nsubjects);
        for i=1:simparams.nsubjects
            ss = 2; % step size (minimum misalignment will be ss*simparams.groupshift)
            shift = randi(idivide(simparams.groupsize,simparams.groupshift*ss,'floor'))*ss*sign(randn(1));
            shift_log(i) = shift;
            active_groups = mod(active_groups + shift,ngroups) + 1;
            g = cell2mat(G(active_groups));
            g_set = unique(g);
            V = uint32(length(g_set));
            X_truth{i} = zeros(simparams.ntrials,simparams.nvoxels);
            for j=1:idivide(simparams.ntrials,2,'floor')
                temp = uint32(randperm(V,simparams.nactvoxels));
                active_voxels = g_set(temp);
                X_truth{i}(j,active_voxels) = 1;
            end
            if VERBOSE
                disp(active_groups)
            end
        end
 
    case 'Different Sparse Groups'
        for i=1:simparams.nsubjects
            X_truth{i} = zeros(simparams.ntrials,simparams.nvoxels);
            active_groups = uint32(randperm(ngroups,simparams.nactgroups));
            g = cell2mat(G(active_groups));
            g_set = unique(g);
            V = uint32(length(g_set));
            for j=1:idivide(simparams.ntrials,2,'floor')
                temp = uint32(randperm(V,simparams.nactvoxels));
                active_voxels = g_set(temp);
                X_truth{i}(j,active_voxels) = 1;
            end
            if VERBOSE
                disp(active_groups)
            end
        end
        
    case 'Identical No Groups'
        % Groups are selected only to derive V, number of voxels.
        active_groups = uint32(randperm(ngroups,simparams.nactgroups));
        g = cell2mat(G(active_groups));
        g_set = unique(g);
        V = uint32(length(g_set));
        % A random set of voxels is selected; patterns for each trial will
        % derive from v_set.
        v_set = uint32(randperm(simparams.nvoxels,V));
        for i=1:simparams.nsubjects
            X_truth{i} = zeros(simparams.ntrials,simparams.nvoxels);
            for j=1:idivide(simparams.ntrials,2,'floor')
                temp = uint32(randperm(V,simparams.nactvoxels));
                active_voxels = v_set(temp);
                X_truth{i}(j,active_voxels) = 1;
            end
        end

    
    case 'Different No Groups'
        % Groups are selected only to derive V, number of voxels.
        active_groups = uint32(randperm(ngroups,simparams.nactgroups));
        g = cell2mat(G(active_groups));
        g_set = unique(g);
        V = uint32(length(g_set));
        for i=1:simparams.nsubjects
            % A random set of voxels is selected; patterns for each trial will
            % derive from v_set.
            v_set = uint32(randperm(simparams.nvoxels,V));
            X_truth{i} = zeros(simparams.ntrials,simparams.nvoxels);
            for j=1:idivide(simparams.ntrials,2,'floor')
                temp = uint32(randperm(V,simparams.nactvoxels));
                active_voxels = v_set(temp);
                X_truth{i}(j,active_voxels) = 1;
            end
        end
end
if VERBOSE
    for i=1:6
        subplot(2,3,i);
        imagesc(X_truth{i})
    end
end

function key_variables()
% Key to variables:
% simparams.nsubjects: number of subjects           
% simparams.nvoxels: number of voxels
% simparams.groupsize: group size (uniform)       
% simparams.groupshift: space between start of each group (offset for overlap).
% simparams.ntrials: Number of trials  
% simparams.nactgroups: number of ``active'' groups.
% simparams.nactvoxels: The number of voxels that will be active on each trial.         
% V: Number of (unique) voxels across all selected groups.
% X_truth: The cell array holding the pattern for each subject.

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
