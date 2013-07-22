function [simparams,X] = soslasso_sim_setup(DataTypeInd,varargin)
%  SOSLASSO_SIM Setup data with various structures, 
%
%  USAGE:
%    [DPrime, Counts] = SOSLASSO_SIM(sigma, lam, DataTypeInd,...);
%
%    Key to Data Types:
%    1 Same Sparse Groups
%    2 Shifted Sparse Groups
%    3 Different Sparse Groups
%    4 Identical No Groups
%    5 Different No Groups

%% Specify Simulation environment
if DataTypeInd == 1 
	DataType = 'Same Sparse Groups';
elseif DataTypeInd == 2 
	DataType = 'Shifted Sparse Groups';
elseif DataTypeInd == 3 
	DataType = 'Different Sparse Groups';
elseif DataTypeInd == 4 
	DataType = 'Identical No Groups';
elseif DataTypeInd == 5 
	DataType = 'Different No Groups';
end

if nargin > 1
    if islogical(varargin{1})
        VERBOSE = varargin{1};
    else
        error('soslasso_sim:Verbose flag needs to be true or false');
    end
else
    P = uint16(16);     % number of people
    N = uint16(1024);   % number of voxels
    M = uint16(64);     % group size
    L = uint16(4);      % group space
    T = uint16(64);     % trials
    a = 1:L:(N-M+1);    % group start ind
    b = (M):L:N;        % group end ind
    K = uint16(length(b));
    G = cell(1,K);
    for i=1:K
        G{i} = a(i):b(i);
    end
    clear a b;
    J = uint16(64);     % K/J groups will be active
    I = idivide(K,J,'floor');
    H = uint16(idivide(J,4,'floor'));    % V/H voxels will be active.
    y = [ones(idivide(T,2,'floor'),1);-ones(idivide(T,2,'ceil'),1)];
end

if nargin > 2
    if islogical(varargin{1})
        VERBOSE = varargin{1};
    else
        error('soslasso_sim:Verbose flag needs to be true or false');
    end
else
    VERBOSE = false;
end


%% Define Data
% Key to variables:
% P: number of people           
% N: number of voxels
% M: group size (uniform)       
% L: space between start of each group (offset for overlap).
% K: number of ``active'' groups.
% T: Number of trials           
% V: Number of (unique) voxels across all selected groups.
% W: The number of voxels that will be active on each trial.
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

X = cell(P,1);
switch DataType 
    case 'Same Sparse Groups'
        active_groups = uint16(randperm(K,I));
        g = cell2mat(G(active_groups));
        g_set = unique(g);
        V = uint16(length(g_set));
        W = idivide(V,H,'floor');
        for i=1:P
            X{i} = zeros(T,N);
            for j=1:idivide(T,2,'floor')
                temp = uint16(randperm(V,W));
                active_voxels = g_set(temp);
                X{i}(j,active_voxels) = 1;
            end
        end
        
    case 'Shifted Sparse Groups'
        active_groups = uint16(randperm(K,I));
        shift_log = zeros(1,P);
        for i=1:P
            ss = 2; % step size (minimum misalignment will be ss*L)
            shift = randi(idivide(M,L*ss,'floor'))*ss*sign(randn(1));
            shift_log(i) = shift;
            active_groups = circshift(active_groups, shift);
            g = cell2mat(G(active_groups));
            g_set = unique(g);
            V = uint16(length(g_set));
            W = idivide(V,H,'floor');
            X{i} = zeros(T,N);
            for j=1:idivide(T,2,'floor')
                temp = uint16(randperm(V,W));
                active_voxels = g_set(temp);
                X{i}(j,active_voxels) = 1;
            end
        end
 
    case 'Different Sparse Groups'
        for i=1:P
            X{i} = zeros(T,N);
            active_groups = uint16(randperm(K,I));
            g = cell2mat(G(active_groups));
            g_set = unique(g);
            V = uint16(length(g_set));
            W = idivide(V,H,'floor');
            for j=1:idivide(T,2,'floor')
                temp = uint16(randperm(V,W));
                active_voxels = g_set(temp);
                X{i}(j,active_voxels) = 1;
            end
        end
        
    case 'Identical No Groups'
        % Groups are selected only to derive V, number of voxels.
        active_groups = uint16(randperm(K,I));
        g = cell2mat(G(active_groups));
        g_set = unique(g);
        V = uint16(length(g_set));
        W = idivide(V,H,'floor');
        % A random set of voxels is selected; patterns for each trial will
        % derive from v_set.
        v_set = uint16(randperm(N,V));
        for i=1:P
            X{i} = zeros(T,N);
            for j=1:idivide(T,2,'floor')
                temp = uint16(randperm(V,W));
                active_voxels = v_set(temp);
                X{i}(j,active_voxels) = 1;
            end
        end

    
    case 'Different No Groups'
        % Groups are selected only to derive V, number of voxels.
        active_groups = uint16(randperm(K,I));
        g = cell2mat(G(active_groups));
        g_set = unique(g);
        V = uint16(length(g_set));
        W = idivide(V,H,'floor');
        for i=1:P
            % A random set of voxels is selected; patterns for each trial will
            % derive from v_set.
            v_set = uint16(randperm(N,V));
            X{i} = zeros(T,N);
            for j=1:idivide(T,2,'floor')
                temp = uint16(randperm(V,W));
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

simparams.P = P; % number of people           
simparams.N = N; % number of voxels
simparams.M = M; % group size (uniform)       
simparams.L = L; % space between start of each group (offset for overlap).
simparams.K = K; % number of ``active'' groups.
simparams.T = T; % Number of trials           
simparams.V = V; % Number of (unique) voxels across all selected groups.
simparams.W = W; % The number of voxels that will be active on each trial.
simparams.DataType = DataType;
simparams.DataTypeInd = DataTypeInd;

