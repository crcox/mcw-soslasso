function simdata = soslasso_sim_makesimdata(X,G,sigma,varargin)
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
end

if nargin > 3
    if islogical(varargin{1})
        VERBOSE = varargin{1};
    else
        error('soslasso_sim:Verbose flag needs to be true or false');
    end
else
    VERBOSE = false;
end

%% Extract parameters
P = length(X);
[T,N] = size(X{1});
K = length(G);
M = length(G{1});
L = G{2}(1) - G{1}(1);

%% Add gaussian noise
X_noise = cell(P,1);
for i=1:length(X)
    X_noise{i} = X{i} + (randn(T,N)*sigma);
end

%% Replicate (for SOS Lasso)
allvox = 1:N;
replication_index = cell2mat(G);
groups = repmat(uint16(0), M*K, 1);
group_arr = repmat(N+1,K,M);

a = 0; %#ok
b = 0;
ii = 0;
for j=1:K % number of groups
	group_arr(j,1:sum(temp)) = (ii+1):(ii+sum(temp));
	ii = ii + sum(temp);
	a = b + 1;
	b = (a + sum(temp)) - 1;
	groups(a:b) = j;
end

Y = cell(P,1);
y = [ones(idivide(T,2,'floor'),1),-ones(idivide(T,2,'ceil'),1)];
Y(:) = deal({y});

simdata.Y = Y;
simdata.G = G;
simdata.X = X;
simdata.X_noise;
simdata.replication_index = replication_index;
simdata.groups = groups;
simdata.group_arr;
simdata.sigma = sigma;
simdata.Date = date;