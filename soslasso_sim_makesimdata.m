function simdata = soslasso_sim_makesimdata(X_truth,G,sigma,varargin)
%  SOSLASSO_SIM_MAKESIMDATA Create data based on the templates.
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

if nargin > 3
    if islogical(varargin{1})
        VERBOSE = varargin{1};
    else
        error('soslasso_sim_makesimdata:Verbose flag needs to be true or false');
    end
else
    VERBOSE = false;
end

%% Add gaussian noise
P = uint32(length(X_truth));
T = uint32(size(X_truth{1},1));
N = uint32(size(X_truth{1},2));
[noise,X] = deal(cell(P,1));

if iscell(sigma)
% Hidden feature: if you want to reuse the noise from a previous simulation,
% which used exactly the same number of subjects (i.e. tasks), then the cell
% array of noise matrices can be passed to the sigma argument, rather than a
% scalar value. This noise will be used, rather than generating new random
% i.i.d. gaussian noise.
    noise = sigma;
    if ~length(X)==length(noise)
        error('soslasso_sim_makesimdata:Cell arrays X_truth and noise do not match.');
    end
    for i=1:length(X)
        X{i} = X_truth{i} + noise{i};
    end
else
    for i=1:length(X)
        noise{i} = (randn(T,N)*sigma);
        X{i} = X_truth{i} + noise{i};
    end
end

%% Define the replicated space
[RepIndex, groups, group_arr] = define_rep_space(G);

%% Create Y vectors.
Y = cell(P,1);
y = [ones(idivide(T,2,'floor'),1);-ones(idivide(T,2,'ceil'),1)];
Y(:) = deal({y});

simdata.Y = Y;
simdata.G = G;
simdata.X_truth = X_truth;
simdata.noise = noise;
simdata.X = X;
simdata.RepIndex= RepIndex;
simdata.groups = groups;
simdata.group_arr = group_arr;
simdata.sigma = sigma;
simdata.Date = date;
