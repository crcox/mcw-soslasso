function simdata = soslasso_sim_makesimdata(X,simparams,varargin)

if nargin > 2
    if islogical(varargin{1})
        VERBOSE = varargin{1};
    else
        error('soslasso_sim:Verbose flag needs to be true or false');
    end
else
    VERBOSE = false;
end

fields = fieldnames(simparams);
for i = 1:length(fields);
    sprintf('%s = simparams.%s',fields{[i i]});
end

%% Add gaussian noise
% sigma = 0.05;
X_noise = cell(P,1);
for i=1:P
    X_noise{i} = X{i} + (randn(size(X{i}))*sigma);
end
% for i=1:6
%    subplot(2,3,i);
%    imagesc(X_noise{i})
% end

%% Replicate (for SOS Lasso)
R = cell(P,1);
R(:) = deal({zeros(T,M*K)});
allvox = 1:N;
groups = repmat(uint16(0), M*K, 1);
group_arr = repmat(N+1,K,M);

% First loop
for i=1
    a = 0; %#ok
    b = 0;
    ii = 0;
    for j=1:K % number of groups
        temp = ismember(allvox,G{j});
        group_arr(j,1:sum(temp)) = (ii+1):(ii+sum(temp));
        ii = ii + sum(temp);
        a = b + 1;
        b = (a + sum(temp)) - 1;
        R{i}(:,a:b) = X_noise{i}(:,temp);
        groups(a:b) = j;
    end
end

% Remaining loops
for i=2:P
    a = 0; %#ok
    b = 0;
    for j=1:K % number of groups
        temp = ismember(allvox,G{j});
        a = b + 1;
        b = (a + sum(temp)) - 1;
        R{i}(:,a:b) = X_noise{i}(:,temp);
    end
end

Y = cell(P,1);
Y(:) = deal({y});

[Xo,gs,ga]= makeA_multitask_efficient(X_noise,G'); % This just allows me to confirm that Nikhil and I are on the same page.
for i=1:P
    if ~all(all(Xo{i}==R{i}))
        error('soslasso_sim:Something is wrong with cell array R.')
    end
end

if ~all(groups==gs)
    error('soslasso_sim:Something is wrong with vector groups.')
end

if ~all(all(group_arr==ga))
    error('soslasso_sim:Something is wrong with matrix group_arr.')
end

simdata.Y = Y;
simdata.X = X;
simdata.X_noise;
simdata.R = R;
simdata.G = G;
simdata.groups = groups;
simdata.group_arr;