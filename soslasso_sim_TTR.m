function [DPrime, Counts, DiscVox, ActiveVoxels, Betahat, X_noise] = soslasso_sim(sigma, lam, DataTypeInd,varargin)
%  SOSLASSO_SIM Simulate different data with different structures, 
%    bury it in Gaussian noise, and recover using soslasso, lasso,
%    and univariate method.
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
if nargin > 3
    if islogical(varargin{1})
        VERBOSE = varargin{1};
    else
        error('soslasso_sim:Verbose flag needs to be true or false');
    end
else
    VERBOSE = false;
end
P = uint16(16);     % number of people
N = uint16(1024);   % number of voxels
M = uint16(64);     % group size
L = uint16(4);      % group space
T = uint16(64);     % trials

A = 1:L:(N-M+1);    % group start ind
B = (M):L:N;        % group end ind
K = uint16(length(A));
G = cell(1,K);
for i=1:K
    G{i} = A(i):B(i);
end
clear A B;

J = uint16(64);     % K/J groups will be active
I = idivide(K,J,'floor');
H = uint16(idivide(J,4,'floor'));    % V/H voxels will be active.

y = [ones(idivide(T,2,'floor'),1);-ones(idivide(T,2,'ceil'),1)];


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

%% Identify active voxels by subject
ActiveVoxels = cellfun(@any,X,'Unif',0);
ACTIVES = ascolumn(any(cell2mat(ActiveVoxels)));
ActiveVoxels = cellfun(@transpose,ActiveVoxels,'Unif',0)';

%% Add gaussian noise
% sigma = 0.05;
X_noise = cell(P,1);
for i=1:P
    X_noise{i} = (X{i}-0.5) + (randn(size(X{i}))*sigma);
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
% lam = 0.001;
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

%% SOSLasso
[Betahat.soslasso,C.soslasso] = overlap_2stage(1,Y,R,X_noise,G,group_arr,groups,lam);
DiscVox.soslasso_overall = ascolumn(any(abs(Betahat.soslasso)>0,2));
DiscVox.soslasso = abs(Betahat.soslasso)>0;

%% Lasso
[Betahat.lasso,C.soslasso,~] = Logistic_Lasso(X_noise, Y, lam);
DiscVox.lasso_overall = ascolumn(any(abs(Betahat.lasso)>0,2));
DiscVox.lasso = abs(Betahat.lasso)>0;

%% Univarite (FDR corrected)
[MEAN_A, MEAN_B] = deal(zeros(P,N));
[p_individual,h_individual] = deal(zeros(N,P));
for i=1:P
    a = X_noise{i}(1:idivide(T,2,'floor'),:);
    b = X_noise{i}(idivide(T,2,'floor')+1:end,:);
    [~,p_individual(:,i)] = ttest2(a,b);
    h_individual(:,i) = fdr_bh(p_individual(:,i));
    MEAN_A(i,:) = mean(a);
    MEAN_B(i,:) = mean(b);
end
[~,p] = ttest2(MEAN_A,MEAN_B);
[hh,~,~] = fdr_bh(p);
DiscVox.univariate_overall = ascolumn(logical(hh));
DiscVox.univariate = abs(h_individual)>0;

%% Compute 
% Separate for each subject
[DPrime.soslasso,Counts.soslasso] = dprime(cell2mat(ActiveVoxels),DiscVox.soslasso);
[DPrime.lasso,Counts.lasso] = dprime(cell2mat(ActiveVoxels),DiscVox.lasso);
[DPrime.univariate,Counts.univariate] = dprime(cell2mat(ActiveVoxels),DiscVox.univariate);

% Aggregating over solutions, how well is the set of voxels active in any
% subject recovered?
[DPrime.soslasso_overall,Counts.soslasso_overall] = dprime(ACTIVES,DiscVox.soslasso_overall);
[DPrime.lasso_overall,Counts.lasso_overall] = dprime(ACTIVES,DiscVox.lasso_overall);
[DPrime.univariate_overall,Counts.univariate_overall] = dprime(ACTIVES,DiscVox.univariate_overall);

% Print mean results.
if VERBOSE
    structfun(@mean,DPrime,'Unif',0)
end

% In citing MALSAR in your papers, you can use the following:
% [Zhou 2012] J. Zhou, J. Chen and J. Ye. MALSAR: Multi-tAsk Learning via 
%     StructurAl Regularization. Arizona State University, 2012. 
%     http://www.public.asu.edu/~jye02/Software/MALSAR.
% 
% If you use LaTeX, you can download the BibTex entry: 
%     http://www.public.asu.edu/~jye02/Software/MALSAR/citation.bib.
% 
% The MALSAR software project has been supported by research grants from 
% the National Science Foundation (NSF). 
