function [DPrime, Counts, DiscVox, ActiveVoxels] = soslasso_sim_recoversignal(simdata,lam,varargin)

if nargin > 2
    if islogical(varargin{1})
        VERBOSE = varargin{1};
    else
        error('soslasso_sim:Verbose flag needs to be true or false');
    end
else
    VERBOSE = false;
end

%% Import and check simdata
Y = simdata.Y;
R = simdata.R;
X_noise = simdata.X_noise;
G = simdata.G;
group_arr = simdata.group_arr;
groups = simdata.groups;

%% Identify active voxels by subject
ActiveVoxels = cellfun(@any,X,'Unif',0);
ACTIVES = ascolumn(any(cell2mat(ActiveVoxels)));
ActiveVoxels = cellfun(@transpose,ActiveVoxels,'Unif',0)';

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
