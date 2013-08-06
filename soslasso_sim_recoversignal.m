function [DPrime,Counts,DiscVox,Betahat] = soslasso_sim_recoversignal(simdata,lambda,varargin)

if nargin > 2
    if islogical(varargin{1})
        VERBOSE = varargin{1};
    else
        error('soslasso_sim_recoversignal:Verbose flag needs to be true or false');
    end
else
    VERBOSE = false;
end

%% Identify active voxels by subject
ActiveVoxels = cellfun(@any,simdata.X_truth,'Unif',0);
temp = any(cell2mat(ActiveVoxels));
if isrow(temp)
    ACTIVES = temp';
else
    ACTIVES = temp;
end
ActiveVoxels = cellfun(@transpose,ActiveVoxels,'Unif',0)';

% implement cross validation
% think about how to store data
% separate structures for final solutions and CV solutions.

P = uint32(length(simdata.X));
T = uint32(size(simdata.X{1},1));
N = uint32(size(simdata.X{1},2));
K = uint32(length(simdata.G));
M = uint32(length(simdata.G{1}));
L = simdata.G{2}(1) - simdata.G{1}(1);

%% SOSLasso
[Betahat.soslasso,C.soslasso] = overlap_2stage(1,...
    simdata.Y,...
    simdata.X,...
    simdata.G,...
    simdata.RepIndex,...
    simdata.group_arr,...
    simdata.groups,...
    lambda);
temp = any(abs(Betahat.soslasso)>0,2);
if isrow(temp)
    DiscVox.soslasso_overall = temp';
else
    DiscVox.soslasso_overall = temp;
end
DiscVox.soslasso = abs(Betahat.soslasso)>0;
[DPrime.soslasso,Counts.soslasso] = dprime(cell2mat(ActiveVoxels),DiscVox.soslasso);
[DPrime.soslasso_overall,Counts.soslasso_overall] = dprime(ACTIVES,DiscVox.soslasso_overall);

%% Lasso
[Betahat.lasso,C.soslasso,~] = Logistic_Lasso(...
    simdata.X,...
    simdata.Y,...
    lambda);
temp = any(abs(Betahat.lasso)>0,2);
if isrow(temp)
    DiscVox.lasso_overall = temp';
else
    DiscVox.lasso_overall = temp;
end
DiscVox.lasso = abs(Betahat.lasso)>0;
[DPrime.lasso,Counts.lasso] = dprime(cell2mat(ActiveVoxels),DiscVox.lasso);
[DPrime.lasso_overall,Counts.lasso_overall] = dprime(ACTIVES,DiscVox.lasso_overall);

%% Univarite (FDR corrected)
[MEAN_A, MEAN_B] = deal(zeros(P,N));
[p_individual,h_individual] = deal(zeros(N,P));
for i=1:P
    a = simdata.X{i}(1:idivide(T,2,'floor'),:);
    b = simdata.X{i}(idivide(T,2,'floor')+1:end,:);
    [~,p_individual(:,i),~,stats] = ttest2(a,b);
    h_individual(:,i) = fdr_bh(p_individual(:,i));
    MEAN_A(i,:) = mean(a);
    MEAN_B(i,:) = mean(b);
end
[~,p] = ttest2(MEAN_A,MEAN_B);
[hh,~,~] = fdr_bh(p);
temp = logical(hh);
if isrow(temp)
    DiscVox.univariate_overall = temp';
else
    DiscVox.univariate_overall = temp;
end
DiscVox.univariate = abs(h_individual)>0;

%% Compute 
% Separate for each subject

[DPrime.univariate,Counts.univariate] = dprime(cell2mat(ActiveVoxels),DiscVox.univariate);

% Aggregating over solutions, how well is the set of voxels active in any
% subject recovered?

[DPrime.univariate_overall,Counts.univariate_overall] = dprime(ACTIVES,DiscVox.univariate_overall);

% Print mean results
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
