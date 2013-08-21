function [soslasso,lasso,univariate] = soslasso_sim_recoversignal(X,Y,ActiveVoxels,lambda,sosdata,varargin)
    if nargin > 5
        if islogical(varargin{1})
            VERBOSE = varargin{1};
        else
            error('soslasso_sim_recoversignal:Verbose flag needs to be true or false');
        end
    else
        VERBOSE = false;
    end
    
    %% SOSLASSO
    [Betahat,C] = soslasso_solve(X,Y,lambda,sosdata);
    [dp,counts] = soslasso_evaluate(ActiveVoxels,Betahat,'Overall',true);
    soslasso.Betahat = Betahat;
    soslasso.C = C;
    soslasso.dp = dp;
    soslasso.counts = counts;
    
    %% LASSO
    [Betahat,C] = lasso_solve(X,Y,lambda);
    [dp,counts] = lasso_evaluate(ActiveVoxels,Betahat,'Overall',true);
    lasso.Betahat = Betahat;
    lasso.C = C;
    lasso.dp = dp;
    lasso.counts = counts;

    %% Univarite (FDR corrected)
	[h,p] = univariate_solve(X,Y,'Overall',true);
    [dp,counts] = univariate_evaluate(ActiveVoxels,h);
    univariate.h = h;
    univariate.C = C;
    univariate.dp = dp;
    univariate.counts = counts;
end



%% SUB-FUNCTIONS
%% SOSLasso
function [Betahat,C] = soslasso_solve(X,Y,lambda,sosdata)
	[Betahat,C] = overlap_2stage(1,Y,X, ...
		sosdata.G,sosdata.RepIndex,sosdata.group_arr,sosdata.groups,lambda);
end

function [dp,counts] = soslasso_evaluate(ActiveVoxels,Betahat,varargin)
	if nargin > 5
		error('soslasso_evaluate: Too many input arguments.')
	end
	Overall = varargin{2};
	nzbeta = abs(Betahat)>0;
	if Overall
		[dp,counts] = dprime(any(ActiveVoxels)',any(nzbeta,2));
	else
		[dp,counts] = dprime(ActiveVoxels',nzbeta);
	end
end

%% Lasso
function [Betahat,C] = lasso_solve(X,Y,lambda)
	[Betahat,C,~] = Logistic_Lasso(X,Y,lambda);
end

function [dp,counts] = lasso_evaluate(ActiveVoxels,Betahat,varargin)
	[dp,counts] = soslasso_evaluate(ActiveVoxels,Betahat,C,varargin);
end

%% Univariate
function [h,p] = univariate_solve(X,Y,varargin)
	P = length(X);
	T = size(X{1},1);
	N = size(X{1},2);
	ani = Y{1}>0;
	if Overall
		[MEAN_X] = zeros(P,N);
		for i=1:P
			MEAN_X = X{i} + MEAN_X;
		end
		MEAN_X = MEAN_X ./ P;
		[~,p] = ttest2(MEAN_X(ani,:),MEAN_X(~ani,:));
		h = fdr_bh(p);
	else
		[p,h] = deal(zeros(N,P));
		for i=1:P
			[~,p(:,i)] = ttest2(X{i}(ani,:),X{i}(~ani,:));
			h(:,i) = fdr_bh(p);
		end
	end
end

function [dp,counts] = univariate_evaluate(ActiveVoxels,h,varargin)
	[dp,counts] = dprime(ActiveVoxels',h);
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
