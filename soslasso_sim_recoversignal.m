function varargout = soslasso_sim_recoversignal(X,Y,ActiveVoxels,method,varargin)

    switch method
        case 'univariate'
            %% Univarite (FDR corrected)
            [h,p] = univariate_solve(X,Y,'Overall',true);
            [dp,counts] = univariate_evaluate(ActiveVoxels,h);
            varargout{1} = h;
            varargout{2} = NaN;
            varargout{3} = dp;
            varargout{4} = counts;

        case 'lasso'
            lambda = varargin{1};
            [Betahat,C] = lasso_solve(X,Y,lambda);
            [dp,counts] = lasso_evaluate(ActiveVoxels,Betahat,'Overall',true);
            varargout{1} = Betahat;
            varargout{2} = C;
            varargout{3} = dp;
            varargout{4} = counts;

        case 'soslasso'
            lambda = varargin{1};
            sosdata = varargin{2};
            [Betahat,C] = soslasso_solve(X,Y,lambda,sosdata);
            [dp,counts] = soslasso_evaluate(ActiveVoxels,Betahat,'Overall',true);
            varargout{1} = Betahat;
            varargout{2} = C;
            varargout{3} = dp;
            varargout{4} = counts;
    end

end



%% SUB-FUNCTIONS
%% SOSLasso
function [Betahat,C] = soslasso_solve(X,Y,lambda,sosdata)
	[Betahat,C,niter] = overlap_2stage(1,Y,X, ...
		sosdata.G,sosdata.RepIndex,sosdata.group_arr,sosdata.groups,lambda);
    fprintf('\n')
    fprintf('niter: %d',niter)
    fprintf('\n')
end

function [dp,counts] = soslasso_evaluate(ActiveVoxels,Betahat,varargin)
	if nargin > 5
		error('soslasso_evaluate: Too many input arguments.')
	end
	Overall = varargin{2};
	nzbeta = Betahat>0;
	if Overall
		[dp_all,counts] = dprime(ActiveVoxels,any(nzbeta,2));
        dp = mean(dp_all);
	else
		[dp,counts] = dprime(ActiveVoxels,nzbeta);
	end
end

%% Lasso
function [Betahat,C] = lasso_solve(X,Y,lambda)
	[Betahat,C,~] = Logistic_Lasso(X,Y,lambda);
end

function [dp,counts] = lasso_evaluate(ActiveVoxels,Betahat,varargin)
	[dp,counts] = soslasso_evaluate(ActiveVoxels,Betahat,varargin{:});
end

%% Univariate
function [h,p] = univariate_solve(X,Y,varargin)
	P = length(X);
	T = size(X{1},1);
	N = size(X{1},2);
	ani = Y{1}>0;
    Overall = varargin{2};
	if Overall
		[MEAN_X] = zeros(T,N);
		for i=1:P
			MEAN_X = X{i} + MEAN_X;
		end
		MEAN_X = MEAN_X ./ P;
		[~,p,~,stats] = ttest2(MEAN_X(ani,:),MEAN_X(~ani,:));
        t = stats.tstat;
		h = all([fdr_bh(p);t>0])';
	else
		[p,h,t] = deal(zeros(N,P));
		for i=1:P
			[~,p(:,i),~,stats] = ttest2(X{i}(ani,:),X{i}(~ani,:));
            t(:,i) = stats.tstat;
			h(:,i) = all([fdr_bh(p);t>0])';
		end
	end
end

function [dp,counts] = univariate_evaluate(ActiveVoxels,h,varargin)
	[dp_all,counts] = dprime(ActiveVoxels,h);
    dp = mean(dp_all);
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
