function [DPrime,Counts,DiscVox,Betahat] = soslasso_sim_recoversignal(X,Y,ActiveVoxels,lambda,varargin)
    if nargin > 2
        if islogical(varargin{1})
            VERBOSE = varargin{1};
        else
            error('soslasso_sim_recoversignal:Verbose flag needs to be true or false');
        end
    else
        VERBOSE = false;
    end
    
    %% SOSLASSO
    [Betahat,C] = soslasso_solve(simdata,soslasso_params);
    [dp,counts] = soslasso_evaluate(ActiveVoxels,Betahat,C,varargin);
    soslasso.Betahat = Betahat;
    soslasso.C = C;
    soslasso.dp = dp;
    soslasso.counts = counts;
    
    %% LASSO
    [Betahat,C] = lasso_solve(simdata,sosdata);
    [dp,counts] = lasso_evaluate(ActiveVoxels,Betahat,C,varargin);
    soslasso.Betahat = Betahat;
    soslasso.C = C;
    soslasso.dp = dp;
    soslasso.counts = counts;


    %% Univarite (FDR corrected)

    % Print mean results
    if VERBOSE
        structfun(@mean,DPrime,'Unif',0)
    end
end
%% SUB-FUNCTIONS
%% SOSLasso
function [Betahat,C] = soslasso_solve(simdata,soslasso_params)
	[Betahat.soslasso,C.soslasso] = overlap_2stage(1,...
    	simdata.Y,...
    	simdata.X,...
    	simdata.G,...
    	simdata.RepIndex,...
    	simdata.group_arr,...
    	simdata.groups,...
    	lambda);
end

function [dp,counts] = soslasso_evaluate(ActiveVoxels,Betahat,C,varargin)
	if nargin > 5
		error('soslasso_evaluate: Too many input arguments.')
	end
	Overall = varargin{2};
	nzbeta = abs(Betahat.soslasso)>0;
	if Overall
		[dp,counts] = dprime(any(ActiveVoxels)',any(nzbeta,2));
	else
		[dp,counts] = dprime(ActiveVoxels',nzbeta);
	end
end


function [Betahat,C] = lasso_solve(X,Y,lambda)
	[Betahat,C,~] = Logistic_Lasso(simdata.X,...
		simdata.Y,...
		lambda);
end

function [dp,counts] = lasso_evaluate(ActiveVoxels,Betahat,C)
	[dp,counts] = soslasso_evaluate(ActiveVoxels,Betahat,C,varargin);
end

function [] = univariate_solve()
	Logistic_Lasso(X,Y,0);
% [MEAN_A, MEAN_B] = deal(zeros(P,N));
% [p_individual,h_individual] = deal(zeros(N,P));
% for i=1:P
%     a = simdata.X{i}(1:idivide(T,2,'floor'),:);
%     b = simdata.X{i}(idivide(T,2,'floor')+1:end,:);
%     [~,p_individual(:,i),~,stats] = ttest2(a,b);
%     h_individual(:,i) = fdr_bh(p_individual(:,i));
%     MEAN_A(i,:) = mean(a);
%     MEAN_B(i,:) = mean(b);
% end
% [~,p] = ttest2(MEAN_A,MEAN_B);
% [hh,~,~] = fdr_bh(p);
% temp = logical(hh);
% if isrow(temp)
%     DiscVox.univariate_overall = temp';
% else
%     DiscVox.univariate_overall = temp;
% end
% DiscVox.univariate = abs(h_individual)>0;

%% Compute 
% Separate for each subject

[DPrime.univariate,Counts.univariate] = dprime(cell2mat(ActiveVoxels),DiscVox.univariate);

% Aggregating over solutions, how well is the set of voxels active in any
% subject recovered?

[DPrime.univariate_overall,Counts.univariate_overall] = dprime(ACTIVES,DiscVox.univariate_overall);
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
