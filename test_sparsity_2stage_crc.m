% testbench to gauge performance of different methods when noise is varied

%% Prepare workspace
% clear; clc; close all; restoredefaultpath;

% This is where things like dprime() reside.
addpath(genpath('~/matlab')) 

% This is where mcw_metadata.mat and original data resides.
addpath(genpath('~/mcw-data'))

% This is where all my specific code resides.
addpath(genpath('~/mcw-soslasso')) 

% This is where all Nikhil's "opt cluster only" code resides.
% addpath(genpath('~/mcw-nikhil')) 

%% Set Parameters
M = 500;
k = 20;
B = 6;

%% Start the clock...
tic

% generate groups
overlap = 1;
isgaussian = 1;
G = cell(0);
if overlap == 0
    for ii  = 1:M
        g = [(ii-1)*B+1:ii*B];
        G = [G;{[g]}];
    end
else
    for ii = 1:M
        if ii == 1;
            g = [1:B];
            lastind = B;
        else
            g = [lastind+1-B+1 : lastind+1];
            lastind = g(end);
        end
        G = [G;{[g]}];
    end
end

if isgaussian
    G = cell(0);
    if overlap == 0
        for ii  = 1:M
            g = [(ii-1)*B+1:ii*B];
            G = [G;{[g]}];
        end
    else
        for ii = 1:M
            if ii == 1;
                g = [1:B];
                lastind = B;
            else
                g = [lastind+1-floor(B/2)+1 : lastind+1+ceil(B/2)];
                lastind = g(end);
            end
            G = [G;{[g]}];
        end
    end
end
n = G{end}; n = n(end);

% test parameters
alpharange = .5;
noise = 0.1;
ntests = 1;
% lamset = linspace(1e-1,2,20);
lamset = .5;


% data parameters
numsubjects = 20;
numvoxels = n;
numtrials = 300;
% alpha = 0.8;
del_lev = 0.0;
nind = 0;
for alpha = alpharange
    nind = 1 + nind;
    
    for test = 1:ntests
        
        % generate data and add noise to measurements
        if ~isgaussian
            [Beta,X] = Generate_Data(numsubjects, numvoxels, numtrials, G, k, alpha, del_lev);
            Y = [ones(numtrials,numsubjects); -ones(numtrials,numsubjects)];
            Y = Y + noise*randn(size(Y));
            temp  = Y;
            Y = cell(0);
            for T = 1:numsubjects
                Y = [Y ; {[temp(:,T)]}];
            end
        else
            [Beta,X,Y] = Generate_Data_Gaussians(numsubjects, numvoxels, numtrials, G, k, alpha, 0);
            temp = Y;
            Y = cell(0);
            for T = 1:numsubjects
                y = temp{T};
                Y = [Y ; {[ y + noise*randn(size(y)) ]} ];
            end
        end
%         [Xo, groups, group_arr] = makeA_multitask_efficient(X,G);
        [RepIndex, groups, group_arr] = define_rep_space(G);

        
        % perform the tests
        i1 = 0;
        for lam1 = lamset
            i1 = 1+i1;
            
            %2 stage
            [Betahat,~] = overlap_2stage(1,Y,X,G,RepIndex,group_arr,groups, lam1);
            errlam(i1) = norm(Beta-Betahat,'fro')^2/numel(Betahat);
            nnz(Betahat)
%             
%             %1 stage (group lasso)
%             [Betahat,~] = overlap_1stage(0,Y,Xo,X,G,group_arr,lam1);
%             errlamg(i1) = norm(Beta-Betahat,'fro')^2/numel(Betahat);
%             
%             %lasso
%             [Betahat,~,] = Least_Lasso(X, Y, lam1);
%             errlaml(i1) = norm(Beta-Betahat,'fro')^2/numel(Betahat);
%             
        end
        
        
        errNESTED(test) = min(errlam);
%         errLASSO(test)  = min(errlaml);
%         errGLASSO(test) = min(errlamg);
        
        fprintf('.')
    end
    
    NESTED(nind) = mean(errNESTED);
%     LASSO(nind)  = mean(errLASSO);
%     GLASSO(nind) = mean(errGLASSO);
    sNESTED(nind)= std(errNESTED);
%     sLASSO(nind) = std(errLASSO);
%     sGLASSO(nind)= std(errGLASSO);
    
%     save COMPARE_METHODS_ALPHA_GAUSSIAN_L21 sLASSO sNESTED sGLASSO LASSO NESTED GLASSO
    
    fprintf('\n alpha  = %f \n',alpha);
end


