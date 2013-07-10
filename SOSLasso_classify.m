function SOSLasso_classify(modelInfo)
%SOSLasso_classify   Run the Sparse Overlapping Sets Lasso algorithm.
%   SOSLasso is a multitask optimization routine that is a generalization of
%   group lasso.  Groups are defined based on the features of several tasks.
%   These groups may overlap. If a group is selected in one subject, some sparse
%   set of features within that group must be used in all other subjects.
%   A sparsity pattern is disovered over groups and within groups.
%
%   SOSLasso_classify(modelInfo) solves the optimization problem based on the
%   conditions and parameters specified in the structured array modelInfo.
%
%   See also MEDIAN, STD, MIN, MAX, VAR, COV, MODE.

%   Nikhil Rao, Christopher Cox, Robert Nowak, Timothy Rogers 
%   Date: 2013/07/08

disp(modelInfo.cv);
if ( modelInfo.cv < 1 || modelInfo.cv > 9 )
    error('cv should be an integer from 1-9 \n');
end

try
    DataDir = modelInfo.DataDir;
catch ME
    warning('No data directory specified. Assuming data resides in current directory...\n (To specify another directory, set modelInfo.DataDir="some/path/")\n');
    DataDir = pwd;
end 

tic;
cv            = modelInfo.cv;                   %    1x1  double
classify      = modelInfo.classify;             %    1x1  double
whatmethod    = modelInfo.ClassificationMethod; %    1x1  double
nums          = modelInfo.subjects;             %    1x14 double
Ambiguous     = modelInfo.Ambiguous;            % 1200x1  logical
Blocks        = modelInfo.Blocks;               % 1200x10 logical
TrueAnimals   = modelInfo.TrueAnimals;          % 1200x1  logical
TrueArtifacts = modelInfo.TrueAnimals;          % 1200x1  logical
Task          = modelInfo.ClassificationTask;   %         char 
                                                % (ani-all,art-all,ani-art)  
switch Task
case 'ani-all'
    Words  = [true(900,1);false(300,1)];        % 1200x1  logical
    Blocks = Blocks(~Ambiguous & Words,:);      %  809x10 logical
    y      = TrueAnimals(~Ambiguous & Words);   %  809x1  logical
    
case 'art-all'
    Words  = [true(900,1);false(300,1)];        % 1200x1  logical
    Blocks = Blocks(~Ambiguous & Words,:);      %  809x10 logical
    y      = TrueArtifacts(~Ambiguous & Words); %  809x1  logical

case 'ani-art'
    Words  = [true(240,1);false(960,1)];        % 1200x1  logical
    Blocks = Blocks(~Ambiguous & Words,:);      %  235x10 logical
    y      = TrueAnimals(~Ambiguous & Words);   %  235x1  logical
end

N          = length(nums);
cv_set     = Blocks(:,cv);
train_set  = ~cv_set;

[trainY,testY] = deal(cell(N,1));
[testY{:}]  = deal(y(cv_set));
[trainY{:}] = deal(y(train_set));
fprintf('cv blocks configured, and y made \n'); toc;

% Note: The data is already free of ambiguous words.
filename = 'FMRI_SMALL_OR_BLURRED.mat';
filepath = fullfile(DataDir,filename);
load(filepath);
%   Name                  Size                   Bytes  Class
% 
%   G                 18912x1                 10239656  cell                
%   G_oldspace        18912x1                 10239656  cell                
%   I                     1x1                        8  double              
%   J                     1x1                        8  double              
%   K                     1x1                        8  double              
%   imin                  1x1                        8  double              
%   jmin                  1x1                        8  double              
%   kmin                  1x1                        8  double              
%   true_nonzero          1x239158             1913264  double   

[trainX,testX] = deal(cell(N,1));
for i = 1:N
    filename = sprintf('Xs%d.mat',i);
    filepath = fullfile(DataDir,filename);
    load(filepath,'Xtemp'); 
    Xtemp = bsxfun(@minus,Xtemp,mean(Xtemp)); % Handles mean centering.
    testX{i} = Xtemp(cv_set,:);
    trainX{i} = Xtemp(train_set,:);
end
numvoxels = size(Xtemp,2);
numsamples = size(Xtemp,1);
clear Xtemp;

fprintf('data loaded (and mean centered) \n'); toc;
fprintf('training and test sets created \n')

%% algorithm parameters SHOULD BE MOVED TO THE modelInfo.
lamset = 2.^[-10:2];
lamerrs = zeros(length(lamset),1); %we will store the errors here

toc;

%% replicte the data matrix

[Xo, groups, group_arr] = makeA_multitask(trainX,G);
fprintf('matrix replicated \n'); toc;

%% iterate over all lambdas
for lamind = 1:length(lamset);

    lam = lamset(lamind);

    switch whatmethod
        case 1 %lasso
            if classify==1
                [Betahat,~,~] = Logistic_Lasso(trainX, trainY, lam);
            else
                [Betahat,~] = Least_Lasso(trainX, trainY, lam);
            end

        case 2
            if classify==1
                [Betahat,~,~] = Logistic_L21(trainX, trainY, lam);
            else
                [Betahat,~] = Least_L21(trainX, trainY, lam);
            end
        case 3
            % lam is the regularizer, rho = reg. for L1 term
            if classify==1
                [Betahat,~] = overlap_2stage(1,trainY,Xo,trainX,G,group_arr,groups,lam);
            else
                [Betahat,~] = overlap_2stage(0,trainY,Xo,trainX,G,group_arr,lam);
            end

        case 4
            if classify==1
                [Betahat,~] = overlap_1stage(1,trainY,Xo,trainX,G,group_arr,lam);
            else
                [Betahat,~] = overlap_1stage(0,trainY,Xo,trainX,G,group_arr,lam);
            end

        otherwise
            error('method not implemented \n')

    end
           
    fprintf('training done, lambda = %f \n', lam); toc;

    %cross validate
    a = cell2mat(testY);
    clear b;
    for person = 1:numpersons
        B = Betahat(:,person);
        x_for_err = testX{person};
        b(:,person) = sign(x_for_err*B);
    end
    b = b(:);
    L = length(a);err = 0;
    for l = 1:L
        err = err + (sign(a(l))~=sign(b(l)));
    end

    lamerrs(lamind) = err;

    fprintf(2,'validation done: lambda = %f \n', lam); toc
end
