% script to perform ONE cv run
clear;clc;
%% specify which cv run to perform
cv = 5;
disp(cv);
classify = 1;
whatmethod = 3;
if (cv < 1 || cv > 9)
    error('cv should be an integer from 1-9 \n');
end
tic;
%% load the data
load FMRI_SMALL_OR_BLURRED;
X = cell(0);
nums = [314 1407 1577 2177 2184 2199 2200 2202 2227,...
    2232 2347 3283 3303];
X = cell(0);
for i = 1:length(nums)
    str = strcat('load Xs',num2str(i),' ;');
    eval(str);
    X = [X; {Xtemp}];
end
fprintf('data loaded \n'); toc;
%% load the metadata, and extract the indices for this cv run
load mcw_metadata
ambwrds = find(Ambiguous(1:900)==0);
Blocks = Blocks(ambwrds,:);
shuffle = [];
setsizes = zeros(10,1);
for i = 1:10;
    inds = find(Blocks(:,i)==1);
    shuffle = [shuffle inds'];
    setsizes(i) = length(inds);
end
startind = sum(setsizes(1:cv))+1;
cvind = shuffle(startind:startind+setsizes(cv+1)-1);
trainind = setdiff(shuffle(setsizes(1)+1:end),cvind);
if ~isempty(intersect(cvind, trainind))
   error('CV set and Training set are incompatible \n') 
end

%% extract a data matrix to get parameters
X1 = X{1};
numpersons = length(X);
numvoxels = size(X1,2);
numsamples = size(X1,1);
clear X1;

%% create y and center the data
y = [ones(94,1) ; -ones(numsamples-94,1)]; %the first 94 points correspond to animals
Y = cell(0);
for ii = 1:numpersons
    Xtemp = X{ii};
    mx = mean(Xtemp);
    mx = repmat(mx,size(Xtemp,1),1);
    Xtemp = Xtemp - mx;
    X{ii} = Xtemp;
    Y = [Y; {[y]}];
end
fprintf('data centered and y made \n'); toc;
%% algorithm parameters
lamset = 2.^[-10:2];
lamerrs = zeros(length(lamset),1); %we will store the errors here

%% create training and test sets
trainX = cell(0);trainY = cell(0);testX = cell(0); testY = cell(0);
for person = 1:numpersons
    Xtemp = X{person};
    Ytemp = Y{person};
    Xtempt = Xtemp(trainind,:); %training data
    Xtempc = Xtemp(shuffle(cvind),:);    % CV data
    Ytempt = Ytemp(trainind,:);
    Ytempc = Ytemp(shuffle(cvind),:);
    trainX = [trainX ; {[Xtempt]}];
    testX  = [testX ; {[Xtempc]}];
    trainY = [trainY ; {[Ytempt]}];
    testY  = [testY ; {[Ytempc]}];
end
fprintf('training and test sets created \n')
toc;

%% replicte the data matrix
[Xo, groups, group_arr] = makeA_multitask(trainX,G);
fprintf('matrix replicated \n'); toc;

%% iterate over all lambdas
lamind = 0;
for lam = lamset
    lamind = 1+lamind;
    
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
    y_for_err = repmat(Ytempc,1,numpersons);
    
    a = y_for_err(:);
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
%% save final errors for later
str = strcat('save classify_ANIMAL_',num2str(whatmethod),'_',num2str(cv),' lamerrs');
eval(str);
fprintf('DATA SAVED \n'); toc;

