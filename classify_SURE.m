% script to test animals vs all

clear
clc;
close all;

%load metadata
load('mcw_metadata.mat');

% load the data
load('GROUPS.size21.offset10.mat');
nums = [314 1407 1577 2177 2184 2199 2200 2202 2227,...
    2232 3283 3303];
X = cell(length(nums),1);
for i = 1:length(nums)
    filename = sprintf('%d.aseg.tlrc.mat',nums(i));
    load(filename,'funcData_tlrc');
    X{i} = funcData_tlrc;
end

% form chunks for training and testing

Words = [true(240,1); false(960,1)] & ~Ambiguous;
% to make chunks, take the same chunks as in Blocks but discard ambiguous
Blocks = Blocks(Words,:);
Blocks5 = false(sum(Words),5);
for i = 1:5
    Blocks5(:,i) = any([Blocks(:,i),Blocks(:,i+5)],2);
end
setsizes = sum(Blocks5);

rand_order = [];
pieces = 5;
order = [2 3 4 5 1];
numcvs = pieces - 1;
for j = 1:length(order)
    i = order(j);
    vv = Blocks(:,2*i-1 : 2*i); vv = sum(vv,2);
    inds = find(vv==1);
    rand_order = [rand_order inds'];
end

% MAKE SURE THERE IS NO VESTIGIAL EMPTY GROUP
g = G{end};
if isempty(g)
    G = G(1:end-1);
end
fprintf('FMRI DATA LOADED \n');

% extract a data matrix to get parameters
X1 = X{1};
numpersons = length(X);
numvoxels = size(X1,2);
numsamples = size(X1,1);
clear X1;

whatmethod = 3 ; % 1 = lasso, 2 = glasso, 3 = soslasso, 4 = oglasso
classify = 1;
flip = 0;

% center the data
numanimals = sum(TrueAnimals(Words));
numartifacts = sum(TrueArtifacts(Words));
y = [ones(numanimals,1) ; -ones(numsamples-numanimals,1)]; %the first 94 points correspond to animals
Y = cell(length(nums),1);
Y(:) = deal({y});
for ii = 1:numpersons
    Xtemp = X{ii};
    mx = mean(Xtemp);
    mx = repmat(mx,size(Xtemp,1),1);
%     sx = std(Xtemp);
%     sx = repmat(sx,size(Xtemp,1),1);
    Xtemp = (Xtemp - mx);
    X{ii} = Xtemp;
end

fprintf('DATA CENTERED AND CELLS CREATED \n')

%% problem parameters for full data (total size = 240)

% lamset = 2.^[-10:10]; % regularizer values
% lamset = linspace(2^8,2^12, 100);
lamset = linspace(500,1000,10);
% rhoset = 1;  % multiplier for the L1 case. 1 corresponds to ONE regularizer

testinds = rand_order(1:setsizes(1));
Beta_CV = cell(0);
endind = setsizes(1);
[TRAIN_ERROR,CV_ERROR,CV_SPARSE] = deal(zeros(numcvs,length(lamset)));
BETAHAT = cell(numcvs,length(lamset));

%% cross validation module
for cv = 1:numcvs  
    % make training and hold out X and Y;
    cvind = rand_order(endind+1:endind+setsizes(cv+1));
    trainind = setdiff(rand_order(length(testinds)+1:numsamples),rand_order(cvind));
    scv = randsample(length(cvind),length(cvind));
    cvind = cvind(scv);
    str = randsample(length(trainind),length(trainind));
    trainind = trainind(str);
    trainX = cell(size(X));trainY = cell(size(Y));testX = cell(size(X)); testY = cell(size(Y));
    for person = 1:numpersons
        trainX{person} = X{person}(trainind,:);
        testX{person}  = X{person}(rand_order(cvind),:);
        trainY{person} = Y{person}(trainind);
        testY{person}  = Y{person}(rand_order(cvind));
    end
    
    % run the method over all lambda values
%    [RepIndex, groups, group_arr] = definerepspace(G);
%    [~, groups, group_arr] = makeA_multitask(X,G);
    
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
                        [Betahat,~] = overlap_2stage(1,trainY,trainX,G,RepIndex,group_arr,groups, lam);
                    else
                        [Betahat,~] = overlap_2stage(0,trainY,trainX,G,RepIndex,group_arr,lam);
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
            
            fprintf('.')
            
            %cross validate
            y_for_err = cell2mat(trainY');
            
            a = y_for_err(:);
            clear b
            for person = 1:numpersons
                B = Betahat(:,person);
                x_for_err = trainX{person};
                b(:,person) = sign(x_for_err*B);
            end
            b = b(:);
            L = length(a);err = 0;
            for l = 1:L
                err = err + (sign(a(l))~=sign(b(l)));
            end
            
            TRAIN_ERROR(cv,lamind) = err;
            
            y_for_err = cell2mat(testY');
            
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

            CV_ERROR(cv,lamind) = err;
            CV_SPARSE(cv,lamind) = nnz(Betahat);
            
            fprintf(2,'%d CV Run done \n',cv);
            BETAHAT{cv,lamind} = sparse(Betahat);
    end
    endind = setsizes(cv+1)+endind;
end
save('DIAGNOSTICS.mat','CV_ERROR','CV_SPARSE','TRAIN_ERROR','BETAHAT');
fprintf('\n ALL CV DONE \n')
ALLERRS = CV_ERROR;
if numcvs>1
    CV_ERROR = mean(CV_ERROR);
end

%% PICK THE BEST REGULARIZATION PARAMETER AND RELEARN MODEL
[value,lamin] = min(CV_ERROR);
lam = lamset(lamin);  % pick the lambda that minimzes the CV error
trainind = rand_order(1+length(testinds):numsamples);  %entire training set
trainX = cell(0);trainY = cell(0);
for person = 1:numpersons
    Xtemp = X{person};
    Ytemp = Y{person};
    Xtempt = Xtemp(trainind,:);
    Ytempt = Ytemp(trainind,:);
    trainX = [trainX ; {[Xtempt]}];
    trainY = [trainY ; {[Ytempt]}];
end

if ((whatmethod==3)||(whatmethod==4))
    [RepIndex, groups, group_arr] = definerepspace(G);
end
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
        if classify==1
            [Betahat,~] = overlap_2stage(1,trainY,trainX,G,RepIndex,group_arr,groups, lam);
        else
            [Betahat,~] = overlap_2stage(0,trainY,trainX,G,RepIndex,group_arr,lam);
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
% Betahat is the final model


%% TEST MODEL PERFORMANCE ON TEST SET

finX = cell(0); finY = cell(0);
for person = 1:numpersons
    Xtemp = X{person};
    Ytemp = Y{person};
    Xtemp = Xtemp(testinds,:);
    Ytemp = Ytemp(testinds,:);
    finX = [finX ; {[Xtemp]}];
    finY = [finY ; {[Ytemp]}];
end
y_for_err = repmat(Ytemp,1,numpersons);
a = y_for_err(:);
clear b;
for person = 1:numpersons
    B = Betahat(:,person);
    x_for_err = finX{person};
    b(:,person) = sign(x_for_err*B);
end
b = b(:);
L = length(a);err = 0;
for l = 1:L
    err = err + (sign(a(l))~=sign(b(l)));
end %this is the overall error
%%
c = clock;
timestamp = sprintf('%d-%02d-%02d_%02d:%02d',c(1:5));
switch whatmethod
    case 1 %lasso
        if classify==1
            str = strcat('save classify_SURE_LASSO_',timestamp,' Betahat X Y rand_order err ALLERRS');
            eval(str);
        else
            str = strcat('save regress_SURE_LASSO_',timestamp,' Betahat X Y rand_order err ALLERRS ');
            eval(str);
        end
    case 2 %glasso
        if classify==1
            str = strcat('save classify_SURE_GLASSO_',timestamp,' Betahat X Y rand_order err ALLERRS');
            eval(str);
        else
            str = strcat('save regress_SURE_GLASSO_',timestamp,' Betahat X Y rand_order err ALLERRS');
            eval(str);
        end
    case 3 %soglasso
        if classify==1
            str = strcat('save classify_SURE_SOSLASSO_',timestamp,' Betahat X Y rand_order err ALLERRS');
            eval(str);
        else
            str = strcat('save regress_SURE_SOSLASSO_',timestamp,' Betahat X Y rand_order err ALLERRS');
            eval(str);
        end
    case 4 %oglasso
        if classify==1
            str = strcat('save classify_SURE_OGLASSO_',timestamp,' Betahat X Y rand_order err ALLERRS');
            eval(str);
        else
            str = strcat('save regress_SURE_OGLASSO_',timestamp,' Betahat X Y rand_order err ALLERRS');
            eval(str);
        end
end
fprintf('DATA SAVED \n')
