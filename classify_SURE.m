% script to test animals vs all

clear
clc;
close all;


% load the data
load FMRI_SURE_SMALL_OLDBG
X = cell(0);
nums = [314 1407 1577 2177 2184 2199 2200 2202 2227,...
    2232 2347 3283 3303];
X = cell(0);
for i = 1:length(nums)
    str = strcat('load Xolds',num2str(i),' ;');
    eval(str);
    X = [X; {Xtemp}];
end

%load metadata and form chunks for training and testing
load mcw_metadata
ambwrds = find(Ambiguous(1:240)==0);
% to make chunks, take the same chunks as in Blocks but discard ambiguous
Blocks = Blocks(1:240,:);
shuffle = [];
setsizes = zeros(5,1);
pieces = 5;
order = [2 3 4 5 1];
numcvs = pieces - 1;
for j = 1:length(order)
    i = order(j);
    vv = Blocks(:,2*i-1 : 2*i); vv = sum(vv,2);
   inds = find(vv==1);
   shuffle = [shuffle inds'];
   setsizes(j) = length(inds);
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
numanimals = 94;
y = [ones(numanimals,1) ; -ones(numsamples-numanimals,1)]; %the first 94 points correspond to animals
Y = cell(0);
for ii = 1:numpersons
    Xtemp = X{ii};
    mx = mean(Xtemp);
    mx = repmat(mx,size(Xtemp,1),1);
%     sx = std(Xtemp);
%     sx = repmat(sx,size(Xtemp,1),1);
    Xtemp = (Xtemp - mx);
    X{ii} = Xtemp;
    Y = [Y; {[y]}];
end

fprintf('DATA CENTERED AND CELLS CREATED \n')

%% problem parameters for full data (total size = 240)

% lamset = 2.^[-10:10]; % regularizer values
lamset = linspace(2^8,2^12, 100);
% rhoset = 1;  % multiplier for the L1 case. 1 corresponds to ONE regularizer

testinds = shuffle(1:setsizes(1));
Beta_CV = cell(0);
endind = setsizes(1);
%% cross validation module
for cv = 1:numcvs
    
    % make training and hold out X and Y;
    cvind = shuffle(endind+1:endind+setsizes(cv+1));
    trainind = setdiff(shuffle(length(testinds)+1:numsamples),shuffle(cvind));
    scv = randsample(length(cvind),length(cvind));
    cvind = cvind(scv);
    str = randsample(length(trainind),length(trainind));
    trainind = trainind(str);
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
    
    % run the method over all lambda values
    if ((whatmethod==3)||(whatmethod==4))
        [~, groups, group_arr] = makeA_multitask(trainX,G);
    end
    y_for_err = repmat(Ytempt,1,numpersons);
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
                        [Betahat,~] = overlap_2stage(1,trainY,trainX,G,group_arr,groups, lam);
                    else
                        [Betahat,~] = overlap_2stage(0,trainY,trainX,G,group_arr,lam);
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
            y_for_err = repmat(Ytempt,1,numpersons);
            
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

            CV_ERROR(cv,lamind) = err;
            CV_SPARS(cv,lamind) = nnz(Betahat);
            
            fprintf(2,'%d CV Run done \n',cv);
            
            save DIAGNOSTICS CV_ERROR CV_SPARS TRAIN_ERROR
    end
    endind = setsizes(cv+1)+endind;
end
fprintf('\n ALL CV DONE \n')
ALLERRS = CV_ERROR;
if numcvs>1
    CV_ERROR = mean(CV_ERROR);
end

%% PICK THE BEST REGULARIZATION PARAMETER AND RELEARN MODEL
[value,lamin] = min(CV_ERROR);
lam = lamset(lamin);  % pick the lambda that minimzes the CV error
trainind = shuffle(1+length(testinds):numsamples);  %entire training set
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
    [Xo, groups, group_arr] = makeA_multitask(trainX,G);
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
            [Betahat,~] = overlap_2stage(1,trainY,trainX,G,group_arr,groups, lam);
        else
            [Betahat,~] = overlap_2stage(0,trainY,trainX,G,group_arr,lam);
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
c = c(2:end-1); % month day hour min
strin = strcat(num2str(c(1)),num2str(c(2)),'_',num2str(c(3)),num2str(c(4)));
switch whatmethod
    case 1 %lasso
        if classify==1
            str = strcat('save classify_SURE_LASSO_',strin,' Betahat X Y shuffle err ALLERRS');
            eval(str);
        else
            str = strcat('save regress_SURE_LASSO_',strin,' Betahat X Y shuffle err ALLERRS ');
            eval(str);
        end
    case 2 %glasso
        if classify==1
            str = strcat('save classify_SURE_GLASSO_',strin,' Betahat X Y shuffle err ALLERRS');
            eval(str);
        else
            str = strcat('save regress_SURE_GLASSO_',strin,' Betahat X Y shuffle err ALLERRS');
            eval(str);
        end
    case 3 %soglasso
        if classify==1
            str = strcat('save classify_SURE_SOSLASSO_',strin,' Betahat X Y shuffle err ALLERRS');
            eval(str);
        else
            str = strcat('save regress_SURE_SOSLASSO_',strin,' Betahat X Y shuffle err ALLERRS');
            eval(str);
        end
    case 4 %oglasso
        if classify==1
            str = strcat('save classify_SURE_OGLASSO_',strin,' Betahat X Y shuffle err ALLERRS');
            eval(str);
        else
            str = strcat('save regress_SURE_OGLASSO_',strin,' Betahat X Y shuffle err ALLERRS');
            eval(str);
        end
end
fprintf('DATA SAVED \n')