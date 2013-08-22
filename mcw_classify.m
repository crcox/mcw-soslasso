function [final_test,final_train,cv_test,cv_train] = mcw_classify(Y,X,GroupInfo,metadata,params)
% MCW_CLASSIFY Solves optimization problem to discover active voxels
%
%  USAGE:
%    [...] = mcw_classify(Y,X,GroupInfo,metadata,params)
%
%  INPUTS (all are required):
%    Y Cell array of target vectors (1 or -1) for each task (i.e. subject).
%    X Cell array of prepped functional data for each task.
%    GroupInfo Structure containing information for managing groups.
%    metadata See MCW_SCRIPT_TEMPLATE
%    params See MCW_SCRIPT_TEMPLATE
%
%  OUTPUTS (optional):
%    final_test Final model result on hold out set.
%    final_train Performance on the final training set.
%    cv_test Cross validated performance during lambda selection.
%    cv_train Training performance during lambda selection.
%
%  See also:
%    MCW_HELP
%    MCW_PREP_DATA
%    MCW_DEFINE_REPSPACE
%    MCW_SCRIPT_TEMPLATE
%
% Nikhil Rao and Chris Cox | University of Wisconsin-Madison | 2013-08-06

lamset     = params.lamset;
numcvs     = params.numcvs;
whatmethod = params.whatmethod;
classify   = params.classify;
numpersons = length(X);
numvoxels  = size(X{1},2); % Number of voxels in each cell must be equal.
numsamples = length(Y{1});
numlambda  = length(lamset);

G         = GroupInfo.G;
groups    = GroupInfo.groups;
group_arr = GroupInfo.group_arr;
RepIndex  = GroupInfo.RepIndex;

CVBlocks       = metadata.CVBlocks; 
FinalTestBlock = metadata.FinalTestBlock; 

%% Mean-center the data
if params.MeanCenter == true
    for ii = 1:numpersons
        X{ii} = bsxfun(@minus,X{ii},mean(X{ii}));
    end
    fprintf('DATA CENTERED\n')
end

%% Normalize variance
if params.NormVariance == true
    for ii = 1:numpersons
        X{ii} = bsxfun(@rdivide,X{ii},std(X{ii}));
    end
    fprintf('DATA CENTERED\n')
end


%% Cross Validation Module
if params.RecoveryMode > 0
    load('recovery.mat','Betahat','C');

else
    Betahat = zeros(numvoxels,numpersons*numcvs*numlambda);
    C = zeros(1,numpersons*numcvs*numlambda);
    ix     = uint32(0);
    
    for lamind = 1:numlambda
        lam = lamset(lamind);
        for cv = 1:numcvs;
            % Compute indexes for storing the betahats
            a  = uint32(ix * numpersons + 1);
            ix = uint32(1 + ix);
            b  = uint32(ix * numpersons); 
            
            CVTestBlock = CVBlocks(:,cv);
            train = ~CVTestBlock & ~FinalTestBlock;
            
            [trainX, trainY] = deal(cell(numpersons,1));
            for person = 1:numpersons
                trainX{person} = X{person}(train,:);
                trainY{person} = Y{person}(train);
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
                    % lam is the regularizer, rho = reg. for L1 term
                    try
                        if classify==1
                            [tempB,tempC] = overlap_2stage(1,trainY,trainX,G,RepIndex,group_arr,groups, lam);
                            Betahat(:,a:b) = tempB;
                            C(a:b) = tempC; 
                        else
                            [Betahat,~] = overlap_2stage(0,trainY,trainX,G,RepIndex,group_arr,lam);
                        end
                    catch ME
                        save('recovery.mat','Betahat','C');
                        rethrow(ME);
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
            fprintf(2,'%d CV Run done \n',cv);
        end
    end
    Betahat = sparse(Betahat);
    save('recovery.mat','Betahat','C');
end

%% Compute D-Prime for all CV runs.
% Betahat is subjects * cv * lambda, populated in that order.
N = numpersons*numcvs*numlambda;
scores = zeros(numsamples,N);

for i = 1:numpersons
    scores(:,i:numpersons:N) = bsxfun(@plus,X{i} * Betahat(:,i:numpersons:N), C(1:numpersons:N));
end

prediction = reshape(scores>0,numsamples*numpersons,numcvs*numlambda);
truth = cell2mat(Y) > 0;
TEST  = repmat(CVBlocks(:,1:numcvs),numpersons,numlambda);
TRAIN = ~TEST & ~repmat(FinalTestBlock,numpersons,numcvs*numlambda);

cv_test.HIT = sum(bsxfun(@and,bsxfun(@and, truth, prediction),TEST));
cv_test.HR  = cv_test.HIT ./ sum(TEST(truth,:));
cv_test.FA  = sum(bsxfun(@and,bsxfun(@and, ~truth, prediction),TEST));
cv_test.FAR = cv_test.HIT ./ sum(TEST(~truth,:));
cv_test.FAR(cv_test.FAR==0) = .0005;cv_test.FAR(cv_test.FAR==1) = .9995;
cv_test.HR(cv_test.HR==0)   = .0005;cv_test.HR(cv_test.HR==1)   = .9995;
cv_test.DPrime = norminv(cv_test.HR) - norminv(cv_test.FAR);

cv_train.HIT = sum(bsxfun(@and,bsxfun(@and, truth, prediction),TRAIN));
cv_train.HR  = cv_train.HIT ./ sum(TRAIN(truth,:));
cv_train.FA  = sum(bsxfun(@and,bsxfun(@and, ~truth, prediction),TRAIN));
cv_train.FAR = cv_train.HIT ./ sum(TRAIN(~truth,:));
cv_train.FAR(cv_train.FAR==0) = .0005;cv_train.FAR(cv_train.FAR==1) = .9995;
cv_train.HR(cv_train.HR==0)   = .0005;cv_train.HR(cv_train.HR==1)   = .9995;
cv_train.DPrime = norminv(cv_train.HR) - norminv(cv_train.FAR);

notes = ' Counts and DPrime are over all subjects at once. \n Betahat is numpersons*numcvs*numlambda, in that order.\n';
save('DIAGNOSTICS.mat','cv_test','cv_train','Betahat','numpersons','numcvs','lamset','notes');
fprintf('\n ALL CV DONE \n')

%% PICK THE BEST REGULARIZATION PARAMETER AND RELEARN MODEL
[~,lamax] = max(mean(reshape(cv_test.DPrime,numcvs,numlambda)));
lam = lamset(lamax);  % pick the lambda that minimzes the CV error
test = FinalTestBlock;
train = ~test;
trainX = cellfun(@(x) x(train,:),X,'Unif',0);
trainY = cellfun(@(y) y(train,:),Y,'Unif',0);

if params.RecoveryMode > 1
    load('recovery_FINAL.mat','Betahat','C');

else
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
                [Betahat,C] = overlap_2stage(1,trainY,trainX,G,RepIndex,group_arr,groups, lam);
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
    Betahat = sparse(Betahat);
    save('recovery_FINAL.mat','Betahat','C');
end

%% TEST MODEL PERFORMANCE ON TEST SET
scores = zeros(numsamples,numpersons);
for i = 1:numpersons
    scores(:,i) = bsxfun( @plus,X{i} * Betahat(:,i),C );
end
prediction = scores(:)>0;
truth = cell2mat(Y) > 0;
TEST  = repmat(test,numpersons,1);
TRAIN = repmat(train,numpersons,1);
final_test.HIT = sum(bsxfun(@and,bsxfun(@and, truth, prediction),TEST));
final_test.HR  = final_test.HIT ./ sum(TEST(truth,:));
final_test.FA  = sum(bsxfun(@and,bsxfun(@and, ~truth, prediction),TEST));
final_test.FAR = final_test.HIT ./ sum(TEST(~truth,:));
final_test.FAR(final_test.FAR==0) = .0005;final_test.FAR(final_test.FAR==1) = .9995;
final_test.HR(final_test.HR==0)   = .0005;final_test.HR(final_test.HR==1)   = .9995;
final_test.DPrime = norminv(final_test.HR) - norminv(final_test.FAR);

final_train.HIT = sum(bsxfun(@and,bsxfun(@and, truth, prediction),TRAIN));
final_train.HR  = final_train.HIT ./ sum(TRAIN(truth,:));
final_train.FA  = sum(bsxfun(@and,bsxfun(@and, ~truth, prediction),TRAIN));
final_train.FAR = final_train.HIT ./ sum(TRAIN(~truth,:));
final_train.FAR(final_train.FAR==0) = .0005;final_train.FAR(final_train.FAR==1) = .9995;
final_train.HR(final_train.HR==0)   = .0005;final_train.HR(final_train.HR==1)   = .9995;
final_train.DPrime = norminv(final_train.HR) - norminv(final_train.FAR);

c = clock; 
fmt = '%s_%s_%s_%s.mat';
if classify==1
    filename_type = 'classify';
else
    filename_type = 'regress';
end
filename_script = 'SURE';
filename_timestamp = sprintf('%d-%02d-%02d_%02d:%02d',c(1:5));

switch whatmethod
    case 1 %lasso
        filename_method = 'LASSO';
        filename = sprintf(fmt,filename_type, filename_script, filename_method, filename_timestamp);
        
    case 2 %glasso
        filename_method = 'GLASSO';
        filename = sprintf(fmt,filename_type, filename_script, filename_method, filename_timestamp);
        
    case 3 %soglasso
        filename_method = 'SOSLASSO';
        filename = sprintf(fmt,filename_type, filename_script, filename_method, filename_timestamp);
        
    case 4 %oglasso
        filename_method = 'OGLASSO';
        filename = sprintf(fmt,filename_type, filename_script, filename_method, filename_timestamp);     
end
save(filename,'Betahat','X','Y','final_train','final_test');
fprintf('DATA SAVED \n')

end
