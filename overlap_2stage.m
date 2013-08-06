function [Xhat,C] = overlap_2stage(loss,Y,X,G,RepIndex,group_arr,groups,lambda)

% function to perform the overlapping SGL optimization, with both LS and LOGIT loss
% INPUTS
% loss    = type of loss function. 0 = least squares, 1 = logistic
% Y       = T X 1 cell array of observations . T  =number of tasks
% X       = T X 1 cell array of data
% Xo      = T X 1 cell array of replicated data matrices
% G       = cell array of groups
% group_arr = output from replication step
% lambda  = regularizer
% OUTPUTS
% Xhat    = debiased output
% C       = the bias term in the output
% 
% CODE REQUIRES MALSAR PACKAGE
%
% Nikhil Rao
% 3/17/13
if loss == 1
    [W, C, ~] = Logistic_L21_2stage_1reg(X, Y, lambda, RepIndex, group_arr, groups);
elseif loss == 0
    [W, ~] = Least_L21_2stage(Xo, Y, lambda, group_arr);
    C = [];
else
    error('loss has to be 0 or 1 \n');
end
% W is the output matrix
%we now need to combine overlapping groups

for iii = 1:length(G)
    temp = G{iii};
    n(iii) = max(temp);
end
n = max(n);
T = length(Y);
Xhat = zeros(n,T);
% identify whether a dummy variable exists and chuck it
dummy = max(max(group_arr));
mask = (group_arr == dummy);
isdummy = 0;
if sum(sum(mask))>1
    isdummy = 1;
end
for ii = 1:length(G)
    t = G{ii};
    s = group_arr(ii,:);
    if isdummy == 1
       dummyind = find(s == dummy);
       s(dummyind) = [];
    end
    Xhat(t,:) = Xhat(t,:) + W(s,:);
end

%debias the solution
temp = Xhat;
for ii = 1:T
    idx = find(temp(:,ii)~=0);
    if ~isempty(idx)
        Xtemp = X{ii};
        Xtemp = Xtemp(:,idx);
        vtemp = Xtemp\Y{ii};
        Xhat(idx,ii) = vtemp;
    end
end
    
end




