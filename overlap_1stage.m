function [Xhat,C] = overlap_1stage(loss,Y,Xo,X,G,group_arr,lambda)

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
    [W, C, ~] = Logistic_L21_1stage(Xo, Y, lambda, group_arr);
elseif loss == 0
    [W, ~] = Least_L21_1stage(Xo, Y, lambda, group_arr);
    C = [];
else
    error('loss has to be 0 or 1 \n');
end
% W is the output matrix
%we now need to combine overlapping groups

n = G{end};
n = n(end);
T = length(Y);
Xhat = zeros(n,T);
for ii = 1:length(G)
    t = G{ii};
    s = group_arr(ii,:);
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




