%%% shrinkage on product norm

function y = SGL_overlap_shrink(X,group_arr,groups,lam1,lam2)

y = zeros(size(X));
% identify whether a dummy variable exists and chuck it
dummy = max(max(group_arr));
mask = group_arr == dummy;
isdummy = false;
if sum(sum(mask))>1
    isdummy = true;
end
% BEGIN EFFICIENT CODE

% step 1: perform soft thresholding
X_soft = sign(X).*max(abs(X) - lam1,0);

%step 2: perform group soft thresholding
M = size(group_arr,1); % number of groups
X_soft = [X_soft; zeros(1,size(X_soft,2))]; % for the dummy
Xtemp = sum(X_soft.^2,2); %xtemp is now a vector
Xtemp = sum(Xtemp(group_arr),2);
Xtemp = sqrt(Xtemp);
Xtemp = max(Xtemp - lam2,0); % this is the multiplying factor
Xtemp = Xtemp./(Xtemp + lam2);
if (size(Xtemp,1)~=M)
    error('something weird is happening with the group shrinkage \n');
end
Xtemp = Xtemp(groups);
Xtemp = repmat(Xtemp,1,size(X,2));
y = X_soft(1:end-1,:).*Xtemp;

% END EFFICIENT CODE


end