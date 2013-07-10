%%% shrinkage on product norm

function y = SGL_overlap_shrink(X,G,groups,lam1,lam2)

y = zeros(size(X));
% identify whether a dummy variable exists and chuck it
dummy = max(max(G));
mask = G == dummy;
isdummy = 0;
if sum(sum(mask))>1
    isdummy = 1;
end
% BEGIN EFFICIENT CODE

% step 1: perform soft thresholding
X_soft = sign(X).*max(abs(X) - lam1,0);

%step 2: perform group soft thresholding
M = size(G,1); % number of groups
X_soft = [X_soft; zeros(1,size(X_soft,2))]; % for the dummy
Xtemp = sum(X_soft.^2,2); %xtemp is now a vector
Xtemp = sum(Xtemp(G),2);
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


% for ii = 1:size(G,1)
%     t = G(ii,:);
%     
%     if isdummy == 1
%        dummyind = find(t == dummy);
%        t(dummyind) = [];
%     end
%     x = X(t,:);
%     x1 = x(:);
%     n2 = norm(x1);
%     n1 = norm(x1,1);
%        
%     %shrink
%     xtemp = sign(x).*max(0,abs(x) - lam1);
%     nx = norm(xtemp(:));
%     if nx~=0
%         y(t,:) = xtemp/nx.*max(0,nx - lam2);
%     else
%         y(t,:) = 0;
%     end
% end


end