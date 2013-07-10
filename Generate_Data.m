function [B,X] = Generate_Data(numsubjects, numvoxels, numtrials, G, k, alpha, deletion_level)

% generates data matrices X and associated sparse matrix B, such that 
% 1 = X_1B and -1 = X_2B
% INPUTS:
% numsubjects = number of columns in B
% numvoxels   = number of rows in B
% numtrials   = number of rows in X
% G           = groups
% k           = number of active groups
% alpha       = active group sparsity level
% deletion_l  = fraction of inactive columns that can be deleted from X

% form approx group sparse matrix
M = length(G);
activeg = randsample(M,k); % active groups
n = G{end}; n = n(end);
B = zeros(numvoxels,numsubjects);
for ii = 1:k
    g = G{activeg(ii)};
    B(g,:) = B(g,:) + 2*rand(length(g),numsubjects) - 1;
end
B = B.*(abs(B)>=alpha); 
% coefficient matrix created

% create a DATA matrix X for each subject, based on his sparsity pattern
X = cell(0);
for j = 1:numsubjects
    
    y = ones(numtrials,1);
    y = [y ; -y];
    b = B(:,j);
    
    for l = 1:2*numtrials
        loc = find(b~=0);
        if ~isempty(loc)
            loc = loc(1);
            
            
            temp = rand(1,numvoxels-1);
            btemp = b; btemp(loc) = [];
            atemp = (y(l) - temp*btemp)/b(loc);
            Xtemp(l,:) = [temp(1:loc-1) atemp temp(loc:end)];
        else
            Xtemp(l,:) = rand(1,numvoxels);
        end
    end
    X = [X ; {Xtemp}];
    
end

if deletion_level>0
    for j =1:numsubjects
        A = X{j};
        idxj = find(B(:,j)~=0);
        availables = setdiff([1:n],idxj);
        to_del = randsample(availables, floor(deletion_level*length(availables)));
        A(:,to_del) = zeros(2*numtrials,length(to_del));
        X(j) = {A};
    end
end   

end
    
        
    
    