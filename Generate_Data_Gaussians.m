function [B,X,Y] = Generate_Data_Gaussians(numtasks, numvariable, nummeasure, G, k, alpha, type)

% generates data matrices X,Y and associated sparse matrix B, such that 
% Y{i} = XB(:,i)
% INPUTS:
% numtasks    = number of columns in B
% numvariable = number of rows in B
% nummeasure  = number of rows in X
% G           = groups
% k           = number of active groups
% alpha       = active group sparsity level
% type        = Logistic if 1, else Least Squares

% form approx group sparse matrix
M = length(G);
activeg = randsample(M,k); % active groups
n = G{end}; n = n(end);
B = zeros(numvariable,numtasks);
for ii = 1:k
    g = G{activeg(ii)};
    B(g,:) = B(g,:) + 2*rand(length(g),numtasks) - 1;
end
B = B.*(abs(B)>=alpha); 
% coefficient matrix created

% create a DATA matrix X for each subject, based on his sparsity pattern
X = cell(0);
Y = cell(0);
for j = 1:numtasks
    
    Xtemp = randn(nummeasure,numvariable)/sqrt(nummeasure);
    X = [X ; {Xtemp}];
    Ytemp = Xtemp*B(:,j);
    
    if type==1
        Ytemp = sign(Ytemp);
    end
    Y = [Y ; {Ytemp}];
end

end
    
        
    
    