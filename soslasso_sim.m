clear;clc;
%% Specify Simulation environment
% DataType = 'Sparse Group';
% DataType = 'Identical Group';
DataType = 'No Group';

P = uint16(16);     % people
N = uint16(1024);   % voxels
M = uint16(64);     % group size
L = uint16(4);      % group space

A = 1:L:(N-M+1);    % group start ind
B = (M):L:N;        % group end ind
K = uint16(length(A));
G = cell(1,K);
for i=1:K
    G{i} = A(i):B(i);
end

J = uint16(64);     % K/J groups will be active
I = idivide(K,J,'floor');

active_groups = uint16(randperm(K,I));

g = cell2mat(G(active_groups));

H = uint16(J/2);     % V/H voxels will be active.
g_set = unique(g);
V = uint16(length(g_set));
W = idivide(V,H,'floor');

%% Define Data
X = zeros(P,N);
active_voxels = zeros(P,W);

switch DataType
    case 'Sparse Group'    
        for i=1:P;
            temp = uint16(randperm(V,W));
            active_voxels(i,:) = g_set(temp);
            X(i,active_voxels(i,:)) = 1;
        end
        
    case 'No Group'
        for i=1:P
            temp = uint16(randperm(N,W));
            active_voxels(i,:) = temp;
            X(i,active_voxels(i,:)) = 1;
        end
        
    case 'Identical Group'
        temp = uint16(randperm(V,W));
        for i=1:P
            active_voxels(i,:) = temp;
            X(i,active_voxels(i,:)) = 1;
        end  
end
imagesc(X)





% In citing MALSAR in your papers, you can use the following:
% [Zhou 2012] J. Zhou, J. Chen and J. Ye. MALSAR: Multi-tAsk Learning via 
%     StructurAl Regularization. Arizona State University, 2012. 
%     http://www.public.asu.edu/~jye02/Software/MALSAR.
% 
% If you use LaTeX, you can download the BibTex entry: 
%     http://www.public.asu.edu/~jye02/Software/MALSAR/citation.bib.
% 
% The MALSAR software project has been supported by research grants from 
% the National Science Foundation (NSF). 