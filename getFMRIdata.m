
%% FIND RANGE OF I J K
clear; clc; close all
nums = [314 1225 1407 1577 2177 2184 2199 2200 2202 2227,...
    2232 2347 3283 3303];

I = []; J = []; K = [];icount = 0; jcount = 0; kcount = 0;
for jj = 1:length(nums)
    
    disp(jj);
    str = strcat('load s',num2str(nums(jj)),'.aseg_sort.mat');
    eval(str);
    
    %check if there are zeros
    if ~isempty(find(ijk(:,1)==0))
        icount = 1+icount;
    end
    if ~isempty(find(ijk(:,2)==0))
        jcount = 1+jcount;
    end
    if ~isempty(find(ijk(:,3)==0))
        kcount = 1+kcount;
    end
    
    ijk = ijk + 1;
    
    I = [I; max(ijk(:,1))];
    J = [J; max(ijk(:,2))];
    K = [K; max(ijk(:,3))];
    
end
I = max(I);J = max(J);K = max(K);

%% MAPPING FROM ijk ---> n = i + I(j-1) + IJ(k-1)

% MAKE NEW DATA MATRICES 
X = cell(0);
Xsure = cell(0);
for jj = 1:length(nums)
    disp(jj);
    str = strcat('load s',num2str(nums(jj)),'.aseg_sort.mat');
    eval(str);
    ijk = 1+ijk;
    numvox = size(funcData,2);
    Xtemp = sparse(1200,I*J*K);
    for ii = 1:numvox
        i = ijk(ii,1); j = ijk(ii,2); k = ijk(ii,3);
        n = i + I*(j-1) + I*J*(k-1);
        Xtemp(:,n) = funcData(:,ii);
        
    end
    XX = Xtemp(1:240,:);
    X = [X; {[Xtemp]}];
    Xsure = [Xsure; {[XX]} ];
    
end

%%
save FMRI_SURE_LOWRES I J K Xsure