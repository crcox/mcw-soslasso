%% FIND RANGE OF I J K
clear; clc; close all
nums = [314 1225 1407 1577 2177 2184 2199 2200 2202 2227 2232 2347 3283 3303];

I = []; J = []; K = [];icount = 0; jcount = 0; kcount = 0;
for i = 1:length(nums)
    filename = sprintf('s%d.aseg.mat',nums(i));
    load(filename,'ijk');
    
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
    
    I = [I; max(ijk(:,1))];
    J = [J; max(ijk(:,2))];
    K = [K; max(ijk(:,3))];
end
I = max(I);J = max(J);K = max(K);

%% MAPPING FROM ijk ---> n = i + I(j-1) + IJ(k-1)

% MAKE NEW DATA MATRICES 
for jj = 1:length(nums)
    disp(jj);
    filename = sprintf('s%d.aseg.mat',nums(jj));
    load(filename,'ijk');
    for ii = 1:length(ijk)
        i = ijk(ii,1); j = ijk(ii,2); k = ijk(ii,3);
        n = i + I*(j-1) + I*J*(k-1);
    end
end

%%
save FMRI_SURE_LOWRES I J K Xsure
