%%% script to make the data cell arrays

clear; clc; close all;
nums = [314 1407 1577 2177 2184 2199 2200 2202 2227,...
    2232 2347 3283 3303]; %remove subject 1225

load mcw_metadata
% load fmri_classification_metadata
I = []; J = []; K = [];icount = 0; jcount = 0; kcount = 0;
for jj = 1:length(nums)
    
    disp(jj);
    str = strcat('load s',num2str(nums(jj)),'.aseg_or_ribbon.blur.mat');
    eval(str);
%     str = strcat('load s',num2str(nums(jj)),'.aseg_sort.mat');
%     eval(str);
    str = strcat('ijk = XYZ_tlrc.s',num2str(nums(jj)),';');
    eval(str);
%     str = strcat('ijk = xyz_tlrc.s',num2str(nums(jj)),';');
%     eval(str);
    ijk = round(ijk);
    
    imin(jj) = min(ijk(:,1));
    jmin(jj) = min(ijk(:,2));
    kmin(jj) = min(ijk(:,3));
    
end
imin = min(imin); jmin=min(jmin);kmin = min(kmin);
for jj = 1:length(nums)
    
    disp(jj);
        str = strcat('load s',num2str(nums(jj)),'.aseg_or_ribbon.blur.mat');
        eval(str);
%     str = strcat('load s',num2str(nums(jj)),'.aseg_sort.mat');
%     eval(str);
        str = strcat('ijk = XYZ_tlrc.s',num2str(nums(jj)),';');
        eval(str);
%     str = strcat('ijk = xyz_tlrc.s',num2str(nums(jj)),';');
%     eval(str);
    ijk = round(ijk);
    ijk(:,1) = ijk(:,1)-imin;
    ijk(:,2) = ijk(:,2)-jmin;
    ijk(:,3) = ijk(:,3)-kmin; % make the minimum = 0
    %check if there are zeros
    if ~isempty(find(ijk(:,1)==0))
        fprintf('I has a zero \n')
    end
    if ~isempty(find(ijk(:,2)==0))
        fprintf('J has a zero \n')
    end
    if ~isempty(find(ijk(:,3)==0))
        fprintf('K has a zero \n')
    end
    
    ijk = ijk + 1; % make the minimum = 1
    
    I = [I; max(ijk(:,1))];
    J = [J; max(ijk(:,2))];
    K = [K; max(ijk(:,3))];  % see how big these numbers get
    
end
I = max(I);J = max(J);K = max(K); %total dimension = IJK

%% MAPPING FROM ijk ---> n = i + I(j-1) + IJ(k-1)

% MAKE NEW DATA MATRICES 
X = cell(0);
Xsure = cell(0);
for jj = 1:length(nums)
    disp(jj);
    str = strcat('load s',num2str(nums(jj)),'.aseg_or_ribbon.blur.mat');
    eval(str);
%     str = strcat('load s',num2str(nums(jj)),'.aseg_sort.mat');
%     eval(str);
    str = strcat('ijk = XYZ_tlrc.s',num2str(nums(jj)),';');
    eval(str);
%     str = strcat('ijk = xyz_tlrc.s',num2str(nums(jj)),';');
%     eval(str);
    ijk = round(ijk);
    ijk(:,1) = ijk(:,1)-imin;
    ijk(:,2) = ijk(:,2)-jmin;
    ijk(:,3) = ijk(:,3)-kmin; % make the minimum = 0
    ijk = 1+ijk;
    
    numvox = size(funcData,2);
    Xtemp = sparse(1200,I*J*K);
    str = strcat('b = ReduxObj_Blur.s',num2str(nums(jj)),'.voxels;');
    eval(str);
    ambvox = find(b==0);
     
    for ii = 1:numvox
        % remove ambiguous voxels 
        if ismember(ii,ambvox)
            continue
        end
        i = ijk(ii,1); j = ijk(ii,2); k = ijk(ii,3);
        n = i + I*(j-1) + I*J*(k-1);
        Xtemp(:,n) = funcData(:,ii);
        
    end
    Xtemp = Xtemp(1:900,:);
    %remove ambiguous words
    ambwrds = find(Ambiguous(1:900)==0);
    Xtemp = Xtemp(ambwrds,:);
    X = [X; {[Xtemp]}];
%     Xsure = [Xsure; {[XX]} ];
    
end
%%
overlap = 1; siz = 21;
G = formFMRIgroups(I,J,K,siz,overlap);
%%
save FMRI_BIG_OR_BLURRED I J K imin jmin kmin G
% the Xs are too big so save separately
for i = 1:length(nums)
    Xtemp = X{i};
%     str = strcat('X',num2str(nums(i)),' = Xtemp;');
%     eval(str);
    str = strcat('save X',num2str(i),' Xtemp');
    eval(str);
    disp(i)
end