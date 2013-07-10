%%% script to eliminate "commonly zero" data points from the matrices

clear; clc; close all;
load FMRI_BIG_OR_BLURRED
nums = [314 1407 1577 2177 2184 2199 2200 2202 2227,...
    2232 2347 3283 3303];
X = cell(0);
for i = 1:length(nums)
    str = strcat('load X',num2str(i),' ;');
    eval(str);
    disp(i)
    X = [X; {Xtemp}];
end
%%
Xsure = X;
clear X;
numpersons = length(Xsure);
nonzero_locations = cell(0);
for ii = 1:numpersons
    disp(ii);
    X = Xsure{ii};
    s = sum(X~=0);
    inds = find(s~=0);
    nonzero_locations = [nonzero_locations; {[inds]}];
end
% nonzero_locations has the locations of data for each subject

%% the actual non zero locations are given by the union of these locations
true_nonzero = [];
for ii = 1:numpersons
    disp(ii);
    true_nonzero = [true_nonzero nonzero_locations{ii}];
    true_nonzero = unique(true_nonzero);
end

%% now, make new data to that only retains the true nonzeros
Xsure_full = Xsure;
clear Xtrue;
Xsure = cell(0);
for ii = 1:numpersons
    disp(ii)
    Xtemp = Xsure_full{ii};
    Xtemp = Xtemp(:,true_nonzero);
    Xsure = [Xsure; {[Xtemp]}];
end

%% we now need to make new groups, based on the old ones. 

% look at each active component, and find out what group(s) actually
% contain it. 
G_full = G;
M = length(G);
clear G;
for ii = 1:M
    disp(ii);
    g = G_full{ii};
    g = intersect(g,true_nonzero);
    G{ii} = g;
end

%% remove empty groups
nonblank_index = [];
for ii = 1:M
    g = G{ii};
    if ~isempty(g)
        nonblank_index = [nonblank_index ii];
    end
end
G = G(nonblank_index);
%%
% remove repeated groups
G = uniquecell(G);
G = G';
%% now, for each of the new groups, map it to the small new space we have
G_oldspace = G;
M = length(G);
clear G;
for ii = 1:M
    disp(ii);
    g = G_oldspace{ii};
    vec = [];
    for i = 1:length(g)
        temp = g(i);
        ind = find(true_nonzero == temp);
        vec = [vec ind];
    end
    G{ii} = vec;
end
G = G';
X = Xsure;
clear Xsure;
%% save shit
save FMRI_SMALL_OR_BLURRED I J K G G_oldspace true_nonzero imin jmin kmin

for i = 1:length(nums)
    Xtemp = X{i};
%     str = strcat('X',num2str(nums(i)),' = Xtemp;');
%     eval(str);
    str = strcat('save Xs',num2str(i),' Xtemp');
    eval(str);
    disp(i)
end