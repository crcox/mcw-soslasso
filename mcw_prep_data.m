%% FIND RANGE OF I J K
clear; clc; close all
nums = [314 1225 1407 1577 2177 2184 2199 2200 2202 2227 2232 2347 3283 3303];
load('mcw_metadata.mat','IJK','XYZ_tlrc','ReduxObj','ReduxObj_Blur');
BLUR = false;
if BLUR
    ReduxObj = ReduxObj_Blur;
end
GroupSize = 21;
%% Compute minimum number of points needed to represent all data in TLRC
temp = min(cell2mat(struct2cell(XYZ_tlrc)));
imin = temp(1); 
jmin = temp(2);
kmin = temp(3);

temp = max(cell2mat(struct2cell(XYZ_tlrc)));
I = temp(1) - imin + 1; % zero base, then one base 
J = temp(2) - jmin + 1;
K = temp(3) - kmin + 1;

%% Inflate slightly, as needed, to accomodate group size (and overlap).
I = I + (GroupSize-mod(I-GroupSize,floor(GroupSize/2)));
J = J + (GroupSize-mod(J-GroupSize,floor(GroupSize/2)));
K = K + (GroupSize-mod(K-GroupSize,floor(GroupSize/2)));

%% MAPPING FROM ijk ---> n = i + I(j-1) + IJ(k-1)
IJK_TLRC = bsxfun(@minus,cell2mat(struct2cell(XYZ_tlrc)),[imin,jmin,kmin])+1;
IND_TLRC = uint32(unique(sub2ind([I,J,K],IJK_TLRC(:,1),IJK_TLRC(:,2),IJK_TLRC(:,3))));
% Breakdown:
% 1. Take coordinates for all subjects as one large matrix, subtract the
% minimums, and add 1. Transformation to 1 based IJK indexes in TLRC space.
% 2. Transform those IJK indexes to literal column indexes. Retain only the
% unique ones, and store as an array of unsigned 32 bit integers.
clear IJK_TLRC;

NZ_TLRC = false(1,I*J*K);
NZ_TLRC(IND_TLRC) = true;

% MAKE NEW DATA MATRICES 
[X,Xsure] = deal(cell(nsubj,1));

switch Task
case 'ani-all'
    Words  = [true(900,1);false(300,1)];               % 1200x1  logical
    y      = TrueAnimals(~Ambiguous & Words);          %  809x1  logical
    
case 'art-all'
    Words  = [true(900,1);false(300,1)] & ~Ambigious;  % 1200x1  logical
    y      = TrueArtifacts(~Ambiguous & Words);        %  809x1  logical

case 'ani-art'
    Words  = [true(240,1);false(960,1)]  & ~Ambigious; % 1200x1  logical
    y      = TrueAnimals(~Ambiguous & Words);          %  235x1  logical
end

for jj = 1:length(nums)
    S = sprintf('s%d',nums(jj));
    filename = sprintf('%d.aseg.mat',nums(jj));
    load(filename,'funcData');
    xyz_tlrc = XYZ_tlrc.(S);
    ijk_tlrc = bsxfun(@minus,xyz_tlrc,[imin,jmin,kmin]) + 1;
    ijk_tlrc = ijk_tlrc(ReduxObj.(S).voxels,:); % remove outlier voxels.
    ind_tlrc = sub2ind([I J K],ijk_tlrc(:,1),ijk_tlrc(:,2),ijk_tlrc(:,3));
    
    funcData_tlrc = zeros(sum(Words),I*J*K);
    funcData_tlrc(:,ind_tlrc) = funcData;
    funcData_tlrc = sparse(funcData_tlrc(:,NZ_TLRC));
    
    filename = sprintf('%d.aseg.tlrc.mat',nums(jj));
    save(filename,'funcData_tlrc','ijk_tlrc','ind_tlrc');
    clear funcData funcData_tlrc ijk xyz_tlrc;
end

%% Define groups
overlap = true; 
offset = (GroupSize-1)/2;

if ~overlap
    irange = 1:(2*offset+1):(I-(2*offset+1));
    jrange = 1:(2*offset+1):(J-(2*offset+1));
    krange = 1:(2*offset+1):(K-(2*offset+1));
else
    irange = 1:offset:(I-GroupSize);
    jrange = 1:offset:(J-GroupSize);
    krange = 1:offset:(K-GroupSize);
end

[IG,JG,KG] = ndgrid(irange,jrange,krange);
IJK_corners = [IG(:), JG(:), KG(:)]; % Grid of group corners.

[IG,JG,KG]= ndgrid(0:(siz-1),0:(siz-1),0:(siz-1));
GroupExtent = [IG(:), JG(:), KG(:)]; % Grid of distances from corner.

M = length(IJK_corners);
G = cell(M,1);

for i = 1:M
    G_ijk = bsxfun(@plus, IJK_corners(i,:), GroupExtent);
    G_ind = sub2ind([I,J,K], G_ijk(:,1),G_ijk(:,2),G_ijk(:,3))';
    G_ind_nz = G_ind(NZ(G_ind));
    G{i} = NZ_cumsum32(G_ind_nz);
end
G = G(~cellfun(@isempty,G));
G = uniquecell(G);
G = G';

%% Define the replicated, non-overlapping space.
[RepIndex, groups, group_arr] = define_rep_space(G);

filename = sprintf('GROUPS.size%d.offset%d.mat',GroupSize,offset);
save(filename,'RepIndex', 'groups', 'group_arr', 'G');


