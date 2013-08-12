function [X,GroupInfo] = mcw_prep_data(metadata,params,varargin)
    %% Parameters
    GroupSize  = params.GroupSize;
    overlap    = params.overlap; 
    offset     = params.offset(GroupSize); % anonymous function. 
    nums       = metadata.subjects;
    ReduxObj   = metadata.ReduxObj;
    Words      = metadata.Words;
    XYZ_tlrc   = metadata.XYZ_tlrc;
    numpersons = length(nums);
    OnlyGroupInfo = false;
    RemoveOutliersLater = false;
    if nargin > 2
        if any(strcmp('OnlyGroupInfo',varargin))
            OnlyGroupInfo = true;
        end
        if any(strcmp('RemoveOutliersLater',varargin))
            RemoveOutliersLater = true;
        end
    end
    
    %% Compute minimum number of points needed to represent all data in TLRC
    if ~RemoveOutliersLater
        % Outliers are removed from XYZ_tlrc immediately.
        for s = 1:numpersons
            S = sprintf('s%d',nums(s));
            XYZ_tlrc.(S) = XYZ_tlrc.(S)(ReduxObj.(S).voxels,:);
        end
    end
    
    %% Concatenate all coordinates
    XYZ_tlrc_mat = cell2mat(struct2cell(XYZ_tlrc)); 
    
    %% Find the range of each dimension that contains all coordinates.
    minXYZ = min(XYZ_tlrc_mat);    RemoveOutliersLater = false;

    maxXYZ = max(XYZ_tlrc_mat);
    maxIJK = maxXYZ - minXYZ + 1;
    I = maxIJK(1);
    J = maxIJK(2);
    K = maxIJK(3);
    clear maxXYZ maxIJK;
    
    %% Inflate slightly, as needed, to accomodate group size (and overlap).
    if params.AccommodateGroups == true
        d = (GroupSize - 2*offset);
        I = I + mod(I,offset) + d;
        J = J + mod(J,offset) + d;
        K = K + mod(K,offset) + d;
        clear d;
    end
    GroupInfo.ijk_tlrc_range = [I J K];

    %% MAPPING FROM ijk ---> n = i + I(j-1) + IJ(k-1)
    IJK_TLRC = bsxfun(@minus,XYZ_tlrc_mat,minXYZ)+1;
    IND_TLRC = uint32(unique(sub2ind([I,J,K],IJK_TLRC(:,1),IJK_TLRC(:,2),IJK_TLRC(:,3)),'stable'));
    IND_TLRC = sort(IND_TLRC); % In Nikhil's it seems to be sorted...
    clear IJK_TLRC;
    % Breakdown:
    % 1. Take coordinates for all subjects as one large matrix, subtract the
    % minimums, and add 1. Transformation to 1 based IJK indexes in TLRC space.
    % 2. Transform those IJK indexes to literal column indexes. Retain only the
    % unique ones, and store as an array of unsigned 32 bit integers.

    NZ_TLRC = false(1,I*J*K);
    NZ_TLRC(IND_TLRC) = true;
    GroupInfo.NZ_TLRC = NZ_TLRC;
    

    % MAKE NEW DATA MATRICES 
    X = cell(numpersons,1);
    if ~OnlyGroupInfo
        
        fprintf('\n');
        for jj = 1:length(nums)
            fprintf('% 6d',nums(jj))
            S = sprintf('s%d',nums(jj));
            filename = sprintf('s%d.aseg.mat',nums(jj));
            load(filename,'funcData');
            xyz_tlrc = XYZ_tlrc.(S);
            ijk_tlrc = bsxfun(@minus,xyz_tlrc,minXYZ) + 1;
            if RemoveOutliersLater
                ijk_tlrc = ijk_tlrc(ReduxObj.(S).voxels,:);
            end
            ind_tlrc = sub2ind([I J K],ijk_tlrc(:,1),ijk_tlrc(:,2),ijk_tlrc(:,3));

            funcData_tlrc = zeros(sum(Words),I*J*K);
            funcData_tlrc(:,ind_tlrc) = funcData(Words,ReduxObj.(S).voxels);
            funcData_tlrc = sparse(funcData_tlrc(:,NZ_TLRC));

            if params.Save == true
                filename = sprintf('%d.aseg.tlrc.mat',nums(jj));
                save(filename,'funcData_tlrc','ijk_tlrc','ind_tlrc');
            end

            X{jj} = funcData_tlrc;
            clear funcData funcData_tlrc ijk xyz_tlrc;
        end
    end
    fprintf('\n')

    %% Define Range 
    if ~overlap
        irange = 1:GroupSize:(I-GroupSize);
        jrange = 1:GroupSize:(J-GroupSize);
        krange = 1:GroupSize:(K-GroupSize);
    else
        irange = 1:offset:(I-GroupSize);
        jrange = 1:offset:(J-GroupSize);
        krange = 1:offset:(K-GroupSize);
    end

    %% Generate Grid of Corners
    [KG,JG,IG] = ndgrid(krange,jrange,irange);
    IJK_corners = [IG(:), JG(:), KG(:)]; % Grid of group corners.

    %% Generate Grid-patch to make groups
    [IG,JG,KG]= ndgrid(0:(GroupSize-1),0:(GroupSize-1),0:(GroupSize-1));
    GroupExtent = [IG(:), JG(:), KG(:)]; % Grid of distances from corner.

    %% Make G (full)
    M = length(IJK_corners);
    G = cell(M,1);
    for i = 1:(M)
        G_ijk = bsxfun(@plus, IJK_corners(i,:), GroupExtent);
        temp = sub2ind([I,J,K], G_ijk(:,1),G_ijk(:,2),G_ijk(:,3))';
        G{i} = uint32(temp(NZ_TLRC(temp)));
    end
    temp = setdiff(IND_TLRC,cell2mat(G'));
    G{M+1} = uint32(temp(NZ_TLRC(temp)));
    clear IND_TLRC;
    %% Remove Empty Groups
    G = G(~cellfun('isempty',G));

    %% Remove Repeated Groups 
    G = uniquecell(G);
    G = G';

    %% Map to small space
    GroupInfo.G_oldspace = G;
    M = length(G);
    G = cell(M,1);

    nth_NZ_TLRC = repmat(uint32(0),1,length(NZ_TLRC));
    nth_NZ_TLRC(NZ_TLRC) = 1:sum(NZ_TLRC);
    for ii = 1:M
        G{ii} = nth_NZ_TLRC(GroupInfo.G_oldspace{ii});
    end
    GroupInfo.G = G;
    clear G;

    %% Define the replicated, non-overlapping space.
    [GroupInfo.RepIndex, GroupInfo.groups, GroupInfo.group_arr] = define_rep_space(GroupInfo.G);
    
    if params.Save == true
        filename = sprintf('GROUPS.size%d.offset%d.mat',GroupSize,offset);
        save(filename,'-struct','GroupInfo');
    end
end
