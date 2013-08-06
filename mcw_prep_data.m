function [X,GroupInfo] = mcw_prep_data(metadata,params)
    %% Parameters
    GroupSize  = params.GroupSize;
    nums       = metadata.subjects;
    ReduxObj   = metadata.ReduxObj;
    Words      = metadata.Words;
    XYZ_tlrc   = metadata.XYZ_tlrc;
    numpersons = length(nums);
    
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
    GroupInfo.ijk_tlrc_range = [I J K];

    %% MAPPING FROM ijk ---> n = i + I(j-1) + IJ(k-1)
    IJK_TLRC = bsxfun(@minus,cell2mat(struct2cell(XYZ_tlrc)),[imin,jmin,kmin])+1;
    IND_TLRC = uint32(unique(sub2ind([I,J,K],IJK_TLRC(:,1),IJK_TLRC(:,2),IJK_TLRC(:,3))));
    clear IJK_TLRC;
    % Breakdown:
    % 1. Take coordinates for all subjects as one large matrix, subtract the
    % minimums, and add 1. Transformation to 1 based IJK indexes in TLRC space.
    % 2. Transform those IJK indexes to literal column indexes. Retain only the
    % unique ones, and store as an array of unsigned 32 bit integers.

    NZ_TLRC = false(1,I*J*K);
    NZ_TLRC(IND_TLRC) = true;

    % MAKE NEW DATA MATRICES 
    X = cell(numpersons,1);
    for jj = 1:length(nums)
        disp(jj)
        S = sprintf('s%d',nums(jj));
        filename = sprintf('s%d.aseg.mat',nums(jj));
        load(filename,'funcData');
        xyz_tlrc = XYZ_tlrc.(S);
        ijk_tlrc = bsxfun(@minus,xyz_tlrc,[imin,jmin,kmin]) + 1;
        ijk_tlrc = ijk_tlrc(ReduxObj.(S).voxels,:); % remove outlier voxels.
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

    [IG,JG,KG]= ndgrid(0:(GroupSize-1),0:(GroupSize-1),0:(GroupSize-1));
    GroupExtent = [IG(:), JG(:), KG(:)]; % Grid of distances from corner.

    M = length(IJK_corners);
    GroupInfo.G = cell(M,1);
    NZ_cumsum32 = uint32(cumsum(NZ_TLRC));
    for i = 1:M
        G_ijk = bsxfun(@plus, IJK_corners(i,:), GroupExtent);
        G_ind = sub2ind([I,J,K], G_ijk(:,1),G_ijk(:,2),G_ijk(:,3))';
        G_ind_nz = G_ind(NZ_TLRC(G_ind));
        if ~isempty(G_ind_nz)
            GroupInfo.G{i} = NZ_cumsum32(G_ind_nz);
        end
    end
    GroupInfo.G = GroupInfo.G(~cellfun(@isempty,GroupInfo.G));
    GroupInfo.G = uniquecell(GroupInfo.G);

    %% Define the replicated, non-overlapping space.
    [GroupInfo.RepIndex, GroupInfo.groups, GroupInfo.group_arr] = define_rep_space(GroupInfo.G);
    
    filename = sprintf('GROUPS.size%d.offset%d.mat',GroupSize,offset);
    save(filename,'-struct','GroupInfo');
end
