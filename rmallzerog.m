function G = rmallzerog(G,true_nonzero)
% RMALLZEROG Remove any indexes that point to empty voxels across all
% subjects. Also, ensure indexes are represented as integers rather than
% doubles.

M = length(G);
for i = 1:M
    G{i} = uint32(intersect(G{i},true_nonzero));
end