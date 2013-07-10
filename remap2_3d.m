
function [map,actives] = remap2_3d(B,I,J)
%function to map the matrix of covariates to ijk space

[numrows,numsubjects ] = size(B);

map = cell(0);
actives = cell(0);
for s = 1:numsubjects
    b = B(:,s);
    M = zeros(numrows,4);
    for l = 1:numrows
        [i,j,k] = map2ijk(l,I,J);
        vec = [i j k b(l)];
        M(l,:) = vec;
        
    end
    map = [map ; {[M]}];
    actinds = find(M(:,4)~=0);
    acts = M(actinds,:);
    actives = [actives,{[acts]}];
    disp(s);
end
end

function [i,j,k] = map2ijk(n,I,J)

R = rem(n,I*J);
i = rem(R,I);
j = floor(R/I) + 1;
k = floor(n/I/J) + 1;

end
