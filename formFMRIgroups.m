function G = formFMRIgroups(I,J,K,siz,overlap)
% function to form groups in n-space given i,j,k

if mod(siz,2)==0
    error('siz has to be odd');
end

offset = (siz-1)/2;


if ~overlap
irange = offset+1:2*offset+1:I-offset;
jrange = offset+1:2*offset+1:J-offset;
% krange = offset+1:2*offset+1:K-offset;
% ^comment  above line and keep below line for nXnX1
krange = 1:K;
else
irange = offset+1:offset:I-offset;
jrange = offset+1:offset:J-offset;
% krange = offset+1:offset:K-offset;
% ^comment  above line and keep below line for nXnX1
krange = 1:K;
end

% FOR EACH I, J, K IN THE RANGE, FORM GROUP AND MAP TO N SPACE

gcount = 0;
for i = irange
    disp(i);
    for j = jrange
        for k = krange
            gcount = 1+gcount;
            v = [i j k];
            center(gcount,:) = v;
            % center is the center of each group
        end
    end
end

M = size(center,1);
G = cell(0);
for ii = 1:M
    IJKS = makegrp(center(ii,:),offset);
    len = size(IJKS,1);
    g = zeros(1,len);
    for l=1:len
        g(l) = map2n(IJKS(l,:),I,J);
    end
    g = sort(g);
    G = [G; {[g]}];
    
end

% add the "rest" in another group
allinds = 1:I*J*K;
allsel = [];
M = length(G);
for grpind = 1:M;
    tempvec = G{grpind};
    allsel = [allsel tempvec];
end
allsel = unique(allsel);
rest = setdiff(allinds,allsel);
G = [G; {[rest]}];

end


function grp = makegrp(center,d)
i = center(1); j = center(2) ; k = center(3);
krange = k-d:k+d; 
krange = k;% have this line for nXnX1 group
jrange = j-d:j+d;
irange = i-d:i+d;  

grp = [];
for p = irange
    for q = jrange
        for r = krange
            vec = [p q r];
            grp = [grp; vec];
        end
    end
end


end

function n = map2n(center,I,J)
i = center(1); j = center(2); k = center(3);
n = i + I*(j-1) + I*J*(k-1);
end



