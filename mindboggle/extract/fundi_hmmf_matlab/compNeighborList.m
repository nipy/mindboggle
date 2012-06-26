function nList = compNeighborList(faces,inds)

m = size(inds,1);

disp(m);

nList = zeros(m,15);

for i = 1:m
    neighs = findNeighbors(faces,inds(i));
    nList(i,1:size(neighs,1)) = neighs';
end

