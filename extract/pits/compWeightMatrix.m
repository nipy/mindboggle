function W = compWeightMatrix(fv)

[m n] = size(fv.vertices);

[mfaces nfaces] = size(fv.faces);

W = sparse(m,m);

for i = 1:mfaces   
    inds = fv.faces(i,:);
    
    W(inds(1),inds(2)) = compEucDist(fv.vertices(inds(1),:),fv.vertices(inds(2),:));
    W(inds(2),inds(1)) = W(inds(1),inds(2));
    
    W(inds(1),inds(3)) = compEucDist(fv.vertices(inds(1),:),fv.vertices(inds(3),:));
    W(inds(3),inds(1)) = W(inds(1),inds(3));
    
    W(inds(2),inds(3)) = compEucDist(fv.vertices(inds(2),:),fv.vertices(inds(3),:));
    W(inds(3),inds(2)) = W(inds(2),inds(3));        
    
    if (mod(i,10000) == 0)
       disp([num2str(100*i/mfaces) '% done']);
    end
    
end
