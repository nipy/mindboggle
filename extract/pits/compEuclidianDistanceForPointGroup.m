function D = compEuclidianDistanceForPointGroup(pitInds,pits)

m = length(pitInds);

D = zeros(m,m);

for i = 1:m
    for j = m:-1:1
        if (j > i && D(i,j) == 0)
            D(i,j) = compEucDist(pits(i,:),pits(j,:));
            D(j,i) = D(i,j);
            
%            if (compEucDist(pits(i),pits(j)) < maxED)
                                
%                 options.end_points = pitInds(j);
%                 [geoD,S] = perform_dijkstra(W,pitInds(i), options);
% 
%                  for k = 1:m
%                    if (geoD(pitInds(k)) < 100000)
%                        D(i,k) = geoD(pitInds(k));
%                        D(k,i) = geoD(pitInds(k));
%                    end                    
%                 end
%            end
        end
    end
end




