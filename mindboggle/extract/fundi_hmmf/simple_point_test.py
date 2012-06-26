function sp = isSimplePointMesh(fv,ind,values)

thr = 0.5;

neighs = findNeighbors(fv.faces,ind);

neighvals = values(neighs);

nOutside = sum(neighvals <= thr);
nInside = sum(neighvals > thr);

if (nOutside == 0 || nInside == 0)
    sp = 0;
    return;
end
if (nOutside == 1)
    sp = 1;
    return;
end
if (nInside == 1)
    sp = 1;
    return;
end

neighsInside = neighs(neighvals > thr);
nSets = size(neighsInside,1);

sets = zeros(nSets,20);
setLabels = [1:nSets]';

values(ind) = -1;

for i = 1:nSets
    currNeighs = findNeighbors(fv.faces,neighsInside(i));
    currNeighs = currNeighs(values(currNeighs)>thr);
    sets(i,1:size(currNeighs,1)) = currNeighs';
    sets(i,size(currNeighs,1) + 1) = neighsInside(i);
end

change = 1;
while (change > 0)
    change = 0;
    for i = 1:nSets
        for j = 1:nSets
           if (j>i && setLabels(i) ~= setLabels(j))
                iInds = unique(sets(i,:));
                iInds = iInds(iInds > 0);
                iSize = size(iInds,2);
                
                jInds = unique(sets(j,:));
                jInds = jInds(jInds > 0);
                jSize = size(jInds,2);
                
                totInds = unique([iInds jInds]);
                totSize = size(totInds,2);
                
                if (totSize < iSize + jSize)
                    
                    minLabel = min(setLabels(j),setLabels(i));
                    
                    setLabels(i) = minLabel;
                    setLabels(j) = minLabel;
                    
                    change = 1;                   
                end
           end            
        end
    end
end

numberOfSeparateSets = size(unique(setLabels),1);
if (numberOfSeparateSets < 2)
    sp = 1;
else
    sp = 0;
end
   

