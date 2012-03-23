function minDist = findMinDistOfExtractedPits(D,Q,pitInds)

minDist = Inf;

[m n] = size(D);

for i = 1:m
    for j = 1:n
        if (Q(pitInds(j))> 0.5 && D(i,j) > 0 && Q(pitInds(i)) > 0.5 && D(i,j) < minDist)
            minDist = D(i,j);
        end
    end
end
