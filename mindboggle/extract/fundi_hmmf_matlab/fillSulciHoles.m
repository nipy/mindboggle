function Sulci = fillSulciHoles(faces,Sulci)

%%%%%%%%%%%
disp(['Number of sulci before pruning: ' num2str(max(Sulci))]);

counter = 0;

while (counter < max(Sulci))
    counter = counter + 1;
    ns = sum(Sulci == counter);
    if (ns < 50)
        Sulci(Sulci == counter) = 0;
        Sulci(Sulci > counter) = Sulci(Sulci > counter) - 1;
        counter = counter - 1;
    end    
end

disp(['Number of sulci after pruning: ' num2str(max(Sulci))]);
%%%%%%%%%%%

Holes = zeros(size(Sulci));

seeds = find(Sulci == 0);

seedSize = size(seeds,1);
counter = 0;

while (seedSize > 0)
    counter = counter + 1;
    TEMP = zeros(size(Sulci));
    
    rseed = round(rand*(seedSize-1)) + 1;
    
    TEMP(seeds(rseed,1),1) = 2;
    newSize = size(find(TEMP>1));
    
    % grow region until no more connected points available
    while(newSize > 0)
        indsList = find(TEMP == 2);
        
        TEMP(TEMP == 2) = 1;
        
        for i = 1:size(indsList,1)
            neighs = findNeighbors(faces,indsList(i,1));
            neighs = neighs(Sulci(neighs) == 0);
            neighs = neighs(TEMP(neighs) == 0);
            TEMP(neighs) = 2;
        end
        newSize = size(find(TEMP>1));
    end
    Holes(TEMP > 0) = counter;
    
    Sulci(Holes > 0) = .5;
    seeds = find(Sulci == 0);
    
    seedSize = size(seeds,1);
    
    %%%%%%%%%%
    % display current sulcus size
    disp(size(find(Holes == counter),1));
    pause(.0001);
    %%%%%%%%%%
end
Sulci(Sulci < 1) = 0;

for i = 1:max(Holes)
    found = 0;
    currInds = find(Holes == i);
    
    if (size(currInds,1) < 10000)        
        j = 0;
        while (found == 0)
            j = j + 1;
            neighs = findNeighbors(faces,currInds(j,1));
            nIdx = max(Sulci(neighs));
            if (nIdx > 0)
                Sulci(currInds) = nIdx;
                found = 1;
                disp(['found ' num2str(i)]);
            end
        end
    end
end




