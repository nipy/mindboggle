function Sulci = extractSulci(faces,Depth,depthTh)

Sulci = zeros(size(Depth));

seeds = find(Depth>depthTh);

seedSize = size(seeds,1);
counter = 0;

% loop until all points included in sulci
while (seedSize > 0)
    counter = counter + 1;
    TEMP = zeros(size(Depth));
    
    %choose random seed point (selection does not affect result)
    rseed = round(rand*(seedSize-1)) + 1;
    
    TEMP(seeds(rseed,1),1) = 2;
    newSize = size(find(TEMP>1));
    
    % grow region until no more connected points available
    while(newSize > 0)        
        indsList = find(TEMP == 2);
        
        TEMP(TEMP == 2) = 1;
        
        for i = 1:size(indsList,1)
            neighs = findNeighbors(faces,indsList(i,1));
            neighs = neighs(Depth(neighs) > depthTh);           
            neighs = neighs(TEMP(neighs) == 0);
            TEMP(neighs) = 2;                            
        end
        newSize = size(find(TEMP>1));
    end
    Sulci(TEMP > 0) = counter;
        
    Depth(Sulci > 0) = depthTh - .5;       
    seeds = find(Depth>depthTh);

    seedSize = size(seeds,1);
    
    %%%%%%%%%%
    % display current sulcus size
    disp(size(find(Sulci == counter),1));
    pause(.0001);
    %%%%%%%%%%
end
