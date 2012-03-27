function checkList = findLocalMaximaWithinRegion(vertices,faces,depthVector,checkList,centroidIndex,pointOfInterestIndex,maxDist,depth)

if (depthVector(pointOfInterestIndex) > depth)
    checkList(centroidIndex) = 0;
    return;
else
    
    if (pointOfInterestIndex ~= centroidIndex)
        checkList(pointOfInterestIndex) = 0;
    end
    [neighbors neighborCount] = findNeighborsWithinMaxDist(vertices,faces,maxDist,centroidIndex,pointOfInterestIndex,checkList);
    %neighbors = findNeighborsWithinMaxDist(faces,pointOfInterestIndex);
    
    for i = 1:neighborCount
        if (checkList(neighbors(i)) ~= 0)
            checkList = findLocalMaximaWithinRegion(vertices,faces,depthVector,checkList,centroidIndex,neighbors(i),maxDist,depth);
            if (checkList(centroidIndex) == 0)
                
                return;
            end
        end
    end
end
