function checkList = findLocalMaxima(vertices,faces,curv,checkList,centroidIndex,pointOfInterestIndex,maxDist,curvature)

if (curv(pointOfInterestIndex) > curvature)
    checkList(centroidIndex) = 0;
    return;
else
    
    if (pointOfInterestIndex ~= centroidIndex)
        checkList(pointOfInterestIndex) = 0;
    end
    %[neighbors neighborCount] = findNeighborsWithinMaxDist(vertices,faces,maxDist,centroidIndex,pointOfInterestIndex,checkList);
    neighbors = findNeighborsWithinMaxDist(faces,pointOfInterestIndex);
    
    for i = 1:neighborCount
        if (checkList(neighbors(i)) ~= 0)
            checkList = findLocalMaxima(vertices,faces,curv,checkList,centroidIndex,neighbors(i),maxDist,curvature);
            if (checkList(centroidIndex == 0))
                return;
            end
        end
    end
end


