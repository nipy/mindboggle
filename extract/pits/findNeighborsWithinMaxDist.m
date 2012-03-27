function [neighbors neighborCount] = findNeighborsWithinMaxDist(vertices,faces,maxDist,centroidIndex,pointOfInterestIndex,checkList)

inds1 = faces(:,1) == pointOfInterestIndex;
inds2 = faces(:,2) == pointOfInterestIndex;
inds3 = faces(:,3) == pointOfInterestIndex;

inds = [faces(inds1,:);faces(inds2,:);faces(inds3,:)];
inds = inds(:);
inds = unique(inds);


neighbors = -1;
neighborCount = 0;

for i = 1:length(inds)
    if (checkList(inds(i)) < 0)
        if (compEucDist(vertices(centroidIndex,:), vertices(inds(i),:)) < maxDist)
            neighborCount = neighborCount + 1;
            neighbors(neighborCount) = inds(i);
        end
    end
end
