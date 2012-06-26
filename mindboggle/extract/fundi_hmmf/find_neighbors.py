function [inds] = findNeighbors(faces,pointOfInterestIndex)

inds1 = faces(:,1) == pointOfInterestIndex;
inds2 = faces(:,2) == pointOfInterestIndex;
inds3 = faces(:,3) == pointOfInterestIndex;


inds = [faces(inds1,:);faces(inds2,:);faces(inds3,:)];
inds = inds(:);
inds = unique(inds);


% inds = zeros(20,1);
% indsCount = 0;
% 
% for i = 1:length(inds1)
%    tempInds = faces(inds1(i),:);
%    indsCount = indsCount + 3;
%    inds(indsCount - 2 : indsCount) = tempInds;   
% end
% for i = 1:length(inds2)
%    tempInds = faces(inds2(i),:);
%    indsCount = indsCount + 3;
%    inds(indsCount - 2 : indsCount) = tempInds;   
% end
% for i = 1:length(inds3)
%    tempInds = faces(inds3(i),:);
%    indsCount = indsCount + 3;
%    inds(indsCount - 2 : indsCount) = tempInds;   
% end
% 
% inds = unique(inds);

inds(inds == pointOfInterestIndex) = 0;
inds = inds(inds > 0);

