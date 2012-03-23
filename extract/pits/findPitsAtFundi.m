function sulcusPits = findPitsAtFundi(allLabels,allDistances,subjectList)

distThreshold = 3;

% combine [2,10,23,26] as "2" for cingulate (remove 10, 23, and 26)
% combine [3,27] as "3" for middlefrontal (remove 27)
% combine [18,19,20] as "18" for inferiorfrontal
%
% with these three new aggregate labels (2,3,18), here is a reasonable collection of major lateral folds (sulci):
%
% 1) superior frontal: [28,3]
% 2) inferior frontal: [3,18]
% 3) sylvian fissure: [18,30]
% 4) superior temporal: [30,15]
% 5) inferior temporal: [15,9]
% 6) precentral: [[28,24],[3,24],[18,24]]
% 7) postcentral: [24,22]
% 8) intraparietal: [29,8]
%
% and just a few medial sulci:
%
% 9) cingulate: [[2,28],[2,17],[25,17]]
% 10) parietooccipital fissure: [5,25]
% 11) calcarine fissure: [5,13]

sulcusPits = zeros(40,11);

allLabels = modifyLabels(allLabels);

for i = 1:size(allLabels,1)
    for j = 1:3
        if (allDistances(i,j) < distThreshold)
            sulc = findCorrectSulcus(allLabels(i,1),allLabels(i,j+1));
            if (sulc > 0)
                sulcusPits(subjectList(i),sulc) = sulcusPits(subjectList(i),sulc) + 1;                
            end
        end
    end
end
