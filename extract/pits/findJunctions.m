function labelsMatrix = findJunctions(allLabels,allDistances,subjectList)

distThreshold =3;

labelsMatrix = zeros(35,35,35);

for i = 1:40
    inds = find(subjectList == i);
    tempLabelsMatrix = zeros(35,35,35);
    
    for j = 1:length(inds)
        if (allDistances(inds(j),2) < distThreshold)
            currInds = sort(allLabels(inds(j),1:3));
            
            if (currInds == [3 24 28])
                disp(i);
                disp(j);
                disp(length(inds))
                disp('^^');
            end
            
            tempLabelsMatrix(currInds(1),currInds(2),currInds(3)) = 1;
        end
    end
    labelsMatrix = labelsMatrix + tempLabelsMatrix;
end

for i = 1:35
    for j = 1:35
        for k = 1:35
            if (labelsMatrix(i,j,k)>20)
                disp(i);
                disp(j);
                disp(k);
                disp(labelsMatrix(i,j,k));
                disp('****');
            end
        end
    end
end


