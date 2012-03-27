function labelsMatrix = analyzeAllLabels(allLabels,allDistances,subjectList)

distThreshold =3;

labelsMatrix = zeros(37,37);

for i = 1:40
   inds = find(subjectList == i); 
   tempLabelsMatrix = zeros(37,37); 
   
   for j = 1:length(inds)
       if (allDistances(inds(j)) < distThreshold)
           inMin = min(allLabels(inds(j),:)) +1;
           inMax = max(allLabels(inds(j),:)) + 1;
           
           tempLabelsMatrix(inMin,inMax) = tempLabelsMatrix(inMin,inMax) +1;
       end
   end   
   tempLabelsMatrix(tempLabelsMatrix > 1) =1;
   labelsMatrix = labelsMatrix + tempLabelsMatrix;
end

for i = 1:37
    for j = 1:37
        if (labelsMatrix(i,j)>35)
            disp(i);
            disp(j);
            disp(labelsMatrix(i,j));
            disp('****');
        end        
    end
end


