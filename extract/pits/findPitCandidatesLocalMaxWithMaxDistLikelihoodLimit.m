%function [vertices faces curv checkList] = xtractPitsSurf(surf_filename,curv_filename,subjectId,side)
function [pitInds pits] = findPitCandidatesLocalMaxWithMaxDistLikelihoodLimit(faces,vertices,depthVector,L,maxDist)

[m n] = size(vertices);

checkList = -1*ones(m,1);

tic
for i = 1:m
    
    if (L(i) < .5)
        checkList(i) = 0;
    end
    
    if (checkList(i) == -1)
        
        connectedNeighbors = findNeighbors(faces,i);
        
        if (max(depthVector(connectedNeighbors)) > depthVector(i))
            checkList(i) = 0;
        end
        
    end
    
    
    if (checkList(i) == -1)
        
        
        
        %disp(i);
        
        tempCheckList = -1*ones(m,1);
        tempCheckList(i) = 1;
        
        %tempCheckList = findLocalMaximaWithinRegion(vertices,faces,depthVector,tempCheckList,i,i,maxDist,depthVector(i));
        
        if (tempCheckList(i) == 1)
            checkList(tempCheckList == 0) = 0;
            checkList(i) = 1;
            
            disp('************');
            disp(depthVector(checkList == 1));
            
            disp('*');
            disp(find(checkList == 1));
            disp('*');
            
            disp(i);
            done = sum((checkList > -1));
            disp([num2str(100*done/m) '% done']);
            disp(done);
            disp(sum(checkList == 1));
            
            disp('************');
            
        else
            checkList(i) = 0;
        end
    end
end

pitInds = find(checkList == 1);

pits = zeros(length(pitInds),3);

for i = 1:length(pitInds)
    pits(i,:) = vertices(pitInds(i),:);
end

%a = outputPits(pits,subjectId,side);

toc


return;

