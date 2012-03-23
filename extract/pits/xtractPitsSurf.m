%function [vertices faces curv checkList] = xtractPitsSurf(surf_filename,curv_filename,subjectId,side)
function [vertices faces curv pitInds pits] = xtractPitsSurf(faces,vertices,curv)

%finds local maxima of curvature. For other measures, such as depth,
%change variable 'curv'

maxDist = 2.0;
curvThreshold = .01;

% Add freesurfer matlab codes to path
%addpath('/Applications/freesurfer/matlab')

% Use freesurfer to read a surface
%[vertices,faces] = freesurfer_read_surf(surf_filename);
%curv = read_curv(curv_filename);

[m n] = size(vertices);

checkList = -1*ones(m,1);

tic

%for i = 1:m
for i = 1:m     
    
    if (curv(i) < curvThreshold)
        checkList(i) = 0;
    end
        
    if (checkList(i) == -1)
        %disp(i);
        
        tempCheckList = -1*ones(m,1);
        tempCheckList(i) = 1;
        
        tempCheckList = findLocalMaxima(vertices,faces,curv,tempCheckList,i,i,maxDist,curv(i));
        
        if (tempCheckList(i) == 1)
            checkList(tempCheckList == 0) = 0;
            checkList(i) = 1;
            
            disp('************');
            disp(curv(checkList == 1));
            
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

