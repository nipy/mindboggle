
function [vertices faces pitInds] = findPitCandidates(faces,vertices,depth)

%finds local maxima of depth.

depthThreshold = -Inf;

[m n] = size(vertices);

checkList = -1*ones(m,1);

tic

for i = 1:m
    
    if (depth(i) < depthThreshold)
        checkList(i) = 0;
    end
    
    if (checkList(i) == -1)
        %disp(i);
        
        checkList(i) = 0;
        neighbors = findNeighbors(faces,i);
        
        if (max(depth(neighbors)) < depth(i))
            checkList(i) = 1;
            checkList(neighbors) = 0;
            
            disp('************');
            disp(depth(checkList == 1));
            
            disp('*');
%             disp(find(checkList == 1));
%             disp('*');
            disp(i);
            done = sum((checkList > -1));
            disp([num2str(100*done/m) '% done']);
            disp(done);
            
            disp('# pit candidates found:')
            disp(sum(checkList == 1));
            disp('************');
            
        end
        
        
        
    end
end

pitInds = find(checkList == 1);

% pits = zeros(length(pitInds),3);
% 
% for i = 1:length(pitInds)
%     pits(i,:) = vertices(pitInds(i),:);
% end

toc


return;

