function [spls, corresp, spls_adj, joints, segments] = remove_irrelevant_extrama(spls, corresp, spls_adj, joints, segments,global_dist,threshold,SHOW_IRRELEVANT_EXTRAMA)


tmp = global_dist(global_dist ~= Inf);
skel_size = max(tmp) - min(tmp);
while true
    edges = zeros(0,3); % i, j, distance
    for i=1:size(joints,1) % irrelevant extrama are leaves on joints.
        joint = int32(joints(i,1));
        links = find( spls_adj(joint,:)==1 );
        links = links( links~=joint );
        for j =links
            if length(find( spls_adj(j,:)==1 ))==2 % is_leaf
                edges(end+1,:) = [double(joint), double(j), abs(global_dist(joint) - global_dist(j))/skel_size];
            end
        end
    end

    if isempty(edges), break, end;
    edges = sortrows(edges,3);
    if edges(1,3) > threshold, break, end; 
    % remove shortest edges
    i = int32(edges(1,1)); j = int32(edges(1,2));
    segments(j) = 0;
    corresp( corresp==j ) = i;
    spls(j,:)=NaN;
    spls_adj(j,:)=0;
    spls_adj(i,j)=0;
    
    links = find( spls_adj(i,:)==1 );
    if length(links) < 4 % not a joint agian        
        segments(i) = NaN;
        ind = find(joints(:,1)==i);
        joints(ind,:) = [];
    end    
end
% draw joints, draw roots; 
if SHOW_IRRELEVANT_EXTRAMA
    figure('Name','Remove irrelevant extrama','NumberTitle','off');set(gcf,'color','white');
    movegui('center');hold on;
    plot_skeleton(spls, spls_adj);
    scatter3(spls(joints(:,1),1),spls(joints(:,1),2),spls(joints(:,1),3),200,'b','filled');
%     scatter3(spls(root_id,1),spls(root_id,2),spls(root_id,3),30,'r','filled');
    axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);view3d rot;
end