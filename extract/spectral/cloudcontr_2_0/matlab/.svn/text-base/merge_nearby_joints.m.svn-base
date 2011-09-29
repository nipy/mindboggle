function [spls, corresp, spls_adj, joints, root_id] = merge_nearby_joints(spls, corresp, spls_adj, joints, root_id)
while true
    for i=joints(:,1)'    
        edges = find( spls_adj(i,:)==1 );
        edges = edges(edges~=i);
        tmp = ismember(edges, joints(:,1)');
        if sum(tmp)~=1, continue, end;

        j = edges(tmp);
        edges = find( spls_adj(j,:)==1 );
        edges = edges(edges~=j);
        tmp = ismember(edges, joints(:,1)');
        if sum(tmp)~=1, continue, end;

      % update the location
        spls( i,: ) = mean( spls( [i,j],: ) );
        spls( j,: ) = NaN;
        % update the correspondents
       corresp(corresp==j ) = i;

        % update the A matrix
        for k=1:size(spls_adj,1)
            if spls_adj( j,k ) == 1, 
                spls_adj( i,k )=1; 
                spls_adj( k,i)=1; 
            end
        end
        % remove the row
        spls_adj( j,: ) = 0;
        spls_adj( :,j ) = 0;

        segments(j) = 0;
        j = find( joints(:,1)==j );
        if root_id == j
            root_id = find( joints(:,1)==i );            
        end
        joints(j,:) = [];
%         break;
    end
    if i == joints(end,1), break, end;
end

%%
figure('Name','Merge nearby joints','NumberTitle','off');set(gcf,'color','white');view3d rot;
movegui('north');
plot_skeleton(spls, spls_adj);
axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);view3d rot;