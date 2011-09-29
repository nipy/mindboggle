function [root_id, global_dist] = find_root_node(spls, spls_adj, joints, show_root)
if isempty(joints)
    root_id = 1;
else
    [c, root_id] = max(joints(:,2));
    root_id = joints(root_id,1);
end

local_dist = zeros( size(spls_adj)); % local distance matrix
for i=1:size(spls_adj,1)
    for j=1:size(spls_adj,2)
        if( spls_adj(i,j)==1 )
            rs = spls(i,:) - spls(j,:);
            local_dist(i,j) = sqrt(dot(rs,rs));
        end
    end
end

global_dist = graphshortestpath(sparse(local_dist), root_id);

% draw joints, draw roots; 
if show_root
    scatter3(spls(root_id,1),spls(root_id,2),spls(root_id,3),400,'r', 'filled');
end