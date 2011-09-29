function [spls, corresp, spls_adj, joints, segments] = remove_inner_nodes(spls, corresp, spls_adj, joints, segments,global_dist, len_thr, angle_thr)
SHOW_INNEREDGECOLLAPSE_PROGRESS = true;

tmp = global_dist(global_dist ~= Inf);
skel_size = max(tmp) - min(tmp);

if SHOW_INNEREDGECOLLAPSE_PROGRESS
    figure;set(gcf,'color','white');hold on; movegui('northeast');
    axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);view3d rot;
end
while true
    edges = zeros(0,3); % i, j, distance
    branches = find(isnan(segments))';
    for i=branches
        links = find( spls_adj(i,:)==1 );
        if length(links) == 2, continue, end;
        links = setdiff(links, i);
        angle = myangle( spls([links(1), i, links(2)],:) );
        if angle<angle_thr, continue, end;
        
        for j = links
%             if segments(j) < 0, continue, end;
            
            dist = abs(global_dist(i) - global_dist(j))/skel_size;            
            if dist>len_thr, continue,end;
            
            edges(end+1,:) = [double(i), double(j), dist];
        end 
    end
    
    if isempty(edges), break, end;
    edges = sortrows(edges, 3);
    % remove shortest edges
    i = int32(edges(1,2)); j = int32(edges(1,1));
    segments(j) = 0;
    spls(j,:)=NaN;
    corresp( corresp==j ) = i;    
    % update the A matrix
    for k=1:size(spls_adj,1)
        if spls_adj(j,k) == 1, 
            spls_adj(i,k)=1; 
            spls_adj(k,i)=1; 
        end
    end
    spls_adj(j,:)=0;
    spls_adj(:,j)=0;  
    
    if SHOW_INNEREDGECOLLAPSE_PROGRESS
        hpts = scatter3(spls(:,1),spls(:,2),spls(:,3), '.r'); 
%         hpts = mypoint3( spls, '.r'); 
        heds = [];
        for i=1:size(spls_adj,1)
            for j=1:size(spls_adj,2)
                if( spls_adj(i,j)==1 )
                    idx = [i;j];
                    heds(end+1) = line( spls(idx,1),spls(idx,2),spls(idx,3));
%                     heds(end+1) = myedge3(spls(i,:), spls(j,:)); %#ok<AGROW>
                end
            end
        end
        axis off; axis equal;set(gcf,'Renderer','OpenGL');
        drawnow;
        delete(hpts);
        delete(heds);
    end
end
if SHOW_INNEREDGECOLLAPSE_PROGRESS
    close;
end
%%
figure('Name','Remove inner nodes','NumberTitle','off');set(gcf,'color','white');
hold on; movegui('northeast');
plot_skeleton(spls, spls_adj);
axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);view3d rot;

function beta = myangle(verts)
u = verts(1,:) - verts(2,:);
v = verts(3,:) - verts(2,:);
du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) );