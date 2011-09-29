function [M, graph] = build_hierarchy(M, global_dist)
% nb: cannot contain any cycles
skelid = find( ~isnan(M.skelver(:,1)) )';% real skeleton nodes.
global prev;
prev = zeros(1,size(M.skelver,1));%previous skeleton node id
while true % do a Breadth-First-Search from root
    for i=1:size(M.skelver,1)
        if prev(i) ~=0, continue, end;
        if isnan( M.skelver(i,1) )
            prev(i) = nan;
            continue;
        end
        if i == M.root_id
            prev(i) = M.root_id;
            continue;
        end
        links = find( M.skel_adj(i,:)==1 );    
        if length(links)==2 %leaf node
            prev(i) = setdiff(links, i);        
            trace_skel(prev(i), i, M,global_dist);
        end
    end
    tmp = find(prev>0);
    if length(tmp) == length(skelid), break, end;
end

%% make root node the first node.
M.prev = zeros(1, length(skelid));
for i=1:length(skelid)
    tmp = prev(skelid(i));
    M.prev(i) = find(skelid==tmp);  
    M.corresp( M.corresp==skelid(i) ) = i;
end
M.root_id = find(skelid==M.root_id);

skelid=1:length(M.prev);
% M.skelid = [M.root_id, find(skelid~=M.root_id)];
M.prev = [M.prev(M.root_id), M.prev(find(skelid~=M.root_id))];
tmp1 = (M.prev==M.root_id);
tmp2 = (M.prev<M.root_id);
M.prev(tmp1) = 1;
M.prev(tmp2) = M.prev(tmp2)+1;
M.prev(1) = 0;
tmp1 = (M.corresp==M.root_id);
tmp2 = (M.corresp<M.root_id);
M.corresp(tmp1) = 1;
M.corresp(tmp2) = M.corresp(tmp2)+1;

M.skelver(isnan(M.skelver(:,1)),:)=[];
a=M.skelver;
a=[a(M.root_id,:);a(skelid~=M.root_id,:)];
M.skelver = a;

%%
skel_adj = eye(length(M.prev));
for i = 1:length(M.prev)
    if M.prev(i)>0
        skel_adj(i, M.prev(i)) = 1;
        skel_adj(M.prev(i), i) = 1;
    end
end
M.skel_adj = skel_adj;

%% build a skeleton graph (build skeleton edges)
graph = zeros(0,2);
for i = 1:length(M.prev)
    if M.prev(i)>0
        graph(end+1,:) = [i,M.prev(i)];
    end
end
%%
figure('Name','Build hierarchy','NumberTitle','off');set(gcf,'color','white');hold on; movegui('southeast');
plot_skeleton(M);
axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);view3d rot;

%%%%%%%%%%%%%%%%%%
function [] = trace_skel(id, nextid, M, global_dist)
global prev;
if id == M.root_id, return, end

links = find( M.skel_adj(id,:)==1 );    
if length(links)>3 %joint node
    tmp = setdiff(links, [id, nextid]);
    [C,I] = min(global_dist(tmp));
    prev(id) = tmp(I(1));
    trace_skel(prev(id), id, M,global_dist);
else %branche node
    prev(id) = setdiff(links, [id, nextid]);
    trace_skel(prev(id), id, M,global_dist);
end