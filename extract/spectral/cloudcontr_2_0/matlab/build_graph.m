function [spls, corresp, spls_adj, graph] = build_graph(spls,corresp,spls_adj)
graph = zeros(0,2);
for i=1:size(spls_adj,1)
    for j=i+1:size(spls_adj,2)
        if( spls_adj(i,j)==1 )
            graph(end+1,:) = [i, j];
        end
    end
end
tmp = find(isnan(spls(:,1)))';
tmp = tmp(length(tmp):-1:1);
for i=tmp
    graph(graph>i) = graph(graph>i) - 1;
    corresp(corresp>i) = corresp(corresp>i) -1;
end
spls(tmp,:) = [];

spls_adj = zeros(size(spls,1));
for a=1:size(graph,1)
    i = graph(a,1); j = graph(a,2);
    spls_adj(i,j) = 1;
    spls_adj(j,i) = 1;
end
%%
figure('Name','Build hierarchy','NumberTitle','off');set(gcf,'color','white');hold on; movegui('southeast');
plot_skeleton(spls,spls_adj);
axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);view3d rot;
