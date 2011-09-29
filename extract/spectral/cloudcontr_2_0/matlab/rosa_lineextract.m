function P = rosa_lineextract(P, RADIUS, bUpdateConnectivity)
%
% we RECOVER CONNECTIVITY after each edge collapse, which is slow, but the
% result sampling is more uniform.
% bUpdateConnectivity: 0 for RECOVER CONNECTIVITY after each edge collapse, which is slow, but the
% result sampling is more uniform.
% bUpdateConnectivity: 1 for update CONNECTIVITY after each edge collapse, which is faster, and the
% result sampling is almost as same as that of the original method (bUpdateConnectivity=0).
%
% @author: Andrea Tagliasacchi
% reformed by jjcao
% @reform-data:     2009-9-3   using 1-ring to construct the connectivity
% graph
% @reform-data:     2010-8-20  decompose the original function into several
% small functions.
SHOW_RESULTS = true;

options.collapse_order = 1;
if (bUpdateConnectivity)
    [P.spls,P.corresp] = farthest_sampling_by_sphere(P.cpts, P.sample_radius);
    P.spls_adj = connect_by_inherit_neigh(P.cpts, P.spls, P.corresp, P.rings);    
    [P.spls, P.spls_adj,P.corresp] = edge_collapse_update(P.cpts, P.spls, P.corresp, P.rings, P.spls_adj, options);
else
    [P.spls,P.corresp] = farthest_sampling_by_sphere(P.cpts, P.sample_radius);
    P.spls_adj = connect_by_inherit_neigh(P.cpts, P.spls, P.corresp, P.rings);
    [P.spls, P.spls_adj,P.corresp] = edge_collapse(P.cpts, P.spls, P.corresp, P.rings, P.spls_adj, options);
end

%%
if (SHOW_RESULTS)
    figure; movegui('northeast');set(gcf,'color','white');hold on;
    plot3( P.spls(:,1), P.spls(:,2), P.spls(:,3), '.r', 'markersize', 5);
    axis off; axis equal;set(gcf,'Renderer','OpenGL');
    plot_connectivity(P.spls, P.spls_adj);
    view3d zoom;
end

