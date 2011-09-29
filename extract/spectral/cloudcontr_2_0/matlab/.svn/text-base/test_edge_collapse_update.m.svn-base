path('toolbox',path);
%% read pcd and show it
clear;clc;close all;    
filename = '../data/2D/n.mat';
load(filename);
P.pts = M.verts;
clear('M');
npts = size(P.pts, 1);
[P.bbox, P.diameter] = GS.compute_bbox(P.pts);

%% call farthest_sampling_by_sphere
P.sample_radius = P.diameter*0.02;
[P.spls,P.corresp] = farthest_sampling_by_sphere(P.pts, P.sample_radius);

%% call connect_by_inherit_neigh
k=5;
kdtree = kdtree_build(P.pts);
P.neigh = zeros(npts, k);
for i = 1:npts
    P.neigh(i,:)  = kdtree_k_nearest_neighbors(kdtree, P.pts(i,:), k)';
end
    
P.spls_adj = connect_by_inherit_neigh(P.pts, P.spls, P.corresp, P.neigh);

%% call edge_collapse
options.collapse_order = 0;
[P.spls, P.skel_adj,P.corresp] = edge_collapse_update(P.pts, P.spls, P.corresp, P.neigh, P.spls_adj, options);
figure; movegui('center');set(gcf,'color','white');
plot3( P.spls(:,1), P.spls(:,2), P.spls(:,3), '.r', 'markersize', 1);
hold on;axis off; axis equal;set(gcf,'Renderer','OpenGL');
GS.plot_connectivity(P.spls, P.skel_adj, 1);
view3d zoom;

kdtree_delete( kdtree );