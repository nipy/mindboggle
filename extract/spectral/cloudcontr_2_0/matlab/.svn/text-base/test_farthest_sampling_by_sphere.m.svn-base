path('toolbox',path);
%% read pcd and show it
clear;clc;close all;    
filename = '../data/2D/n.mat';
load(filename);
P.pts = M.verts;
clear('M');
nverts = size(P.pts, 1);
[P.bbox, P.diameter] = GS.compute_bbox(P.pts);

%% call farthest_sampling_by_sphere
P.sample_radius = P.diameter*0.02;
[P.spls,P.corresp] = farthest_sampling_by_sphere(P.pts, P.sample_radius);

if ~all(P.corresp~=0)
    warning('some points have no correspondence');
end

figure; movegui('northeast');set(gcf,'color','white');hold on;
plot3( P.pts(:,1), P.pts(:,2), P.pts(:,3), '.r', 'markersize', 1);
axis off; axis equal;set(gcf,'Renderer','OpenGL');
plot3( P.spls(:,1), P.spls(:,2), P.spls(:,3), '*g');

