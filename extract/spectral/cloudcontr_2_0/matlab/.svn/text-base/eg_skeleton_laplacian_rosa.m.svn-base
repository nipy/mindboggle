% function [] = eg_point_cloud_curve_skeleton(filename)
% extract curve skeleton from a point cloud or triangular mesh
% update: 2010-8-19
% update: 2010-7-12
% create: 2009-4-26
% by: JJCAO, deepfish @ DUT
%
%% setting
clear;clc;close all;
path('toolbox',path);
options.USING_POINT_RING = GS.USING_POINT_RING;
extension='.off';

%% Step 0: read file (point cloud & local feature size if possible), and
% normalize the modle.
filename = '../data/simplejoint_v4770';% which file we should run on

tic
P.filename = [filename extension];% point set
[P.pts,P.faces] = read_mesh(P.filename);
P.npts = size(P.pts,1);
if exist([filename '_fe.txt'],'file') % result of Tamal K Dey's NormFet
    P.radis = load([filename '_fe.txt']);
else
    P.radis = ones(P.npts,1);
end

P.pts = GS.normalize(P.pts);
[P.bbox, P.diameter] = GS.compute_bbox(P.pts);
disp(sprintf('read point set:'));
toc

%% Step 1: build local 1-ring
% build neighborhood, knn?
tic
P.k_knn = GS.compute_k_knn(P.npts);
if options.USING_POINT_RING
    P.rings = compute_point_point_ring(P.pts, P.k_knn, []);
else    
    P.frings = compute_vertex_face_ring(P.faces);
    P.rings = compute_vertex_ring(P.faces, P.frings);
end
disp(sprintf('compute local 1-ring:'));
toc

%% Step 1: Contract point cloud by Laplacian
tic
[P.cpts, t, initWL, WC, sl] = contraction_by_mesh_laplacian(P, options);
disp(sprintf('Contraction:'));
toc
%% step 2: Point to curve ¨C by cluster ROSA2.0
tic
P.sample_radius = P.diameter*0.02;
P = rosa_lineextract(P,P.sample_radius, 1);
disp(sprintf('to curve:'));
toc
%% show results
figure('Name','Original point cloud and its contraction');movegui('northeast');set(gcf,'color','white')
scatter3(P.pts(:,1),P.pts(:,2),P.pts(:,3),30,'.','MarkerEdgeColor', GS.PC_COLOR);  hold on;
scatter3(P.cpts(:,1),P.cpts(:,2),P.cpts(:,3),30,'.r'); axis off;axis equal;set(gcf,'Renderer','OpenGL');
camorbit(0,0,'camera'); axis vis3d; view(0,90);view3d rot;


figure('Name','Original point cloud and its skeleton'); movegui('center');set(gcf,'color','white');
scatter3(P.pts(:,1),P.pts(:,2),P.pts(:,3),20,'.','MarkerEdgeColor', GS.PC_COLOR);  hold on;
showoptions.sizep=400;showoptions.sizee=2;
plot_skeleton(P.spls, P.spls_adj, showoptions);
axis off;axis equal;set(gcf,'Renderer','OpenGL');view(0,90);view3d rot;
%% save results
default_filename = sprintf('%s_contract_t(%d)_nn(%d)_WL(%f)_WH(%f)_sl(%f)_skeleton.mat',...
    P.filename(1:end-4), t, P.k_knn, initWL, WC, sl);
save(default_filename,'P');