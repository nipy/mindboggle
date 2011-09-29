clear;clc;close all; 
path('toolbox',path);

filename = '../data/simplejoint_v4770';
extension='.off';

M.filename = [filename extension];
[M.verts,M.faces] = read_mesh(M.filename);
M.nverts = size(M.verts,1);

atria = nn_prepare(M.verts); % OpenTSTool,find kNN
M.k_knn = GS.compute_k_knn(M.nverts);

M.rings = compute_point_point_ring(M.verts, M.k_knn);

%% show results
figure;movegui('northeast');set(gcf,'color','white'); hold on;axis off;axis equal;set(gcf,'Renderer','OpenGL');
scatter3(M.verts(:,1),M.verts(:,2),M.verts(:,3),30,'.','MarkerEdgeColor', GS.PC_COLOR); 
camorbit(0,0,'camera'); axis vis3d; view3d rot;
%view(0,90);