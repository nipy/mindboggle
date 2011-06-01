function [output]=applyHierarchical(adj_matrix,options)
% The method apply the hierarchical clustering algorithm to an adjacency 
% matrix in order to separate the different cluster of the graph.
% The graph is assumed to be weighted.
%
% Input:
% - adj_matrix: adjacency matrix (weighted)
% - options: structured variable containing hte options for the hierarchical 
%            algorithm. It includes hier.distance, hier.linkage 
%            hier.dendrogram_Z_threshold, hier.grapn_number (optional)
%
% Output:
% - output: structured variable containing the node cluster label and a new
%           version of the adjacency matrix where the nodes are ordered so
%           that nodes belonging to the same cluster are close to each
%           others.


if ~isfield(options,'hier.distance')
    options.hier.distance = 'correlation';
end
if ~isfield(options,'hier.linkage_method')
    options.hier.linkage_method = 'average';
end
if ~isfield(options,'hier.dendrogram_Z_threshold')
    options.hier.dendrogram_Z_threshold = 0.9;
end
if ~isfield(options,'grapn_number')
    options.grapn_number = 000;
end

% 1) Apply hierarchical clustering algorithm
distance_vector=pdist(adj_matrix,options.hier.distance);
tree_struct = linkage(distance_vector,options.hier.linkage_method);
cluster_indices=cluster(tree_struct,'MaxClust',options.n_cluster);

sorted_indices = [];
for cont = 1:max(cluster_indices)
    new_elements = find(cluster_indices==cont);
    sorted_indices = [sorted_indices; new_elements];
end

        
% Sorting output matrices
% sorting Adjacency matrix
sorted_adj = adj_matrix;
indices = 1:1:length(sorted_indices)';
for cont =1:length(sorted_indices)
    if indices(cont)~=sorted_indices(cont)
        elem1 = find(indices==sorted_indices(cont));
        elem2 = cont;
        sorted_adj=move_elements(sorted_adj,elem1,elem2);
        indices(elem1) = indices(cont);
        indices(cont) = sorted_indices(cont);
    end
end
% Sorting cluster indices
sorted_cluster_indices = cluster_indices(sorted_indices);

% % Code for imaging
% figure
% subplot(121)
% imagesc(sorted_adj,[0,0.05])
% subplot(122)
% stem(sorted_cluster_indices)

output.sorted.adj_matrix = sorted_adj;
output.sorted.cluster_indices = sorted_cluster_indices;
output.unsorted.cluster_indices = cluster_indices;

% % Code for imaging
% figure
% set(gcf,'Units','inches','Position',[2 2 13 6])
% dendrogram(tree_struct,0,'colorthreshold',options.hier.dendrogram_Z_threshold);
% set(gca,...
%     'fontSize',12)
% xlabel('Cluster')
% ylabel('Tree')
% title('Graph dendrogram')
% if options.saveImages
%     saveas(gcf,[options.images_basename,'_graph',num2str(options.grapn_number),'_dendrogram'],'png')
% end