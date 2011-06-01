function [output]=applyKmeans(adj_matrix,options)
% The method apply the k-means algorithm to an adjacency matrix in order to
% separate the different cluster of the graph.
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

% Defining options that are not defined in the input data
if ~isfield(options,'kmeans.distance')
    options.kmeans.distance = 'correlation';
end
if ~isfield(options,'kmeans.start')
    options.kmeans.start = 'sample';
end
if ~isfield(options,'kmeans.replicates')
    options.kmeans.replicates = 5;
end
if ~isfield(options,'kmeans.emptyAction')
    options.kmeans.emptyAction = 'drop';
end
if ~isfield(options,'kmeans.onlinePhase')
    options.kmeans.onlinePhase = 'on';
end
if ~isfield(options,'kmeans.estimatorOptions')
    options.kmeans.estimatorOptions = statset('Display','final','MaxIter',200);
end

[cluster_indices, centroids, within_centroid_distance, centroid_distances] = kmeans(adj_matrix,options.n_cluster,...
    'distance',options.kmeans.distance,...
    'start',options.kmeans.start,...
    'replicates',options.kmeans.replicates,...
    'emptyaction',options.kmeans.emptyAction,...
    'onlinephase',options.kmeans.onlinePhase,...
    'options',options.kmeans.estimatorOptions);

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
sorted_centroids = centroids;
sorted_centroid_distances = zeros(size(centroid_distances));
for cont = 1:max(cluster_indices);
    sorted_centroid_distances(:,cont)=centroid_distances(sorted_indices,cont);
end

figure
subplot(221)
imagesc(sorted_adj,[0,0.05])
subplot(223)
stem(sorted_cluster_indices)
subplot(222)
imagesc(sorted_centroid_distances)

output.sorted.adj_matrix = sorted_adj;
output.sorted.cluster_indices = sorted_cluster_indices;
output.sorted.centroids = sorted_centroids;
output.sorted.centroid_distances = sorted_centroid_distances;

output.unsorted.cluster_indices = cluster_indices;
output.unsorted.centroids = centroids;
output.unsorted.centroid_distances = centroid_distances;

output.within_centroid_distance = within_centroid_distance;
