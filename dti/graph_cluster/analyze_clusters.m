function [cluster_out]=analyze_clusters(adj_matrix,cluster_indices)
% The function analyze a graph using the brain connectivity toolbox
% Input:
% - adj_matrix: the graph adjacency matrix
% - cluster_indices: the node cluster label (provide a vector of ones if no
%   cluster has been done.
% Output:
% - cluster_out: a structured variable containing for each cluster the
%   analysis performed.

nC = max(cluster_indices);
for c=1:nC
    disp([' - analysis of cluster ',num2str(c)])
    sub_graph = extract_graph(adj_matrix,find(cluster_indices==c));
    graph_analysis = characterize_graph(sub_graph);
    field_list =fieldnames(graph_analysis);
    for cont = 1:length(field_list)
        fieldID = field_list{cont};
        cluster_out(c).(fieldID)=graph_analysis.(fieldID);
    end
    module_analysis_results = module_characterization(adj_matrix,cluster_indices==c);
    field_list =fieldnames(module_analysis_results);
    for cont = 1:length(field_list)
        fieldID = field_list{cont};
        cluster_out(c).(fieldID)=module_analysis_results.(fieldID);
    end
end
% graph_analysis = characterize_graph(adj_matrix);
% field_list =fieldnames(graph_analysis);
% for cont = 1:length(field_list)
%     fieldID = field_list{cont};
%     graph_out(c).(fieldID)=graph_analysis.(fieldID);
% end





function [sub_graph]=extract_graph(adj_matrix,vertex_set)
sub_graph = adj_matrix(vertex_set,vertex_set);

function [output]=characterize_graph(adj_matrix)
% The method compute a set of measures on the provided graph to
% characterize it.

if size(adj_matrix,1)<=1
    % The graph is empty
    
    %                                          Basic consepts and notations
    % ---------------------------------------------------------------------
    v.basic.weighted_degree = 0;                                   % Degree (vertex measure)
    v.basic.shortest_path_length = 0;                       % Shortest path length (vertex measure)
    v.basic.weighted_triangle_mean = 0;               % Weighted triangle mean (vertex measure)
    
    %                                               Measures of integration
    % ---------------------------------------------------------------------
    v.integration.characteristic_path_length = 0;
    v.integration.eccentricity = 0;
    v.integration.radius = 0;
    v.integration.diameter = 0;
                                                                           % Characteristic path length
                                                                           % Eccentricity (vertex measure)
                                                                           % radius
                                                                           % diameter
    v.integration.global_efficiency = 0;                    % Global efficiency
    
    %                                               Measures of segregation
    % ---------------------------------------------------------------------
    v.segregation.clustering_coefficient= 0; 
    % the clustering coefficient is defined for each node. The graph 
    % clustering coefficient is defined as the average clustering 
    % coefficient among the nodes.
    v.segregation.transitivity = 0;                            % Transitivity coefficient
    v.segregation.local_efficiency = 0;                   % Local efficiency (vertex measure)
    v.segregation.modularity_Q = 0;                        % modularity coefficient
    
    %                                               Measures of segregation
    % ---------------------------------------------------------------------
    % closeness centrality (do implement)
    v.centrality.btw_centrality = 0;                          % Betweenness centrality (Vertex measure)
    v.centrality.within_module_degree_zscore = 0;     % Within module z-score (vertex measure)
    v.centrality.participation_coeff = 0;               % Participation coefficient (vertex measure)
    
    %                                                        Network motifs
    % ---------------------------------------------------------------------
    % Defined only for directed graphs
    
    %                                                Measures of resilience
    % ---------------------------------------------------------------------
    v.resilience.cum_deg_dist = 0;
    v.resilience.deg_dist_grid = 0;                                        % Degree distribution (graph measure)
    v.resilience.aver_neighbor_deg = 0;                                    % Average neighbor degree (vertex measure)
    v.resilience.assortativity_coef = 0;                                   % Assortativity coefficient
    
    %                                                       Small Worldness
    % ---------------------------------------------------------------------
    v.small_worldness = 0;
else
    % The graph is not empty
    
    %                                          Basic consepts and notations
    % ---------------------------------------------------------------------
    v.basic.weighted_degree = sum(adj_matrix,2);                           % Degree (vertex measure)
    v.basic.shortest_path_length = distance_wei(adj_matrix);               % Shortest path length (vertex measure)
    v.basic.weighted_triangle_mean = weighted_triangles(adj_matrix);       % Weighted triangle mean (vertex measure)
    
    %                                               Measures of integration
    % ---------------------------------------------------------------------
    [v.integration.characteristic_path_length,v.integration.eccentricity,...
        v.integration.radius, v.integration.diameter] = charpath(v.basic.shortest_path_length);
                                                                           % Characteristic path length
                                                                           % Eccentricity (vertex measure)
                                                                           % radius
                                                                           % diameter
    v.integration.global_efficiency = efficiency_wei_und(adj_matrix);      % Global efficiency
    
    %                                               Measures of segregation
    % ---------------------------------------------------------------------
    v.segregation.clustering_coefficient=mean(clustering_coef_wu(adj_matrix)); 
    % the clustering coefficient is defined for each node. The graph 
    % clustering coefficient is defined as the average clustering 
    % coefficient among the nodes.
    v.segregation.transitivity = transitivity_wu(adj_matrix);              % Transitivity coefficient
    v.segregation.local_efficiency = efficiency_wei_und(adj_matrix,1);     % Local efficiency (vertex measure)
    [Ci,v.segregation.modularity_Q] = modularity_und(adj_matrix);          % modularity coefficient
    
    %                                               Measures of segregation
    % ---------------------------------------------------------------------
    % closeness centrality (do implement)
    v.centrality.btw_centrality = betweenness_wei(adj_matrix);             % Betweenness centrality (Vertex measure)
    v.centrality.within_module_degree_zscore = module_degree_zscore(adj_matrix,Ci);     % Within module z-score (vertex measure)
    v.centrality.participation_coeff = participation_coef(adj_matrix,Ci);  % Participation coefficient (vertex measure)
    
    %                                                        Network motifs
    % ---------------------------------------------------------------------
    % Defined only for directed graphs
    
    %                                                Measures of resilience
    % ---------------------------------------------------------------------
    [v.resilience.cum_deg_dist,v.resilience.deg_dist_grid] = degree_dist_wei(adj_matrix);            % Degree distribution (graph measure)
    v.resilience.aver_neighbor_deg = average_neighbor_degree(adj_matrix);               % Average neighbor degree (vertex measure)
    v.resilience.assortativity_coef = assortativity_wei(adj_matrix);                    % Assortativity coefficient
    
    %                                                       Small Worldness
    % ---------------------------------------------------------------------
    Cw = v.segregation.clustering_coefficient;
    Lw = v.integration.characteristic_path_length;
    Cw_rand = zeros(100,1);
    Lw_rand = zeros(100,1);
    disp('Generating random graphs')
    for cont=1:100
        rand_Adj = randmio_und_connected2(adj_matrix,10);
        Cw_rand(cont) = mean(clustering_coef_wu(rand_Adj));
        D_matrix = distance_wei(rand_Adj);
        [Lw_temp,e,radius_temp,diameter_temp] = charpath(D_matrix);
        Lw_rand(cont) = Lw_temp;
    end
    v.small_worldness = (Cw/mean(Cw_rand))/(Lw/mean(Lw_rand));
%    v.small_worldness = 0;
end
output=v;

function [output]=module_characterization(adj_matrix,module)

% Module analysis
output.module_degree = module_degree_zscore(adj_matrix,double(module)); % measure of centrality
output.module_participation = participation_coef(adj_matrix,double(module));