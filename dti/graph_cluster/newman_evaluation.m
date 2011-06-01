% STUDY OF NEWMAN CLUSTER METHOD
% Code to evaluate the performances of the Newman clustering method.
% Questions:
% - method robustness
% - evaluation of cluster quality
% - define the cluster number

%clear all
%close all
%clc
%options.n_cluster = 3;

% Options
data_file = '/data/export/home/denis/DTI_connection_graphs/graph_data/correlationGraphs.mat';
result_file = '/data/export/home/denis/DTI_connection_graphs/graph_cluster/newman_evaluation_study';
options.saveImages = true;
options.images_basename = ['newman'];

load(data_file,'CNT_graphs') % Also MDD_graphs
graph_set = CNT_graphs;
clear CNT_graphs

load(result_file) % Results have already been prepared


% Q1: method robustness
% Apply the same method to all graphs and see if the method perform
% similarly among the different graphs
[nG,nV]=size(Ci_table);

% Q1a: analyzing all graphs
cluster_table = old_Ci_table;
figure
set(gcf,'Units','inches','Position',[2 2 13 6])
subplot(121)
imagesc(cluster_table),colorbar
set(gca,...
    'fontS',12,...
    'XTick',[0:10:nV],...
    'YTick',[0:1:nG],...
    'FontSize',12)
title('Graph set cluster indices')
xlabel('Vertex')
ylabel('Graph')

% Q1b: automatic matching
matched_cluster_table = Ci_table;
subplot(122)
imagesc(matched_cluster_table),colorbar
set(gca,...
    'fontS',12,...
    'XTick',[0:10:nV],...
    'YTick',[0:1:nG],...
    'FontSize',12)
title('Matched clusters indices')
xlabel('Vertex')
ylabel('Graph')
if options.saveImages
    saveas(gcf,[options.images_basename,'_cluster_indices'],'png')
end

% % % Q1b: manual graph matching
% % % After a manual match of the cluster
% % matched_cluster_table = cluster_table;
% % for g = 1:nG
% %     figure(1)
% %     clf
% %     imagesc(matched_cluster_table,[0 options.n_cluster]),colorbar
% %     hold on
% %     plot([0.5,nV+0.5,nV+0.5,0.5,0.5],[g-0.5,g-0.5,g+0.5,g+0.5,g-0.5],'r--')
% %     title('Graph set cluster indices')
% %     xlabel('Vertex')
% %     ylabel('Graph')    
% %     disp(['Graph #',num2str(g)])
% %     for c=1:max(cluster_table(g,:))
% %         ind = find(cluster_table(g,:)==c);
% %         new_cluster_index = input(['New index for cluster ',num2str(c),': ']);
% %         matched_cluster_table(g,ind)=new_cluster_index;
% %     end
% % end
% % figure(1)
% % clf
% % imagesc(matched_cluster_table,[0 options.n_cluster]),colorbar
% % title('Graph set cluster indices')
% % xlabel('Vertex')
% % ylabel('Graph')
% % pause
% % close

disp(' ')
disp(' ------------------------------')
disp(' | Cluster method consistency |')
disp(' ------------------------------')

% Q1c: analysis of the matched graphs
% Overlap
nC = max(max(matched_cluster_table));
cluster_overlap = zeros(nC,1);
disp(' ')
disp(' Clustering overlap among the different graphs:')
for c=1:nC
    intersection_elements = 1:1:nV;
    union_elements = [];
    for g=1:nG
        cluster_elements = find(matched_cluster_table(g,:)==c);
        intersection_elements = intersect(intersection_elements,cluster_elements);
        union_elements = union(union_elements, cluster_elements);
    end
    cluster_overlap(c) = length(intersection_elements)/length(union_elements);
    disp(['  - cluster ',num2str(c),': ',num2str(cluster_overlap(c))])
end
figure
set(gcf,'Units','inches','Position',[2 2 6 6])
stem(cluster_overlap)
xlim([0.5 nC+0.5])
ylim([0 1])
set(gca,...
    'fontS',12,...
    'XTick',[0:1:nC],...
    'YTick',[0:0.2:1],...
    'FontSize',12)
xlabel('Cluster')
ylabel('overlap')
title('Cluster overlap among graphs')
if options.saveImages
    saveas(gcf,[options.images_basename,'_cluster_overlap'],'png')
end


% Element frequency
disp(' ')
disp(' Computing cluster element frequency')
element_probability = zeros(nC,nV);
for v=1:nV
    for c=1:nC
        element_probability(c,v) = length(find(matched_cluster_table(:,v)==c))/nG;
    end
end
% Element frequency analysis
element_max_probability = max(element_probability);
prob_x_vett = 0:0.05:1;
prob_y_vett = zeros(size(prob_x_vett));
for p=1:length(prob_y_vett)
    prob_y_vett(p)=length(find(element_max_probability>prob_x_vett(p)))./nV;
end
figure
set(gcf,'Units','inches','Position',[2 2 13 6])
subplot(121)
imagesc(element_probability,[0 1]),colorbar
set(gca,...
    'fontS',12,...
    'XTick',[0:10:nV],...
    'YTick',[0:1:nC],...
    'FontSize',12)
xlabel('Vertex')
ylabel('Cluster')
title('Assignment probability')
subplot(122)
set(gca,...
    'fontS',12,...
    'XTick',[0:0.2:1],...
    'YTick',[0:0.2:1],...
    'FontSize',12)
plot(prob_x_vett,prob_y_vett,'b-')
axis([0 1.1 0 1.1])
xlabel('maximum probability')
ylabel('amount of vertex')
title('Maximum probability distribution')
if options.saveImages
    saveas(gcf,[options.images_basename,'_cluster_consistency'],'png')
end



% Q2: cluster quality
disp(' ')
disp(' -------------------')
disp(' | Cluster quality |')
disp(' -------------------')
% Cluster density
disp(' ')
disp(' Computing global density')
global_graph_density = squeeze(sum(sum(graph_set,2),1)./nV);
% Internal cluster density
disp(' ')
disp(' Computing internal density')
local_density_table = zeros(nG,nC);
normalized_local_density_table = zeros(nG,nC);
for c=1:nC
    for g=1:nG
        ind = find(matched_cluster_table(g,:)==c);
        cluster_adj_matrix = graph_set(ind,ind,g);
        local_density_table(g,c)=sum(sum(cluster_adj_matrix))./length(ind);
        normalized_local_density_table(g,c)=local_density_table(g,c)./global_graph_density(g);
    end
end
figure
set(gcf,'Units','inches','Position',[2 2 6 6])
imagesc(normalized_local_density_table,[0 2]),colorbar
set(gca,...
    'fontS',12,...
    'XTick',[0:1:nC],...
    'YTick',[0:1:nG],...
    'FontSize',12)
xlabel('Cluster')
ylabel('Graph')
title('Normalized internal cluster density')
if options.saveImages
    saveas(gcf,[options.images_basename,'_internal_density'],'png')
end

% External cluster density
disp(' ')
disp(' Computing external density')
external_density_table = zeros(nG,nC);
for c=1:nC
    for g=1:nG
        ind = find(matched_cluster_table(g,:)==c);
        graph_adj_matrix = graph_set(:,:,g);
        partial_sum = 0;
        for v1 = 1:nV
            for v2 = 1:nV
                if (isempty(intersect(ind,v1))&&(~isempty(intersect(ind,v2))))||(isempty(intersect(ind,v2))&&(~isempty(intersect(ind,v1))))
                    partial_sum = partial_sum + graph_adj_matrix(v1,v2);
                end
            end
        end
        external_density_table(g,c) = partial_sum/length(ind);
    end
end
max_bound = max([max(max(external_density_table)), max(max(local_density_table))]);
figure
set(gcf,'Units','inches','Position',[2 2 13 6])
subplot(121)
imagesc(local_density_table,[0 max_bound]),colorbar
set(gca,...
    'fontS',12,...
    'XTick',[0:1:nC],...
    'YTick',[0:1:nG],...
    'FontSize',12)
xlabel('Cluster')
ylabel('Graph')
title('Internal cluster density')
subplot(122)
imagesc(external_density_table,[0 max_bound]),colorbar
set(gca,...
    'fontS',12,...
    'XTick',[0:1:nC],...
    'YTick',[0:1:nG],...
    'FontSize',12)
xlabel('Cluster')
ylabel('Graph')
title('External cluster density')
if options.saveImages
    saveas(gcf,[options.images_basename,'_cluster_densities'],'png')
end


figure
set(gcf,'Units','inches','Position',[2 2 6 6])
imagesc(local_density_table./external_density_table,[0 2]),colorbar
set(gca,...
    'fontS',12,...
    'XTick',[0:1:nC],...
    'YTick',[0:1:nG],...
    'FontSize',12)
xlabel('Cluster')
ylabel('Graph')
title('Ratio between internal and external cluster density')
if options.saveImages
    saveas(gcf,[options.images_basename,'_densities_ratio'],'png')
end
