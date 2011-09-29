function [spls, A, corresp] = edge_collapse_update(pts, spls, corresp, neigh, spls_adj, options)
% collapse the adjacency matrix spls_adj until the described graph is 1-dimensional.
% this inherits some ideas from [Li et al. 2001] “Decomposing  polygon meshes for interactive
% applications”, which remove all triangles of the graph to form a 1-D skeleton.
%
% 和edge_collapse做了几个对比试验，结果是一致的。
%
% pts: samples
% spls: downsamples
% corresp: correspondence between pts and spls, array of |pts|*1
% neigh: neighbors of pts, a |pts|*? cell
% spls_adj: connection matrix of downsamples (spls)
% collapse_order: 0 for Euclidean distance only; 1 for degree, then distance.
% A: new connection matrix of downsamples (spls) after edge collapse
%
% @author: JJCAO
% Changed from Andrea Tagliasacchi's code
% @data:   2010-5-7
% @version 2.0

options.null = 0;
collapse_order = getoptions(options, 'collapse_order', '0'); 

%% visual debug conditions
SHOW_COLLAPSE_PROGRESS = true;
SHOW_RESULTS = true;

if SHOW_COLLAPSE_PROGRESS || SHOW_RESULTS
    close all;
    figure(1); movegui('northwest');set(gcf,'color','white');hold on;
    plot3( spls(:,1), spls(:,2), spls(:,3), '.r', 'markersize', 5);
    axis off; axis equal;set(gcf,'Renderer','OpenGL');
end

%pause(5);
%% recover the set of edges on triangles & count triangles
A = spls_adj;
A(A>0) = 1;

degrees = ones(size(spls,1),1);
for i=1:length(spls)
    ns = find(A(i,:)==1);
    degrees(i) = length(ns)-1;
end

tricount = 0;
skeds = zeros(0,4);% idx1, idx2, average degree of two end points, distance
for i=1:length(spls)
    ns = find(A(i,:)==1);
    ns = ns( ns>i );%touch every triangle only once (if a edge belong to two triangles, it appears twice!)
    lns = length(ns);
    for j=1:lns
        for k=j+1:lns
            if A(ns(j),ns(k)) == 1
                tricount = tricount+1;
                skeds(end+1,1:3) = [i,ns(j), 0.5*(lns+degrees(ns(j))) ];                                                        %#ok<AGROW>
                skeds(end,4) = euclidean_distance(spls(i,:), spls(ns(j),:) ); 
                skeds(end+1,1:3) = [ns(j),ns(k), 0.5*(degrees(ns(j)) +degrees(ns(k)) )];                                                %#ok<AGROW>
                skeds(end,4) = euclidean_distance( spls(ns(j),:), spls(ns(k),:) );
                skeds(end+1,1:3) = [i,ns(k), 0.5*(lns+degrees(ns(k))) ];                                                     %#ok<AGROW>
                skeds(end,4) = euclidean_distance( spls(ns(k),:), spls(i,:) );
            end
        end
    end
end
    
%% --- EDGE COLLAPSE 
while true
    %DEBUG: display connectivity and vertexes at every iteration
    if SHOW_COLLAPSE_PROGRESS
        heds = [];
        for i=1:size(A,1)
            for j=1:size(A,2)
                if( A(i,j)==1 )
                    idx = [i;j];
                    heds(end+1) = line( spls(idx,1),spls(idx,2),spls(idx,3), 'LineWidth', 2, 'Color', 'b');
                end
            end
        end
        drawnow update
        pause(0.2);
        delete(heds);
    end

    disp(sprintf('decimating skeletal graph, remaining #%d loops:', size(skeds,1)));
    
    %--- STOP CONDITION
    % no more triangles? then the structure is 1D
    if size(skeds,1) == 0, break, end;
    
    %--- DECIMATION STEP + UPDATES
    % collapse the edge with minimum cost, remove the second vertex
    if collapse_order == 1 % cost is degree + distance
        mind = min( skeds(:,3) );
        tmpIdx = find(skeds(:,3)==mind);
        tmpSkeds = skeds(tmpIdx,4);

        [IGN, idx] = min( tmpSkeds );
        edge = skeds(tmpIdx(idx),1:2);
        skeds(tmpIdx(idx),:)=[];
    else % cost is distance
        [IGN, idx] = min( skeds(:,4) );
        edge = skeds(idx,1:2);
        skeds(idx,:)=[];
    end
    disp( sprintf('edge to be delete: %d, %d', edge(1),edge(2)) );

    % update the location
    spls( edge(2),: ) = mean( spls( edge,: ) );
    spls( edge(1),: ) = NaN;
    % update the A matrix
    for k=1:size(A,1)
        if A(edge(1),k) == 1, 
            A(edge(2),k)=1; 
            A(k,edge(2))=1; 
        end
    end
    % remove the row
    A(edge(1),:) = 0;
    A(:,edge(1)) = 0;
    % update the correspondents
    corresp(corresp==edge(1) ) = edge(2);
    
    %%
    % 1) remove skeds connect edge(1) and neighbor of edge(2), called 12; 
    % and skeds connect edge(2) and12; 
    % and update skeds contain edge(1) and non-neighbor of edge(2) to edge(2) and the non-neighbor of edge(2)
    tmpIdx = skeds( skeds(:,1)==edge(2), 2);% index connected with edge(2) in skeds.
    tmpIdx = [tmpIdx; skeds( skeds(:,2)==edge(2), 1)];
    
    [rows,cols] = find(skeds(:,1:2)==edge(1));
    toBeRemoved =  zeros(0,1);   
    for i = 1:length(rows)
        col = 1 + mod(cols(i),2);
        if ismember( skeds(rows(i), col), tmpIdx )%remove
            toBeRemoved(end+1) = rows(i);
        else
            skeds(rows(i), cols(i)) = edge(2);
        end
    end
    if ~isempty(toBeRemoved)
        skeds(toBeRemoved,: ) = [];
    end
    
    %% 2) remove skeds which contain edge(2) and nolonger a edge of a triangle and
    % add new triangle edges containing edge(2).
    % 2.1) find all triangles contain edge(2)
    ns = find( A(edge(2),:)==1 );
    ns = ns( ns~=edge(2) );
    lns = length(ns);
    % triangles contain edges(2) include both the two kinds of tmpEdges 
    tmpEdges = zeros(0,2); % edges contain edge(2)
    tmpEdges1 = zeros(0,2); % edges not contain edge(2)
    for j=1:lns
        for k=j+1:lns
            if A(ns(j),ns(k)) == 1
                tmpEdges(end+1,:) = [edge(2),ns(j)];  
                tmpEdges(end+1,:) = [edge(2),ns(k)];
                tmpEdges1(end+1,:) = [ns(j),ns(k)];                
            end
        end
    end
    
    % 2.2) remove all edges do not belong to tmpEdges
    [rows,cols] = find(skeds(:,1:2)==edge(2));
    toBeRemoved =  zeros(0,1);
    tobedel = zeros(0,1); 
    for j = 1:length(rows)
        col = 1 + mod(cols(j),2);
        tmp = find( tmpEdges(:,2) == skeds(rows(j),col) );
        if tmp % is an edge of tmpEdge, then remove it from tmpEdge    
            for k = 1:length(tmp)
                if ~ismember (tmp(k), tobedel)
                    tobedel = [tobedel; tmp(k)];
                end
            end
        else
            toBeRemoved(end+1) = rows(j);
        end
    end
    if ~isempty(toBeRemoved)
        skeds( toBeRemoved,: ) = [];
    end
    if ~isempty(tobedel)
        tmpEdges(tobedel,:) = [];
    end
    
    % 2.3) add triangle edges new formed
    tmpEdges = [tmpEdges; tmpEdges1];
    for j = 1:size(tmpEdges, 1)
        tedge = tmpEdges(j,:);
        [rows,cols] = find(skeds(:,1:2)==tedge(1));
        bin = false;
        for k = 1:length(rows)
            col = 1 + mod(cols(k),2);
            if skeds(rows(k),col)==tedge(2) % the edge already in skeds
                bin = true;
                break;
            end
        end
        if ~bin % add edge to skeds
            ns = find(A(tedge(1),:)==1);
            degrees(tedge(1)) = (length(ns)-1)*0.5;
            ns = find(A(tedge(2),:)==1);
            degrees(tedge(2)) = (length(ns)-1)*0.5;
            skeds(end+1,1:2) = tedge;
            skeds(end,3) = 0.5*(degrees(tedge(1))+degrees(tedge(2)));
            skeds(end,4) = euclidean_distance(spls(tedge(1),:), spls(tedge(2),:) );             
        end
    end
    
    %% 3) update distance and degree of edges contain edge(2)
    [rows,cols] = find(skeds(:,1:2)==edge(2)); 
    ns = find(A(edge(2),:)==1);
    degrees(edge(2)) = (length(ns)-1)*0.5;
    
    for j = 1:length(rows)
        col = 1 + mod(cols(j),2);   
        k = skeds(rows(j),col);
        ns = find(A(k,:)==1);
        degrees(k) = (length(ns)-1)*0.5;
        
        skeds(rows(j),3) = 0.5*(degrees(edge(2))+degrees(k));
        skeds(rows(j),4) = euclidean_distance(spls(edge(2),:), spls(k,:) ); 
    end
end
function dist = euclidean_distance(p1, p2)
v=p1-p2;
dist = sqrt(dot(v,v));