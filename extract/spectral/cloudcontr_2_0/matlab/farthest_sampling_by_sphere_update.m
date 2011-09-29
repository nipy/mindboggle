function [spls, corresp] = farthest_sampling_by_sphere(pts, RADIUS)
% use balls of size RADIUS to downsample the point cloud in a furthest sample fashion
% Speed up: we only update the distance far than %90!
% pts: points to be sampled
% neighs: neigh
% RADIUS: sample radius
% spls: sampled points
% corresp: correspondence between pts and spls, array of |pts|*1
%
% @author: JJCAO
% Changed from Andrea Tagliasacchi's code
% @data:   2010-5-7
% @version 2.0

%% visual debug conditions
SHOW_SAMPLING_PROGRESS = true;
SHOW_RESULTS = false;

%%
if SHOW_SAMPLING_PROGRESS || SHOW_RESULTS
    close all;
    figure(1); movegui('northwest');set(gcf,'color','white');hold on;
    plot3( pts(:,1), pts(:,2), pts(:,3), '.r', 'markersize', 1);
    axis off; axis equal;set(gcf,'Renderer','OpenGL');
end

% pause(5);
%%--- FURTHEST POINT DOWNSAMPLE THE CLOUD
tic
kdtree = kdtree_build( pts );
spls = zeros( 0, 3 );
corresp = zeros( length(pts), 1 );
mindst = nan( length(pts), 1 ); % mindst(i) is the min distance of pts(i) to the sample piont corresp(i) 
mindst_vect=zeros(0,2); % distance front to be updated for speeding up the algorithm 

for k=1:length(pts)
    if corresp(k)~=0, continue, end;
        
    %--- query all the points for distances
    mindst(k) = inf; % make sure picked first
    mindst_vect(end+1,:)=[inf k];
    
    %--- initialize the priority queue
    while ~all(corresp~=0) %~isempty( find(corresp==0, 1) )
        while true
            maxIdx_vect = max_idx( mindst_vect(:,1) );
            len=size(mindst_vect,1);
            if len==0
                break;
            end
            maxIdx=mindst_vect(maxIdx_vect,2);
            if mindst(maxIdx) == 0 && len>=1
                mindst_vect(maxIdx_vect,:)=[];
            else
                break;
            end
        end
        
        if mindst(maxIdx) == 0
            break
        end

        % query its delta-neighborhood
        [nIdxs, nDsts] = kdtree_ball_query( kdtree,pts(maxIdx,:), RADIUS );%original
        % if maxIdx and all its neighborhood has been marked, skip ahead
        if all( corresp(nIdxs) ~= 0 )
            mindst(maxIdx) = 0; 
            continue;
        end;

        % create new node and update (closest) distances
        spls(end+1,:) = pts(maxIdx,:); %#ok<AGROW>
        if fast
            ind=(nDsts>0.9*RADIUS);
            corresp(nIdxs(~ind))=size(spls,1);
            mindst(nIdxs(~ind))=0;
            nIdxs=nIdxs(ind);
            nDsts=nDsts(ind);
        end  
        for i=1:length(nIdxs)
            if isnan(mindst(nIdxs(i)))
                mindst(nIdxs(i)) = nDsts(i);
                corresp(nIdxs(i)) = size(spls,1);
                mindst_vect(end+1,:)=[nDsts(i) nIdxs(i)];
            elseif mindst(nIdxs(i))>nDsts(i)
                mindst(nIdxs(i)) = nDsts(i);
                corresp(nIdxs(i)) = size(spls,1);
                mindst_vect(find(mindst_vect(:,2)==nIdxs(i)),1)=nDsts(i);        
            end
        end

        if SHOW_SAMPLING_PROGRESS == true
            figure(1); plot3( pts(maxIdx,1), pts(maxIdx,2), pts(maxIdx,3), '*g');
        end
    end
end
toc
kdtree_delete( kdtree );

if SHOW_RESULTS && ~SHOW_SAMPLING_PROGRESS
    plot3( spls(:,1), spls(:,2), spls(:,3), '*g');
    hold on;axis off; axis equal;set(gcf,'Renderer','OpenGL');
end
