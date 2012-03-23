function a = plotResults(fv,cmap,Q,enhance,transparency,pitInds,angle)

%transparency =1;

pitInds = pitInds(cmap(pitInds) > 0.5);

[maxCol maxColInd] = max(cmap);
disp(maxCol);
cmap(maxColInd) = 1;

figure('position',[10 10 1400 600])
p = patch(fv);

if (enhance)
    cmapNew = cmap;
    
    for i = 1:length(cmap)
        if (cmap(i) > .1)
            inds = findNeighbors(fv.faces,i);
            cmapNew(inds) = cmap(i);
        end
    end
    cmap = cmapNew;            
end

%set(p,'FaceColor','blue','EdgeColor','none');

set(p,'FaceVertexCData',cmap,'FaceColor','interp','EdgeColor','none');
colormap(jet);
%colormap(hot);

if (transparency)
    alphaMap = cmap;
    
    %alphaMap(alphaMap < .25) = 1.5*alphaMap(alphaMap < .25) - 1.5*.25;
    %alphaMap(alphaMap < 0) = 0;
    %alphaMap(alphaMap > .8) = .8;
    alphaMap(1) = 0;
    alphaMap(2) = 1;
    
    alpha(alphaMap);                
    disp(min(alphaMap))
    
    %alpha(.3);
    
    scatterMap = zeros(length(Q),1);    
    scatterMap(Q > .5) = 1;
    
    hold on;
    scatter3(fv.vertices(scatterMap ==1,1),fv.vertices(scatterMap ==1,2),fv.vertices(scatterMap ==1,3),'k','filled','SizeData',300);    
    
    scatter3(fv.vertices(pitInds,1),fv.vertices(pitInds,2),fv.vertices(pitInds,3),'b','filled','SizeData',50);
    
    %alpha('color')
    
end

daspect([1 1 1])
view(angle,0); axis off;
%axis tight
camlight
lighting gouraud;
colorbar;

zoom(1.15);
a =0 ;
