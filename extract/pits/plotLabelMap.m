function a = plotLabelMap(fv,cmap,A,enhance,angle,pitInds,specialIndex)

%KKI2009-15: specialIndex = 78288
%KKI2009-16: specialIndex = 64171

transparency =0;

cmap = remapLabels(cmap);

[maxCol maxColInd] = max(cmap);
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


if (size(A,1) > 1)
    isonormals(A,p);
end
%set(p,'FaceColor','red','EdgeColor','none');
set(p,'FaceVertexCData',cmap,'FaceColor','interp','EdgeColor','none','AmbientStrength',1,'SpecularStrength',.1);
colormap(jet);
%colormap(hot);

alpha(1.0);
hold on;
scatter3(fv.vertices(pitInds,1),fv.vertices(pitInds,2),fv.vertices(pitInds,3),'k','filled','SizeData',250);
scatter3(fv.vertices(specialIndex,1),fv.vertices(specialIndex,2),fv.vertices(specialIndex,3),'k','filled','SizeData',2050,'LineWidth',10);
   

if (transparency)
    alphaMap = cmap;
    
     alphaMap(alphaMap < .1) = .1;
     alphaMap(alphaMap > .5) = 1;
     alphaMap(1) = 0;
    
    alpha(alphaMap);
    disp(min(alphaMap))
    
    scatterMap = cmap;
    scatterMap(scatterMap < .5) = 0;
    scatterMap(scatterMap > .3) = 1;
    
    hold on;
    scatter3(fv.vertices(scatterMap ==1,1),fv.vertices(scatterMap ==1,2),fv.vertices(scatterMap ==1,3),'filled');
    
    

%alpha('color')

end

daspect([1 1 1])
view(angle,0); axis off
camlight
lighting flat;
colorbar;
zoom(1.15);
a =0 ;
