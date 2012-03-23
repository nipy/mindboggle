function a = plotSurfWithColormapEnhance(fv,cmap,A)


figure;
p = patch(fv);

cmapNew = cmap;

for i = 1:length(cmap)    
    if (cmap(i) > .1)
        inds = findNeighbors(fv.faces,i);
        cmapNew(inds) = cmap(i);
    end
end


if (size(A,1) > 1)
    isonormals(A,p);
end
%set(p,'FaceColor','red','EdgeColor','none');
set(p,'FaceVertexCData',cmapNew,'FaceColor','interp','EdgeColor','none');


daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud;
colorbar;
a =0 ;
