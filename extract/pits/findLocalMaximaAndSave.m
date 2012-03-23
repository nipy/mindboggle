function a = findLocalMaximaAndSave(path)

addpath('/Applications/freesurfer/matlab');

[V,F] = read_surf([path '/surf/lh.pial']);
fvLH.faces = F+1;
fvLH.vertices = V;
DepthLH = read_curv([path '/surf/lh.sulc']);

[verticesLH facesLH pitIndsLH] = findPitCandidates(fvLH.faces,fvLH.vertices,DepthLH);

save([path '/miccaiData/pitCandidatesLH']);

[V,F] = read_surf([path '/surf/rh.pial']);
fvRH.faces = F+1;
fvRH.vertices = V;
DepthRH = read_curv([path '/surf/rh.sulc']);

[verticesRH facesRH pitIndsRH] = findPitCandidates(fvRH.faces,fvRH.vertices,DepthRH);

save([path '/miccaiData/pitCandidatesRH']);

a = 0;
