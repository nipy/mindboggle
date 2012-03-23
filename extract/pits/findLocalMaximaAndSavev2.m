function a = findLocalMaximaAndSavev2(path,LLH,LRH)

maxDist = 2.0;

addpath('/Applications/freesurfer/matlab');

[V,F] = read_surf([path '/surf/lh.pial']);
fvLH.faces = F+1;
fvLH.vertices = V;

[V,F] = read_surf([path '/surf/lh.inflated']);
fvInflatedLH.faces = F+1;
fvInflatedLH.vertices = V;

DepthLH = read_curv([path '/surf/lh.sulc']);

%[verticesLH facesLH pitIndsLH] = findPitCandidates(fvLH.faces,fvLH.vertices,DepthLH);
[pitIndsLH pitsLH] = findPitCandidatesLocalMaxWithMaxDistLikelihoodLimit(fvInflatedLH.faces,fvInflatedLH.vertices,DepthLH,LLH,maxDist);

save([path '/miccaiData/pitCandidatesLH']);

[V,F] = read_surf([path '/surf/rh.pial']);
fvRH.faces = F+1;
fvRH.vertices = V;

[V,F] = read_surf([path '/surf/rh.inflated']);
fvInflatedRH.faces = F+1;
fvInflatedRH.vertices = V;

DepthRH = read_curv([path '/surf/rh.sulc']);

%[verticesRH facesRH pitIndsRH] = findPitCandidates(fvRH.faces,fvRH.vertices,DepthRH);
[pitIndsRH pitsRH] = findPitCandidatesLocalMaxWithMaxDistLikelihoodLimit(fvInflatedRH.faces,fvInflatedRH.vertices,DepthRH,LRH,maxDist);


save([path '/miccaiData/pitCandidatesRH']);

a = 0;
