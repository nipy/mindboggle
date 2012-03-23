function a = extractPitsAndSave(datapath,lambda)

addpath('/Applications/freesurfer/matlab');

[V,F] = read_surf([datapath '/surf/lh.inflated']);
fvInflatedLH.faces = F + 1;
fvInflatedLH.vertices = V;

[V,F] = read_surf([datapath '/surf/rh.inflated']);
fvInflatedRH.faces = F + 1;
fvInflatedRH.vertices = V;

load([datapath '/miccaiData/pitCandidatesRH.mat']);
extractPitsHMMF(fvLH,fvInflatedLH,DepthLH,pitIndsLH,datapath,lambda,1);
extractPitsHMMF(fvRH,fvInflatedRH,DepthRH,pitIndsRH,datapath,lambda,0);

a = 0;
