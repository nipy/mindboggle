function a = compInputData(datapath)

addpath('/Applications/freesurfer/matlab');

%%%%%%%%%%%%%%%%%%%%%

[V,F] = read_surf([datapath '/surf/lh.inflated']);
fvInflatedLH.faces = F + 1;
fvInflatedLH.vertices = V;

[V,F] = read_surf([datapath '/surf/lh.pial']);
fvLH.faces = F+1;
fvLH.vertices = V;

DepthLH = read_curv([datapath '/surf/lh.sulc']);

%%%%%%%%%%%%%%%%%%%%%

[V,F] = read_surf([datapath '/surf/rh.inflated']);
fvInflatedRH.faces = F + 1;
fvInflatedRH.vertices = V;

[V,F] = read_surf([datapath '/surf/rh.pial']);
fvRH.faces = F+1;
fvRH.vertices = V;

DepthRH = read_curv([datapath '/surf/rh.sulc']);

%%%%%%%%%%%%%%%%%%%%%

compInputDataSave(fvLH,fvInflatedLH,DepthLH,datapath,1);
compInputDataSave(fvRH,fvInflatedRH,DepthRH,datapath,0);

a =0;
