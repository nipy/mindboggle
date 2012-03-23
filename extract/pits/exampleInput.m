function a = exampleInput()

addpath('/Applications/freesurfer/matlab')
[V,F] = read_surf('lh.white');
fvWhite.vertices = V;
fvWhite.faces = F+1;
[V,F] = read_surf('lh.pial');
fv.faces = F+1;
fv.vertices = V;
[V,F] = read_surf('lh.inflated');
fvInflated.faces = F + 1;
fvInflated.vertices = V;
Depth = read_curv('lh.sulc');

a =0;
