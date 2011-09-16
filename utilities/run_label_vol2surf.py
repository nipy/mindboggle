
import os

voldir = "/projects/mindboggle/data/sulci_volumes/"
surfdir = "/Applications/freesurfer/subjects/Brainvisa62/"

files = os.listdir(surfdir)
for f in files: 
    for hemi in ['lh','rh']:
        if hemi == 'lh':
            hemicap = 'L'
        else:
            hemicap = 'R'
        vol = voldir+hemicap+'Bottom_'+f+'_base2008_manual.ima.nii.gz ' 
        surf = surfdir+f+'/surf/'+hemi+'.pial '
        curv = surfdir+f+'/surf/'+hemi+'.curv '
        outpt = f+'_manual.'+hemi+'.vtk '
        c = 'python label_vol2surf.py ' + vol + surf + curv + outpt
        print(c)
        os.system(c) 
