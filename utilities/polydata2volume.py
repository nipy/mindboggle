import pyvtk
import nibabel as nb
import numpy as np
from apply_utils import apply_affine

volume_file = 'ribbon.nii.gz'

surface_file = '/projects/mindboggle/results/workingdir/Mindboggle_workflow/_hemi_lh_subject_id_KKI2009-11/Convert_surface/mapflow/_Convert_surface0/lh.pial_converted.labels.max.vtk'

use_freesurfer_surfaces = 1
if use_freesurfer_surfaces:
    trans = 128  # translation to middle of FreeSurfer conformed space

# Load image volume
vol = nb.load(volume_file)
vol_shape = vol.shape
xfm = vol.get_affine()

# Load surface
VTKReader = pyvtk.VtkData(surface_file)
xyz =  VTKReader.structure.points
xyz = np.array(xyz)

"""
# Create a new volume (permuted and flipped)
xfm2 = np.array([[-1,0,0,128],
                [0,0,-1,128],
                [0,1,0,128],
                [0,0,0,1]],dtype=float)
xyz = apply_affine(xyz[:,0], xyz[:,1], xyz[:,2], xfm2)

V = np.zeros(vol_shape)
for vertex in xyz:
    V[vertex[0], vertex[1], vertex[2]] = 1
"""

# Create a new volume (permuted and flipped)
V = np.zeros(vol_shape)
for vertex in xyz:
    if use_freesurfer_surfaces:
        V[-vertex[0]+trans, -vertex[2]+trans, vertex[1]+trans] = 1
    else:
        V[vertex[0], vertex[1], vertex[2]] = 1
    
# Save the image with the same affine transform
output_img = 'test.nii.gz'
img = nb.Nifti1Image(V, xfm)
img.to_filename(output_img)


