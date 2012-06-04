#!/usr/bin/python

"""
Save the vertices of a FreeSurfer surface mesh as an image volume.


Authors:  
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

def polydata2volume(surface_file, volume_file, output_file,
                    use_freesurfer_surfaces):
    """
    Save the vertices of a FreeSurfer surface mesh as an image volume.
    """

    from os import path
    import numpy as np
    import nibabel as nb

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

    # Create a new volume (permuted and flipped)
    V = np.zeros(vol_shape)
    for vertex in xyz:
        if use_freesurfer_surfaces:
            V[-vertex[0]+trans, -vertex[2]+trans, vertex[1]+trans] = 1
        else:
            V[vertex[0], vertex[1], vertex[2]] = 1
        
    # Save the image with the same affine transform
    img = nb.Nifti1Image(V, xfm)
    img.to_filename(output_file)


    """
    # Create a new volume (permuted and flipped)
    from apply_utils import apply_affine
    xfm2 = np.array([[-1,0,0,128],
                    [0,0,-1,128],
                    [0,1,0,128],
                    [0,0,0,1]],dtype=float)
    xyz = apply_affine(xyz[:,0], xyz[:,1], xyz[:,2], xfm2)

    V = np.zeros(vol_shape)
    for vertex in xyz:
        V[vertex[0], vertex[1], vertex[2]] = 1
    """

