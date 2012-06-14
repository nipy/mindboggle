#!/usr/bin/python

"""
Transform FreeSurfer surfaces to the original ("native") space of an image.

Command: python <this file name> <path to subject> <output directory>

Example: python surfaces_to_native_space.py 
                /Applications/freesurfer/subjects/ output/

For volumes, see http://surfer.nmr.mgh.harvard.edu/fswiki/FsAnat-to-NativeAnat

Authors:  

* Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com
* After https://gist.github.com/1459725 by Satrajit Ghosh

(c) 2011  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os, sys
import numpy as np
import nibabel as nb
from tvtk.api import tvtk

transform_surfaces = 1
write_vtk = 1
plot_output = 0
plot_volume = 0
plot_output_files = 0

if plot_output:
    from mayavi import mlab


# Surfaces and hemispheres to transform
surfaces = ['pial','smoothwm']
hemis = ['lh','rh']

# Check inputs
if len(sys.argv) < 3:
    print("Usage: python transform_surfaces_to_native_space.py <path to FreeSurfer subject> <output directory>")
    exit(-1)
else: 
    subject_path = sys.argv[1]
    output_path = sys.argv[2]
native_volume_mgz = subject_path + '/mri/orig/001.mgz'
conformed_volume_mgz = subject_path + '/mri/brain.mgz'
if os.path.exists(native_volume_mgz) and os.path.exists(conformed_volume_mgz):
    pass
else:
    print(native_volume_mgz + " or " + conformed_volume_mgz + " doesn't exist.")
    exit(-1)

# Create and apply a transform matrix to take FreeSurfer's surface meshes to native space:
if transform_surfaces:

    # Create the transform:
    M = np.array([[-1,0,0,128],
                  [0,0,1,-128],
                  [0,-1,0,128],
                  [0,0,0,1]],dtype=float)
    native = nb.freesurfer.load(native_volume_mgz)
    conformed = nb.freesurfer.load(conformed_volume_mgz)
    affine_native = native.get_affine()
    affine_conformed = conformed.get_affine()
    xfm = np.dot(affine_conformed, np.linalg.inv(M))

    # Apply the above transform to FreeSurfer's surface meshes:
    for surface in surfaces:
        for hemi in hemis:
            freesurfer_surface = subject_path + '/surf/' + hemi + '.' + surface
            vtk_surface = output_path + '/' + hemi + '.' + surface + '.native.vtp'
            if os.path.exists(freesurfer_surface):
                surf = nb.freesurfer.read_geometry(freesurfer_surface)
                xfmd_surf = np.dot(xfm, np.hstack((surf[0], np.ones((surf[0].shape[0],1)))).T)[:3,:].T
                mesh = tvtk.PolyData(points=xfmd_surf, polys=surf[1])

                # Write mesh to vtk file:
                if write_vtk:
                    w = tvtk.XMLPolyDataWriter(input=mesh, file_name=vtk_surface); w.write()

                # Plot mesh:
                if plot_output and not plot_output_files:
                    mlab.pipeline.surface(mesh)

# Plot the brain volume and surface meshes transformed to native space:
if plot_output:

    # Plot a transformed brain volume:
    if plot_volume:
        # *Reslice a conformed volume in native ("transformed") space:
        transformed_volume_nii = output_path + '/transformed.nii.gz'
        args = ['mri_vol2vol --mov', conformed_volume_mgz, '--targ', native_volume_mgz, '--regheader --o', transformed_volume_nii]
        #args = ['mri_convert -rl',native_volume_mgz,'-rt nearest',conformed_volume_mgz,transformed_volume_nii]
        print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()
        native_data = nb.load(transformed_volume_nii).get_data()
    else:
        native_data = native.get_data()
    # Poke thanks to Isaiah Norton -- 
    # poking is simply telling an actor that you want a transformation 
    # applied to your renderer to match it with the physical world.
    mt = tvtk.Matrix4x4()
    for i in range(4):
        for j in range(4):
            mt.set_element(i,j,affine_native[i,j])
    plv = mlab.pipeline.volume(mlab.pipeline.scalar_field(native_data))
    # Poking is simply telling an actor that you want a transformation
    # applied to your renderer to match it with the physical world:
    plv.actors[0].poke_matrix(mt)
    plv.update_pipeline()

    # Plot the transformed surfaces:
    if plot_output_files:
        for surface in surfaces:
            for hemi in hemis:
                 vtk_surface = output_path + '/' + hemi + '.' + surface + '.native.vtp'
                 mesh_reader = tvtk.XMLPolyDataReader(file_name=vtk_surface)
                 mesh2 = mesh_reader.output
                 mlab.pipeline.surface(mesh2)
    mlab.show()

