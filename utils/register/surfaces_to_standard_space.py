#!/usr/bin/python

"""
Apply FreeSurfer's "Talairach" transform matrix 
to take a conformed volume and surfaces to standard space.

Command: python <this file name> 
                <FreeSurfer subjects directory>
                <FreeSurfer subject> <output directory>

Example: python surfaces_to_standard_space.py 
                /Applications/freesurfer/subjects/ bert output/

Batch Example:
  for s in ss:
    args=['python transform_surfaces_to_standard_space.py',
          'Perrot62_sulci/freesurfer5.1_output_plus_surface_features/', s,
          'transformed_surfaces_Perrot62/' ]
    print(" ".join(args)); os.system(" ".join(args))

This program uses FreeSurfer's mris_convert. 
Example: Apply Talairach xfm to white surface, save as binary:
         mris_convert -t bert lh.white lh.white.tal


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2011  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os, sys
import nibabel as nb
from tvtk.api import tvtk

write_vtk = 1

surfaces = ['pial','smoothwm']
hemis = ['lh','rh']

# Check inputs
if len(sys.argv) < 4:
    print("python transform_surfaces_to_standard_space.py <FreeSurfer subjects directory> <FreeSurfer subject> <output directory>")
    exit(-1)
else:
    subject_path = sys.argv[1] + '/'
    subject = sys.argv[2]
    output_path = sys.argv[3] + '/'

# Apply FreeSurfer's "Talairach" transform matrix to take a surface to standard space:
if os.path.exists(subject_path) and os.path.exists(output_path):
    pass
else:
    print(subject_path + " or " + output_path + " doesn't exist.")
    exit(-1)

# Apply FreeSurfer's "Talairach" transform matrix to take conformed surfaces to standard space:
for surface in surfaces:
    for hemi in hemis:

        input_surface = subject_path + subject + '/surf/' + hemi + '.' + surface
        transformed_surface = hemi + '.' + surface + '.talairach'
        args = ['mris_convert -t', subject, input_surface, transformed_surface]
        print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()

        freesurfer_surface = subject_path + subject + '/surf/' + transformed_surface
        if os.path.exists(freesurfer_surface):
            surf = nb.freesurfer.read_geometry(freesurfer_surface)
            mesh = tvtk.PolyData(points=surf[0], polys=surf[1])
        else:
            print(freesurfer_surface + " doesn't exist.")
            exit(-1)

        # Write mesh to vtk file:
        if write_vtk:
            vtk_surface = output_path + subject + '.' + hemi + '.' + surface + '.talairach.vtp'
            w = tvtk.XMLPolyDataWriter(input=mesh, file_name=vtk_surface); w.write()


