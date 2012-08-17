#!/usr/bin/python
"""
Settings for Mindboggle's Nipype pipeline.

Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os, sys

#=============================================================================
# User settings
#=============================================================================
subjects = ['HLN-12-3'] #['MMRR-21-1']
results_path = '/projects/Mindboggle/results'  # Where to save output

#=============================================================================
# System settings
#=============================================================================
#-----------------------------------------------------------------------------
# Debugging options
#-----------------------------------------------------------------------------
do_load_vtk_surfaces = False  # Load VTK surfaces (not FreeSurfer surfaces)
do_save_folds = True  # Save folds as VTK file
do_save_likelihoods = False  # Save likelihood values as VTK file
do_save_fundi = True  # Save fundi as VTK file
do_init_fs_labels = False  # Initialize with a FreeSurfer classifier atlas
do_fill_volume_labels = True  # Fill (gray matter) volumes with surface labels
do_evaluate_surface_labels = 0 #False  # Compute surface overlap of auto vs. manual labels
do_evaluate_volume_labels = 1 #False  # Compute volume overlap of auto vs. manual labels
#-----------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------
subjects_path = os.environ['SUBJECTS_DIR']  # FreeSurfer subjects directory
base_path = os.environ['MINDBOGGLE_HOME']  # Mindboggle home directory
code_path = os.environ['MINDBOGGLE_CODE']  # Mindboggle code directory
temp_path = os.path.join(results_path, 'workingdir')  # Where to save temp files
info_path = os.path.join(code_path, 'info')
templates_path = os.path.join(base_path, 'data', 'templates')
atlases_path = subjects_path
#label_string_old = 'labels.DKT31'
#label_string = 'labels.DKT25'
label_string = 'labels.DKT31'
hemis = ['lh','rh']

# Add to PYTHONPATH
sys.path.append(code_path)

if do_init_fs_labels:
    if do_combine_atlas_labels:
        label_type = label_string
    else:
        label_type = 'labels.fs'
else:
    label_type = 'labels.max'


# Load atlas list as subjects
import utils.io_file as iof
atlas_list_file = os.path.join('info/atlases101.txt')
subjects = iof.read_list_strings(atlas_list_file)

