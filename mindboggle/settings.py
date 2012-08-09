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
do_combine_atlas_labels = False  # Combine atlas labels
do_load_vtk_surfaces = False  # Load VTK surfaces (not FreeSurfer surfaces)
do_init_fs_labels = 0 #False  # Initialize with a FreeSurfer classifier atlas
do_fill_volume_labels = 1 #True  # Fill (gray matter) volumes with surface labels
do_evaluate_labels = 0 #False  # Compute volume overlap of auto vs. manual labels
do_generate_graph = 1 # True
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
label_string = 'labels.DKT31'
#label_string = 'labels.DKT25'
hemis = ['lh','rh']

# Add to PYTHONPATH
sys.path.append(code_path)

"""
# Load atlas list as subjects
import utils.io_file as iof
atlas_list_file = os.path.join(atlases_path, 'atlases101.txt')
subjects = iof.read_list_strings(atlas_list_file)
"""
