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
subjects = ['HLN-12-3']
output_path = '/projects/Mindboggle/output'  # Where to save output

#=============================================================================
# System settings
#=============================================================================
#-----------------------------------------------------------------------------
# Labeling protocol used by Mindboggle:
# 'DKT31': 'Desikan-Killiany-Tourville (DKT) protocol with 31 labeled regions
# 'DKT25': 'fundus-friendly' version of the DKT protocol following fundi
#-----------------------------------------------------------------------------
protocol = 'DKT31'
#-----------------------------------------------------------------------------
# Initialize labels with:
# 'free': the standard FreeSurfer classifier atlas trained on the DK protocol
# <FUTURE: 'freeDKT31': a FreeSurfer-style classifier atlas trained on the DKT protocol>
# 'max': maximum probability (majority vote) labels from multiple atlases
#-----------------------------------------------------------------------------
init_labels = 'free'
#-----------------------------------------------------------------------------
# Labeling source:
# 'manual': manual edits
# <FUTURE: 'adjusted': manual edits that have been automatically aligned with major fundi>
#-----------------------------------------------------------------------------
label_method = 'manual'
#-----------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------
subjects_path = os.environ['SUBJECTS_DIR']  # FreeSurfer subjects directory
base_path = os.environ['MINDBOGGLE_HOME']  # Mindboggle home directory
code_path = os.environ['MINDBOGGLE_CODE']  # Mindboggle code directory
temp_path = os.path.join(output_path, 'workspace')  # Where to save temp files
info_path = os.path.join(code_path, 'info')
templates_path = os.path.join(base_path, 'data', 'templates')
atlases_path = subjects_path
hemis = ['lh','rh']  # Prepend ('lh.') indicating left and right surfaces
sys.path.append(code_path)  # Add to PYTHONPATH
#-----------------------------------------------------------------------------
# Debugging options
#-----------------------------------------------------------------------------
input_vtk_surfaces = False  # Load VTK surfaces (not FreeSurfer surfaces)
fill_volume_labels = True  # Fill (gray matter) volumes with surface labels
evaluate_surface_labels = 0 #False  # Compute surface overlap of auto vs. manual labels
evaluate_volume_labels = 1 #False  # Compute volume overlap of auto vs. manual labels

load_all_atlases = False  # Load atlas list as subjects
if load_all_atlases:
    import utils.io_file as iof
    atlas_list_file = os.path.join('info/atlases101.txt')
    subjects = iof.read_list_strings(atlas_list_file)

