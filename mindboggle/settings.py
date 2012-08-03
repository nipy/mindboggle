#!/usr/bin/python
"""
Settings for Mindboggle's Nipype pipeline.

Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

#=============================================================================
# User settings
#=============================================================================
subjects = ['HLN-12-1'] #['MMRR-21-1']
results_path = '/projects/Mindboggle/results'  # Where to save output


#=============================================================================
# System settings
#=============================================================================
#-----------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------
import os, sys
subjects_path = os.environ['SUBJECTS_DIR']  # FreeSurfer subjects directory
base_path = os.environ['MINDBOGGLE_HOME']  # Mindboggle home directory
code_path = os.environ['MINDBOGGLE_CODE']  # Mindboggle code directory
temp_path = os.path.join(results_path, 'workingdir')  # Where to save temp files
# Atlases and templates
templates_path = os.path.join(base_path, 'data', 'templates')
atlases_path = os.path.join(base_path, 'data', 'atlases')
label_string = 'labels.DKT25'
# Add to PYTHONPATH
sys.path.append(code_path)
#[sys.path.append(x) for x in ['measure','extract','label','propagate','utils']]
