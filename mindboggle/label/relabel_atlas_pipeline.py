#!/usr/bin/python
"""
This is a Nipype pipeline for combining brain surface atlas labels.

It will combine the surface labels in .annot format,
and convert the labels to VTK format.

https://mail.nmr.mgh.harvard.edu/pipermail//freesurfer/2010-June/014620.html
      
`mris_translate_annotation <subject> <hemi> <in annot> <translation file> <out annot>`
      
``translation file``: text file that lists the labels (one per line)
you want to group, and the new label you want to create.  You have to use
the RGB codes; each line will provide the input and output RGB values::

    221     220     60      223     220     60
    221     220     160     223     220     60
    221     220     100     223     220     60


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

import os, sys

#=============================================================================
# Setup: import libraries, set file paths, and initialize main workflow
#=============================================================================
#-----------------------------------------------------------------------------
# Steps to run
#-----------------------------------------------------------------------------
do_load_vtk_surfaces = False
do_combine_atlas_labels = 1
do_convert_atlas_annot = 1
copy_to_fs_atlases = 1
copy_to_mb_atlases = 1
#-----------------------------------------------------------------------------
# From settings.py
#-----------------------------------------------------------------------------
output_path = '/projects/Mindboggle/output'  # Where to save output
subjects_path = os.environ['SUBJECTS_DIR']  # FreeSurfer subjects directory
base_path = os.environ['MINDBOGGLE_HOME']  # Mindboggle home directory
code_path = os.environ['MINDBOGGLE_CODE']  # Mindboggle code directory
temp_path = os.path.join(output_path, 'workingdir')  # Where to save temp files
info_path = os.path.join(code_path, 'info')
sys.path.append(code_path) # Add to PYTHONPATH
label_string_old = 'labels.DKT31.manual'
label_string = 'labels.DKT25.manual'
relabel_file = os.path.join(info_path, 'labels.surface.DKT31to25.txt')
hemis = ['lh','rh']
#-----------------------------------------------------------------------------
# Subjects to process
#-----------------------------------------------------------------------------
import utils.io_file as iof
atlas_list_file = os.path.join(info_path, 'atlases101.txt')
subjects = iof.read_columns(atlas_list_file, 1)[0]
#subjects = ['MMRR-3T7T-2-1','MMRR-3T7T-2-2']
#-----------------------------------------------------------------------------
# Import system and nipype Python libraries
#-----------------------------------------------------------------------------
import os
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import Function as Fn
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import DataGrabber, DataSink
#-----------------------------------------------------------------------------
# Import Mindboggle Python libraries
#-----------------------------------------------------------------------------
from utils.io_vtk import annot_to_vtk
from label.relabel import relabel_annot_file
#-----------------------------------------------------------------------------
# Initialize main workflow
#-----------------------------------------------------------------------------
flow = Workflow(name='Atlas_relabeling_workflow')
flow.base_dir = temp_path
if not os.path.isdir(temp_path):  os.makedirs(temp_path)

#=============================================================================
#   Inputs and outputs
#=============================================================================
#-----------------------------------------------------------------------------
# Iterate inputs over subjects, hemispheres
# (surfaces are assumed to take the form: lh.pial or lh.pial.vtk)
#-----------------------------------------------------------------------------
info = Node(name = 'Inputs',
            interface = IdentityInterface(fields=['subject', 'hemi']))
info.iterables = ([('subject', subjects), ('hemi', hemis)])
#-----------------------------------------------------------------------------
# Location and structure of the surface inputs
#-----------------------------------------------------------------------------
surf = Node(name = 'Surfaces',
            interface = DataGrabber(infields=['subject', 'hemi'],
                                    outfields=['surface_files']))
surf.inputs.base_directory = subjects_path
surf.inputs.template = os.path.join('%s', 'surf', '%s.%s')
surf.inputs.template_args['surface_files'] = [['subject', 'hemi', 'pial']]
flow.connect([(info, surf, [('subject','subject'), ('hemi','hemi')])])
#-----------------------------------------------------------------------------
# Outputs
#-----------------------------------------------------------------------------
datasink = Node(DataSink(), name = 'Results')
datasink.inputs.base_directory = output_path
datasink.inputs.container = 'results'
if not os.path.isdir(output_path):  os.makedirs(output_path)

#=============================================================================
#   Combine .annot labels and convert to VTK
#=============================================================================
atlasflow = Workflow(name='Atlas_workflow')
atlasflow.base_dir = temp_path

#-------------------------------------------------------------------------
#   Combine atlas .annot labels
#-------------------------------------------------------------------------
if do_combine_atlas_labels:
    combine_labels = Node(name='Combine_atlas_labels',
                          interface = Fn(function = relabel_annot_file,
                                         input_names = ['hemi',
                                                        'subject',
                                                        'annot_name',
                                                        'new_annot_name',
                                                        'relabel_file'],
                                         output_names = ['new_annot_name']))
    combine_labels.inputs.annot_name = label_string_old
    combine_labels.inputs.new_annot_name = label_string
    combine_labels.inputs.relabel_file = relabel_file
    flow.connect([(info, atlasflow, [('hemi','Combine_atlas_labels.hemi'),
                                     ('subject','Combine_atlas_labels.subject')])])

#-----------------------------------------------------------------------------
#   Convert .annot labels to VTK format
#-----------------------------------------------------------------------------
if do_convert_atlas_annot:
    atlas_vtk = Node(name = 'Convert_atlas_labels',
                     interface = Fn(function = annot_to_vtk,
                                    input_names = ['surface_file',
                                                   'hemi',
                                                   'subject',
                                                   'subjects_path',
                                                   'annot_name'],
                                    output_names = ['vtk_file']))
    atlasflow.add_nodes([atlas_vtk])
    flow.connect([(info, atlasflow,
                        [('hemi','Convert_atlas_labels.hemi'),
                         ('subject','Convert_atlas_labels.subject')])])
    atlas_vtk.inputs.subjects_path = subjects_path
    atlasflow.connect([(combine_labels, atlas_vtk,
                     [('new_annot_name','annot_name')])])
    if do_load_vtk_surfaces:
        flow.connect([('surf', 'Convert_atlas_labels.atlas_vtk',
                       [('surface_files','Convert_atlas_labels.surface_file')])])
    else:
        flow.connect([(convertsurf, atlasflow,
                         [('vtk_file','Convert_atlas_labels.surface_file')])])
    flow.connect([(atlasflow, datasink,
                     [('Convert_atlas_labels.vtk_file','atlas_labels')])])

##############################################################################
if __name__== '__main__':
    flow.run()

#-------------------------------------------------------------------------
# Copy results to atlas label directories
#-------------------------------------------------------------------------
if copy_to_fs_atlases or copy_to_mb_atlases:

    for s in subjects:
        for h in hemis:

            srcs = []

            if do_combine_atlas_labels:
                src = os.path.join(output_path, datasink.inputs.container,
                                   'combine_labels',
                                   '_hemi_' + h + '_subject_' + s,
                                   h + '.' + label_string + '.manual.vtk')
                srcs.append(src)

            if do_convert_atlas_annot:
                src = os.path.join(output_path, datasink.inputs.container,
                                   'atlas_vtk',
                                   '_hemi_' + h + '_subject_' + s,
                                   h + '.' + label_string + '.manual.vtk')
                srcs.append(src)

            for src in srcs:
                if copy_to_fs_atlases:
                    tgt = os.path.join(atlases_path, s, 'label')
                    cmd = ' '.join(['cp', src, tgt])
                    print(cmd); os.system(cmd)
                if copy_to_mb_atlases:
                    tgt = os.path.join(base_path, 'data', 'atlases',
                                       'freesurfer', s, 'label')
                    cmd = ' '.join(['cp', src, tgt])
                    print(cmd); os.system(cmd)

