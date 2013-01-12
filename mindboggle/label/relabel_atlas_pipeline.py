#!/usr/bin/env python
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
import os

#=============================================================================
# Setup: import libraries, set file paths, and initialize main workflow
#=============================================================================
#-----------------------------------------------------------------------------
# Steps to run
#-----------------------------------------------------------------------------
do_load_vtk_surfaces = False
do_combine_atlas_labels = 1
do_convert_atlas_annot = 1
do_convert_original_atlas_annot = 1
copy_to_output = 1
#-----------------------------------------------------------------------------
# From settings.py
#-----------------------------------------------------------------------------
output_path = '/projects/Mindboggle/output'  # Where to save processing output
subjects_path = os.environ['SUBJECTS_DIR']  # FreeSurfer subjects directory
copy_path = subjects_path
base_path = '/projects/Mindboggle/mindboggle'  # Mindboggle home directory
info_path = '/projects/Mindboggle/mindboggle/mindboggle/info'  # info directory
temp_path = os.path.join(output_path, 'workspace')  # Where to save temp files
label_string_old = 'labels.DKT31.manual'
label_string = 'labels.DKT25.manual'
relabel_file = os.path.join(info_path, 'labels.surface.DKT31to25.txt')
hemis = ['lh','rh']
#-----------------------------------------------------------------------------
# Subjects to process
#-----------------------------------------------------------------------------
from mindboggle.utils.io_file import read_columns
atlas_list_file = os.path.join(info_path, 'atlases101.txt')
subjects = read_columns(atlas_list_file, 1)[0]
subjects = ['OASIS-TRT-20-11']
#-----------------------------------------------------------------------------
# Import system and nipype Python libraries
#-----------------------------------------------------------------------------
from nipype.pipeline.engine import Workflow, Node
from nipype.interfaces.utility import Function as Fn
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import DataGrabber, DataSink
#-----------------------------------------------------------------------------
# Import Mindboggle Python libraries
#-----------------------------------------------------------------------------
from mindboggle.utils.io_vtk import freeannot_to_vtk, freesurface_to_vtk
from mindboggle.label.relabel import relabel_annot_file
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

#-------------------------------------------------------------------------------
# Convert surfaces to VTK
#-------------------------------------------------------------------------------
if not do_load_vtk_surfaces:
    convertsurf = Node(name = 'Surf_to_VTK',
                       interface = Fn(function = freesurface_to_vtk,
                                      input_names = ['surface_file'],
                                      output_names = ['vtk_file']))
    flow.connect([(surf, convertsurf, [('surface_files','surface_file')])])

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
    atlasflow.add_nodes([combine_labels])
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
                     interface = Fn(function = freeannot_to_vtk,
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
#    flow.connect([(atlasflow, datasink,
#                     [('Convert_atlas_labels.vtk_file','atlas_labels')])])


if do_convert_original_atlas_annot:
    orig_atlas_vtk = Node(name = 'Convert_original_atlas_labels',
                     interface = Fn(function = freeannot_to_vtk,
                                    input_names = ['surface_file',
                                                   'hemi',
                                                   'subject',
                                                   'subjects_path',
                                                   'annot_name'],
                                    output_names = ['vtk_file']))
    atlasflow.add_nodes([orig_atlas_vtk])
    flow.connect([(info, atlasflow,
                        [('hemi','Convert_original_atlas_labels.hemi'),
                         ('subject','Convert_original_atlas_labels.subject')])])
    orig_atlas_vtk.inputs.subjects_path = subjects_path
    orig_atlas_vtk.inputs.annot_name = label_string_old
    if do_load_vtk_surfaces:
        flow.connect([('surf', 'Convert_original_atlas_labels.atlas_vtk',
                       [('surface_files','Convert_original_atlas_labels.surface_file')])])
    else:
        flow.connect([(convertsurf, atlasflow,
                         [('vtk_file','Convert_original_atlas_labels.surface_file')])])
#    flow.connect([(atlasflow, datasink,
#                     [('Convert_original_atlas_labels.vtk_file','atlas_labels')])])

##############################################################################
if __name__== '__main__':
    flow.run()

#-------------------------------------------------------------------------
# Copy results to atlas label directories
#-------------------------------------------------------------------------
if copy_to_output:

    for s in subjects:
        for h in hemis:

            if do_convert_atlas_annot:
                src = os.path.join(temp_path, #output_path, datasink.inputs.container,
                                   'Atlas_relabeling_workflow',
                                   'Atlas_workflow',
                                   '_hemi_' + h + '_subject_' + s,
                                   'Convert_atlas_labels',
                                   h + '.pial.' + label_string + '.vtk')
                tgt = os.path.join(copy_path, s, 'label', h + '.' + label_string + '.vtk')
                cmd = ' '.join(['cp', src, tgt])
                print(cmd); os.system(cmd)

            if do_convert_original_atlas_annot:
                src = os.path.join(temp_path, #output_path, datasink.inputs.container,
                                   'Atlas_relabeling_workflow',
                                   'Atlas_workflow',
                                   '_hemi_' + h + '_subject_' + s,
                                   'Convert_original_atlas_labels',
                                   h + '.pial.' + label_string_old + '.vtk')
                tgt = os.path.join(copy_path, s, 'label', h + '.' + label_string_old + '.vtk')
                cmd = ' '.join(['cp', src, tgt])
                print(cmd); os.system(cmd)

