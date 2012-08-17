#!/usr/bin/python
"""
This is a Nipype pipeline for combining brain surface atlas labels.

It will convert surfaces from VTK to annot format,
combine the surface labels to create a new label set,
and convert the labels back to VTK format.

For more information about Mindboggle,
see the website: http://www.mindboggle.info
and read the README.

For information on Nipype: http://www.nipy.org/nipype/
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3159964/


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os, sys

#=============================================================================
# Setup: import libraries, set file paths, and initialize main workflow
#=============================================================================
#-----------------------------------------------------------------------------
# Steps to run
#-----------------------------------------------------------------------------
do_convert_atlas_annot = 0
do_combine_atlas_labels = 0
do_writelabels = 0
do_writeannot = 0
#-----------------------------------------------------------------------------
# From settings.py
#-----------------------------------------------------------------------------
results_path = '/projects/Mindboggle/results'  # Where to save output
subjects_path = os.environ['SUBJECTS_DIR']  # FreeSurfer subjects directory
base_path = os.environ['MINDBOGGLE_HOME']  # Mindboggle home directory
code_path = os.environ['MINDBOGGLE_CODE']  # Mindboggle code directory
temp_path = os.path.join(results_path, 'workingdir')  # Where to save temp files
info_path = os.path.join(code_path, 'info')
atlases_path = subjects_path
label_string_old = 'labels.DKT31'
label_string = 'labels.DKT25'
hemis = ['lh','rh']
sys.path.append(code_path) # Add to PYTHONPATH
# Load atlas list as subjects
import utils.io_file as iof
atlas_list_file = os.path.join(info_path, 'atlases101.txt')
subjects = iof.read_list_strings(atlas_list_file)
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
from utils.io_vtk import surf_to_vtk, annot_to_vtk
from utils.io_file import read_list_2strings
from label.relabel_surface import relabel_surface
from label.surface_labels_to_volume import write_label_file,\
     label_to_annot_file
#-----------------------------------------------------------------------------
# Initialize main workflow
#-----------------------------------------------------------------------------
outflow = Workflow(name='Outer_workflow')
outflow.base_dir = temp_path
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
# Outputs
#-----------------------------------------------------------------------------
datasink = Node(DataSink(), name = 'Results')
datasink.inputs.base_directory = results_path
datasink.inputs.container = 'output'
if not os.path.isdir(results_path):  os.makedirs(results_path)

#=============================================================================
#   Convert surfaces from VTK to annot, combine labels, and back to VTK
#=============================================================================
atlasflow = Workflow(name='Atlas_relabeling_workflow')
atlasflow.base_dir = temp_path

#-----------------------------------------------------------------------------
#   Convert VTK labels to .annot format
#-----------------------------------------------------------------------------
if do_convert_atlas_annot:
    # Input annot
    atlas_annot = Node(name = 'Atlas_annot',
                       interface = DataGrabber(infields=['subject','hemi'],
                                               outfields=['atlas_annot_file']))
    atlas_annot.inputs.base_directory = atlases_path
    atlas_annot.inputs.template = '%s/label/%s.' + label_string_old + '.manual.annot'
    atlas_annot.inputs.template_args['atlas_annot_file'] = [['subject','hemi']]
    atlasflow.add_nodes([atlas_annot])
    outflow.connect([(info, atlasflow, [('subject','Atlas_annot.subject'),
                                        ('hemi','Atlas_annot.hemi')])])
    # Convert .annot to vtk function
    atlas_vtk = Node(name = 'Convert_atlas_labels',
                     interface = Fn(function = annot_to_vtk,
                                    input_names = ['surface_file',
                                                   'hemi',
                                                   'subject',
                                                   'subjects_path',
                                                   'annot_name'],
                                    output_names = ['vtk_file']))
    atlas_vtk.inputs.annot_name = label_string_old + '.manual.annot'
    atlas_vtk.inputs.subjects_path = subjects_path
    atlasflow.add_nodes([atlas_vtk])
    outflow.connect([(info, atlasflow,
                        [('hemi','Convert_atlas_labels.hemi'),
                         ('subject','Convert_atlas_labels.subject')])])
    if do_load_vtk_surfaces:
        outflow.connect([('surf', 'Convert_atlas_labels.atlas_vtk',
                          [('surface_files','Convert_atlas_labels.surface_file')])])
    else:
        outflow.connect([(convertsurf, atlasflow,
                         [('vtk_file','Convert_atlas_labels.surface_file')])])
    outflow.connect([(atlasflow, datasink,
                     [('Convert_atlas_labels.vtk_file','atlas_labels')])])
    # Copy results to atlases label directory
    for s in subjects:
        for h in hemis:
            src = os.path.join(results_path, 'output', 'atlas_labels',
                '_hemi_' + h + '_subject_' + s,
                h + '.' + label_string_old + '.manual.vtk')
            tgt = os.path.join(atlases_path, s, 'label')
            os.system(' '.join(['cp', src, tgt]))

#-------------------------------------------------------------------------
#   Combine atlas labels
#-------------------------------------------------------------------------
if do_combine_atlas_labels:

    # After converting .annot to vtk files, define the atlas_old node
    atlas_old = Node(name = 'Atlas_old',
                     interface = DataGrabber(infields=['subject','hemi'],
                                             outfields=['atlas_old_file']))
    atlas_old.inputs.base_directory = atlases_path
    atlas_old.inputs.template = '%s/label/%s.' + label_string_old + '.manual.vtk'
    atlas_old.inputs.template_args['atlas_old_file'] = [['subject','hemi']]
    outflow.connect([(info, atlasflow, [('subject','Atlas_old.subject'),
                                        ('hemi','Atlas_old.hemi')])])
    # Combine labels
    combine_atlas = Node(name='Combine_atlas_labels',
                         interface = Fn(function = relabel_surface,
                                        input_names = ['vtk_file',
                                                       'relabel_list',
                                                       'new_string'],
                                        output_names = ['relabeled_vtk']))
    combine_atlas.inputs.new_string = label_string + '.manual.vtk'
    combine_atlas.inputs.relabel_list = os.path.join(info_path,
                                                     'labels.DKT31to25.txt')
    atlasflow.connect([(atlas_old, combine_atlas,
                        [('atlas_old_file', 'vtk_file')])])
    outflow.connect([(atlasflow, datasink,
                        [('Combine_atlas_labels.relabeled_vtk',
                          'combine_atlas_labels')])])
    # Copy results to atlases label directory
    for s in subjects:
        for h in hemis:
            src = os.path.join(results_path, 'output', 'combine_atlas_labels',
                                             '_hemi_' + h + '_subject_' + s,
                                             h + '.' + label_string + '.manual.vtk')
            tgt = os.path.join(atlases_path, s, 'label')
            os.system(' '.join(['cp', src, tgt]))

# After combining labels, define the atlas node
atlas = Node(name = 'Atlas',
             interface = DataGrabber(infields=['subject','hemi'],
                                     outfields=['atlas_file']))
atlas.inputs.base_directory = atlases_path

atlas.inputs.template = '%s/label/%s.' + label_string_old + '.manual.vtk'
atlas.inputs.template_args['atlas_file'] = [['subject','hemi']]
atlasflow.add_nodes([atlas])

outflow.connect([(info, atlasflow, [('subject','Atlas.subject'),
                                    ('hemi','Atlas.hemi')])])

#-------------------------------------------------------------------------
# Write .label files for surface vertices
#-------------------------------------------------------------------------
if do_writelabels:
    writelabels = MapNode(name='Write_label_files',
                          iterfield = ['label_number', 'label_name'],
                          interface = Fn(function = write_label_file,
                                         input_names = ['hemi',
                                                        'surface_file',
                                                        'label_number',
                                                        'label_name',
                                                        'scalar_name'],
                                         output_names = ['label_file']))

    ctx_labels_file = os.path.join(info_path, label_string + '.txt')
    ctx_label_numbers, ctx_label_names = read_list_2strings(ctx_labels_file)

    writelabels.inputs.label_number = ctx_label_numbers
    writelabels.inputs.label_name = ctx_label_names
    atlasflow.add_nodes([writelabels])
    writelabels.inputs.scalar_name = 'Labels'
    outflow.connect([(info, atlasflow, [('hemi', 'Write_label_files.hemi')])])
    atlasflow.connect([(atlas, writelabels, [('atlas_file','surface_file')])])

#-------------------------------------------------------------------------
# Write .annot file from .label files
#-------------------------------------------------------------------------
if do_writeannot:
    writeannot = MapNode(name='Write_annot_file',
                         iterfield = ['hemi'],
                         interface = Fn(function = label_to_annot_file,
                                        input_names = ['hemi',
                                                       'subjects_path',
                                                       'subject',
                                                       'label_files',
                                                       'colortable',
                                                       'annot_name'],
                                        output_names = ['annot_name',
                                                        'annot_file']))
    writeannot.inputs.annot_name = label_string + '.manual'
    writeannot.inputs.subjects_path = subjects_path
    writeannot.inputs.colortable = os.path.join(info_path, label_string + '.txt')
    atlasflow.add_nodes([writeannot])
    outflow.connect([(info, atlasflow,
                        [('hemi', 'Write_annot_file.hemi'),
                         ('subject', 'Write_annot_file.subject')])])
    atlasflow.connect([(writelabels, writeannot,
                        [('label_file','label_files')])])
    outflow.connect([(writeannot, datasink,
                      [('annot_file','Write_annot_file.annot')])])

    # Copy results to atlases label directory
    for s in subjects:
        for h in hemis:
            src = os.path.join(atlases_path, s, 'label',
                               h + '.' + label_string + '.manual.annot')
            #atlases_path2 = os.path.join(base_path, 'data/atlases/freesurfer/')
            tgt = os.path.join(atlases_path, s, 'label')
            cmd = ' '.join(['cp', src, tgt])
            print(cmd); os.system(cmd)

##############################################################################
if __name__== '__main__':

    outflow.run(plugin='Linear') #updatehash=False)
