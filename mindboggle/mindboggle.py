#!/usr/bin/python
"""
This is Mindboggle's nipype software pipeline!

Examples
--------
python mindboggle.py <output path> <1 or more subject names>
python mindboggle.py output HLN-12-1 HLN-12-2

.. note::
  Mindboggle assumes a file tree like FreeSurfer's,
  and for label initialization, assumes that subjects have been processed
  by FreeSurfer (autorecon -all), so subject names correspond to directory
  names in FreeSurfer's subjects directory.

For more information about Mindboggle,
see the website: http://www.mindboggle.info and
read the documentation: http://mindboggle.info/software/documentation.html

For information on Nipype (http://www.nipy.org/nipype/):
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3159964/


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

##############################################################################
#
#   Mindboggle workflow:
#   * Multi-atlas labeling
#   * Feature extraction
#   * Feature-based labeling
#   * Shape measurement
#
#   Followed by:
#
#   Label Volume workflow:
#   * Volume-filling labels
#   * Label evaluation
#
##############################################################################

#=============================================================================
#  User settings
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
hemis = ['lh','rh']  # Prepend ('lh.') indicating left and right surfaces
#-----------------------------------------------------------------------------
# Debugging options
#-----------------------------------------------------------------------------
input_vtk = False  # Load VTK surfaces (not FreeSurfer surfaces)
fill_volume = True  # Fill (gray matter) volumes with surface labels
include_free_measures = 0#True  # Include FreeSurfer's thickness and convexity
evaluate_surface_labels = 0 #False  # Compute surface overlap of auto vs. manual labels
evaluate_volume_labels = 1 #False  # Compute volume overlap of auto vs. manual labels

#=============================================================================
#  Setup: import libraries, set file paths, and initialize main workflow
#=============================================================================
#-----------------------------------------------------------------------------
# Import system and nipype Python libraries
#-----------------------------------------------------------------------------
import os, sys
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import Function as Fn
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import DataGrabber, DataSink
#-----------------------------------------------------------------------------
# Import Mindboggle Python libraries
#-----------------------------------------------------------------------------
"""
from mindboggle.utils.io_vtk import load_scalar, write_scalar_subset,\
    write_mean_scalar_table
from utils import io_vtk
from utils import io_file
from utils import io_free
from label import multiatlas_labeling
from label import relabel
from measure import measure_functions
from extract.fundi_hmmf import extract_folds
from extract.fundi_hmmf import extract_fundi
from label import evaluate_labels
"""
import mindboggle as mb
from mindboggle.utils.io_vtk import load_scalar, write_scalar_subset, \
     write_mean_scalar_table
from mindboggle.utils.io_file import read_columns, write_table_means
from mindboggle.utils.io_free import labels_to_annot, labels_to_volume, \
     surf_to_vtk, annot_to_vtk, vtk_to_label_files
from mindboggle.label.multiatlas_labeling import register_template,\
     transform_atlas_labels, majority_vote_label
from mindboggle.label.relabel import relabel_volume
from mindboggle.measure.measure_functions import compute_depth, \
     compute_curvature, mean_value_per_label
from mindboggle.extract.fundi_hmmf.extract_folds import extract_folds
from mindboggle.extract.fundi_hmmf.extract_fundi import extract_fundi
from mindboggle.label.evaluate_labels import measure_surface_overlap, \
     measure_volume_overlap
#-----------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------
subjects_path = os.environ['SUBJECTS_DIR']  # FreeSurfer subjects directory
atlases_path = subjects_path
templates_path = os.path.join(subjects_path, 'Mindboggle_templates')
temp_path = os.path.join(output_path, 'workspace')  # Where to save temp files
info_path = os.path.join(mindboggle.get_info()['pkg_path'], 'info')
#-----------------------------------------------------------------------------
# Initialize main workflow
#-----------------------------------------------------------------------------
mbflow = Workflow(name='Mindboggle')
mbflow.base_dir = temp_path
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
                                    outfields=['surface_files', 'sphere_files']))
surf.inputs.base_directory = subjects_path
surf.inputs.template = '%s/surf/%s.%s'
surf.inputs.template_args['surface_files'] = [['subject', 'hemi', 'pial']]
surf.inputs.template_args['sphere_files'] = [['subject', 'hemi', 'sphere']]
if include_free_measures:
    surf.inputs.template_args['thickness_files'] = [['subject', 'hemi', 'thickness']]
    surf.inputs.template_args['convexity_files'] = [['subject', 'hemi', 'curv.pial']]
mbflow.connect([(info, surf, [('subject','subject'), ('hemi','hemi')])])
#-----------------------------------------------------------------------------
# Outputs
#-----------------------------------------------------------------------------
sink = Node(DataSink(), name = 'Results')
sink.inputs.base_directory = output_path
sink.inputs.container = 'results'
if not os.path.isdir(output_path):  os.makedirs(output_path)
#-----------------------------------------------------------------------------
# Convert surfaces to VTK
#-----------------------------------------------------------------------------
if not input_vtk:
    convertsurf = Node(name = 'Free_to_VTK',
                       interface = Fn(function = surf_to_vtk,
                                      input_names = ['surface_file'],
                                      output_names = ['vtk_file']))
    mbflow.connect([(surf, convertsurf, [('surface_files','surface_file')])])
#-----------------------------------------------------------------------------
# Evaluation inputs: location and structure of atlas surfaces
#-----------------------------------------------------------------------------
if evaluate_surface_labels:
    atlas = Node(name = 'Atlases',
                 interface = DataGrabber(infields=['subject','hemi'],
                                         outfields=['atlas_file']))
    atlas.inputs.base_directory = atlases_path

    atlas.inputs.template = '%s/label/%s.labels.' +\
                            protocol + '.' + label_method + '.vtk'
    atlas.inputs.template_args['atlas_file'] = [['subject','hemi']]

    mbflow.connect([(info, atlas, [('subject','subject'),('hemi','hemi')])])
#-----------------------------------------------------------------------------
# Load data
#-----------------------------------------------------------------------------
ctx_labels_file = os.path.join(info_path, 'labels.surface.' + protocol + '.txt')
ctx_label_numbers, ctx_label_names, RGBs = read_columns(ctx_labels_file,
                                                n_columns=3, trail=True)

##############################################################################
#
#   Multi-atlas labeling workflow
#
##############################################################################
atlasflow = Workflow(name='Label_initialization')

#=============================================================================
#   Initialize labels with a classifier atlas (default to FreeSurfer labels)
#=============================================================================
if init_labels == 'free':
    freelabels = Node(name = 'Free_to_VTK',
                    interface = Fn(function = annot_to_vtk,
                                   input_names = ['surface_file',
                                                  'hemi',
                                                  'subject',
                                                  'subjects_path',
                                                  'annot_name'],
                                   output_names = ['vtk_file']))
    atlasflow.add_nodes([freelabels])
    if input_vtk:
        mbflow.connect([(surf, atlasflow,
                         [('surface_files',
                           'Free_to_VTK.surface_file')])])
    else:
        mbflow.connect([(convertsurf, atlasflow,
                         [('vtk_file',
                           'Free_to_VTK.surface_file')])])
    mbflow.connect([(info, atlasflow,
                     [('hemi', 'Free_to_VTK.hemi'),
                      ('subject', 'Free_to_VTK.subject')])])
    freelabels.inputs.subjects_path = subjects_path
    freelabels.inputs.annot_name = 'aparc.annot'

#=============================================================================
#   Initialize labels using multi-atlas registration
#=============================================================================
elif init_labels == 'max':
    #-------------------------------------------------------------------------
    # Register surfaces to average template
    #-------------------------------------------------------------------------
    template = 'OASIS-TRT-20'

    register = Node(name = 'Register_template',
                    interface = Fn(function = register_template,
                                   input_names = ['hemi',
                                                  'sphere_file',
                                                  'transform',
                                                  'templates_path',
                                                  'template'],
                                   output_names = ['transform']))
    atlasflow.add_nodes([register])
    mbflow.connect([(info, atlasflow, [('hemi', 'Register_template.hemi')]),
                    (surf, atlasflow, [('sphere_files',
                                        'Register_template.sphere_file')])])
    register.inputs.transform = 'sphere_to_' + template + '_template.reg'
    register.inputs.templates_path = os.path.join(templates_path, 'freesurfer')
    register.inputs.template = template + '.tif'
    #-------------------------------------------------------------------------
    # Register atlases to subject via template
    #-------------------------------------------------------------------------
    # Load atlas list
    atlas_list_file = os.path.join(info_path, 'atlases.txt')
    atlas_list = read_columns(atlas_list_file, 1)[0]

    transform = MapNode(name = 'Transform_labels',
                        iterfield = ['atlas'],
                        interface = Fn(function = transform_atlas_labels,
                                       input_names = ['hemi',
                                                      'subject',
                                                      'transform',
                                                      'subjects_path',
                                                      'atlas',
                                                      'atlas_string'],
                                       output_names = ['output_file']))
    atlasflow.add_nodes([transform])
    mbflow.connect([(info, atlasflow,
                     [('hemi', 'Transform_labels.hemi'),
                      ('subject', 'Transform_labels.subject')])])
    atlasflow.connect([(register, transform, [('transform', 'transform')])])
    #transform.inputs.transform = 'sphere_to_' + template + '_template.reg'
    transform.inputs.subjects_path = subjects_path
    transform.inputs.atlas = atlas_list
    transform.inputs.atlas_string = 'labels.' + protocol + '.' + label_method
    #-----------------------------------------------------------------------------
    # Majority vote label
    #-----------------------------------------------------------------------------
    vote = Node(name='Label_vote',
                interface = Fn(function = majority_vote_label,
                               input_names = ['surface_file',
                                              'annot_files'],
                               output_names = ['maxlabel_file',
                                               'labelcounts_file',
                                               'labelvotes_file',
                                               'consensus_vertices']))
    atlasflow.add_nodes([vote])
    if input_vtk:
        mbflow.connect([(surf, atlasflow,
                         [('surface_files', 'Label_vote.surface_file')])])
    else:
        mbflow.connect([(convertsurf, atlasflow,
                         [('vtk_file', 'Label_vote.surface_file')])])
    atlasflow.connect([(transform, vote, [('output_file', 'annot_files')])])
    mbflow.connect([(atlasflow, sink,
                     [('Label_vote.maxlabel_file', 'labels.@max'),
                      ('Label_vote.labelcounts_file', 'labels.@counts'),
                      ('Label_vote.labelvotes_file', 'labels.@votes')])])

##############################################################################
#
#   Surface measurement workflow
#
##############################################################################
measureflow = Workflow(name='Surface_measurement')

#=============================================================================
#   Surface measurements
#=============================================================================
#-----------------------------------------------------------------------------
# Measure surface depth
#-----------------------------------------------------------------------------
depth = Node(name='Depth',
             interface = Fn(function = compute_depth,
                            input_names = ['command',
                                           'surface_file'],
                            output_names = ['depth_file']))
depth_command = os.path.join(ccode_path,'measure', 'surface_measures',
                             'bin', 'travel_depth', 'TravelDepthMain')
depth.inputs.command = depth_command
#-----------------------------------------------------------------------------
# Measure surface curvature
#-----------------------------------------------------------------------------
curvature = Node(name='Curvature',
                 interface = Fn(function = compute_curvature,
                                input_names = ['command',
                                               'surface_file'],
                                output_names = ['mean_curvature_file',
                                                'gauss_curvature_file',
                                                'max_curvature_file',
                                                'min_curvature_file',
                                                'min_curvature_vector_file']))
curvature_command = os.path.join(ccode_path, 'measure', 'surface_measures',
                                 'bin', 'curvature', 'CurvatureMain')
curvature.inputs.command = curvature_command
#-----------------------------------------------------------------------------
# Add and connect nodes, save output files
#-----------------------------------------------------------------------------
measureflow.add_nodes([depth, curvature])
if input_vtk:
    mbflow.connect([(surf, measureflow,
                     [('surface_files','Depth.surface_file')])])
    mbflow.connect([(surf, measureflow,
                     [('surface_files','Curvature.surface_file')])])
else:
    mbflow.connect([(convertsurf, measureflow,
                     [('vtk_file', 'Depth.surface_file')])])
    mbflow.connect([(convertsurf, measureflow,
                     [('vtk_file', 'Curvature.surface_file')])])
mbflow.connect([(measureflow, sink,
                 [('Depth.depth_file', 'measures.@depth')])])
mbflow.connect([(measureflow, sink,
                 [('Curvature.mean_curvature_file',
                   'measures.@mean_curvature'),
                  ('Curvature.gauss_curvature_file',
                   'measures.@gauss_curvature'),
                  ('Curvature.max_curvature_file',
                   'measures.@max_curvature'),
                  ('Curvature.min_curvature_file',
                   'measures.@min_curvature'),
                  ('Curvature.min_curvature_vector_file',
                   'measures.@min_curvature_vectors')])])

##############################################################################
#
#   Feature extraction workflow
#
##############################################################################
featureflow = Workflow(name='Feature_extraction')

#=============================================================================
#   Feature extraction
#=============================================================================
#-----------------------------------------------------------------------------
# Extract folds
#-----------------------------------------------------------------------------
fraction_folds = 0.5
min_fold_size = 50

folds = Node(name='Folds',
             interface = Fn(function = extract_folds,
                            input_names = ['depth_file',
                                           'fraction_folds',
                                           'min_fold_size'],
                            output_names = ['folds',
                                            'n_folds',
                                            'neighbor_lists']))
featureflow.add_nodes([folds])
mbflow.connect([(measureflow, featureflow,
                 [('Depth.depth_file','Folds.depth_file')])])
folds.inputs.fraction_folds = fraction_folds
folds.inputs.min_fold_size = min_fold_size
#-----------------------------------------------------------------------------
# Extract fundi (curves at the bottoms of folds)
#-----------------------------------------------------------------------------
thr = 0.5
min_distance = 5.0

fundi = Node(name='Fundi',
             interface = Fn(function = extract_fundi,
                            input_names = ['folds',
                                           'n_folds',
                                           'neighbor_lists',
                                           'depth_file',
                                           'mean_curvature_file',
                                           'min_curvature_vector_file',
                                           'min_fold_size',
                                           'min_distance',
                                           'thr'],
                            output_names = ['fundi']))
featureflow.connect([(folds, fundi, [('folds','folds'),
                                     ('n_folds','n_folds'),
                                     ('neighbor_lists','neighbor_lists')])])
mbflow.connect([(measureflow, featureflow,
                 [('Depth.depth_file','Fundi.depth_file'),
                  ('Curvature.mean_curvature_file',
                   'Fundi.mean_curvature_file'),
                  ('Curvature.min_curvature_vector_file',
                   'Fundi.min_curvature_vector_file')])])
fundi.inputs.min_fold_size = min_fold_size
fundi.inputs.min_distance = min_distance
fundi.inputs.thr = thr
#-----------------------------------------------------------------------------
# Extract medial surfaces
#-----------------------------------------------------------------------------
#medial = Node(name='Medial',
#                 interface = Fn(function = extract_midaxis,
#                                input_names = ['depth_file',
#                                               'mean_curvature_file',
#                                               'gauss_curvature_file'],
#                                output_names = ['midaxis']))
#-----------------------------------------------------------------------------
# Write folds and fundi to VTK files
#-----------------------------------------------------------------------------
do_save_folds = True
do_save_fundi = True
if do_save_folds:
    save_folds = Node(name='Save_folds',
                      interface = Fn(function = write_scalar_subset,
                                     input_names = ['input_vtk',
                                                    'output_vtk',
                                                    'new_scalars',
                                                    'filter_scalars'],
                                     output_names = ['output_vtk']))
    featureflow.add_nodes([save_folds])
    mbflow.connect([(measureflow, featureflow,
                     [('Depth.depth_file','Save_folds.input_vtk')])])
    save_folds.inputs.output_vtk = 'folds.vtk'
    featureflow.connect([(folds, save_folds, [('folds','new_scalars')])])
    featureflow.connect([(folds, save_folds, [('folds','filter_scalars')])])
    mbflow.connect([(featureflow, sink,
                     [('Save_folds.output_vtk','features.@folds')])])
if do_save_fundi:
    save_fundi = save_folds.clone('Save_fundi')
    featureflow.add_nodes([save_fundi])
    mbflow.connect([(measureflow, featureflow,
                     [('Depth.depth_file','Save_fundi.input_vtk')])])
    save_fundi.inputs.output_vtk = 'fundi.vtk'
    featureflow.connect([(fundi, save_fundi, [('fundi','new_scalars')])])
    featureflow.connect([(folds, save_fundi, [('folds','filter_scalars')])])
    mbflow.connect([(featureflow, sink,
                     [('Save_fundi.output_vtk','features.@fundi')])])

##############################################################################
#
#   Shape analysis workflow
#
##############################################################################
shapeflow = Workflow(name='Shape_analysis')
if include_free_measures:
    column_names = ['labels', 'depth', 'mean_curvature', 'gauss_curvature',
                    'max_curvature', 'min_curvature', 'thickness', 'convexity']
else:
    column_names = ['labels', 'depth', 'mean_curvature', 'gauss_curvature',
                    'max_curvature', 'min_curvature']
measure_file_vars = [x + '_file' for x in column_names[1::]]

#-----------------------------------------------------------------------------
# Labeled surface patch shapes
#-----------------------------------------------------------------------------
input_names = ['filename', 'column_names', 'labels', 'shape_files']
input_names.extend(measure_file_vars)
"""
labeltable = Node(name='Label_table',
                      interface = Fn(function = write_table_means,
                                     input_names = input_names,
                                     output_names = ['filename']))
shapeflow.add_nodes([labeltable])
labeltable.inputs.filename = 'label_shapes.txt'
labeltable.inputs.column_names = measures
"""
#-----------------------------------------------------------------------------
# Segmented sulcus fold shapes
#-----------------------------------------------------------------------------
foldtable = Node(name='Fold_table',
                 interface = Fn(function = write_mean_scalar_table,
                                input_names = input_names,
                                output_names = ['filename']))
foldtable.inputs.column_names = column_names
#foldtable = labeltable.clone('Fold_table')
shapeflow.add_nodes([foldtable])
foldtable.inputs.filename = 'fold_shapes.txt'
mbflow.connect([(featureflow, shapeflow,
                 [('Folds.folds','Fold_table.labels')])])
mbflow.connect([(measureflow, shapeflow,
                 [('Depth.depth_file','Fold_table.depth_file')])])
mbflow.connect([(measureflow, shapeflow,
                 [('Curvature.mean_curvature_file',
                   'Fold_table.mean_curvature_file')])])
mbflow.connect([(measureflow, shapeflow,
                 [('Curvature.gauss_curvature_file',
                   'Fold_table.gauss_curvature_file')])])
mbflow.connect([(measureflow, shapeflow,
                 [('Curvature.max_curvature_file',
                   'Fold_table.max_curvature_file')])])
mbflow.connect([(measureflow, shapeflow,
                 [('Curvature.min_curvature_file',
                   'Fold_table.min_curvature_file')])])
if include_free_measures:
    mbflow.connect([(surf, shapeflow,
                     [('thickness_files', 'Fold_table.thickness_file')])])
    mbflow.connect([(surf, shapeflow,
                     [('convexity_files', 'Fold_table.convexity_file')])])
#-----------------------------------------------------------------------------
# Segmented fundus shapes
#-----------------------------------------------------------------------------
fundustable = foldtable.clone('Fundus_table')
shapeflow.add_nodes([fundustable])
fundustable.inputs.filename = 'fundus_shapes.txt'
mbflow.connect([(featureflow, shapeflow,
                 [('Fundi.fundi','Fundus_table.labels')])])
mbflow.connect([(measureflow, shapeflow,
                 [('Depth.depth_file','Fundus_table.depth_file')])])
mbflow.connect([(measureflow, shapeflow,
                 [('Curvature.mean_curvature_file',
                   'Fundus_table.mean_curvature_file')])])
mbflow.connect([(measureflow, shapeflow,
                 [('Curvature.gauss_curvature_file',
                   'Fundus_table.gauss_curvature_file')])])
mbflow.connect([(measureflow, shapeflow,
                 [('Curvature.max_curvature_file',
                   'Fundus_table.max_curvature_file')])])
mbflow.connect([(measureflow, shapeflow,
                 [('Curvature.min_curvature_file',
                   'Fundus_table.min_curvature_file')])])
if include_free_measures:
    mbflow.connect([(surf, shapeflow,
                     [('thickness_files', 'Fundus_table.thickness_file')])])
    mbflow.connect([(surf, shapeflow,
                     [('convexity_files', 'Fundus_table.convexity_file')])])

##############################################################################
#
#   Surface label evaluation
#
##############################################################################
if evaluate_surface_labels:

    #-------------------------------------------------------------------------
    # Evaluate surface labels
    #-------------------------------------------------------------------------
    eval_surf_labels = Node(name='Evaluate_surface_labels',
                            interface = Fn(function = measure_surface_overlap,
                                           input_names = ['command',
                                                          'labels_file1',
                                                          'labels_file2'],
                                           output_names = ['overlaps']))
    mbflow.add_nodes([eval_surf_labels])
    surface_overlap_command = os.path.join(ccode_path,
        'measure', 'surface_measures','bin',
        'surface_overlap', 'SurfaceOverlapMain')
    eval_surf_labels.inputs.command = surface_overlap_command
    mbflow.connect([(atlas, eval_surf_labels, [('atlas_file','labels_file1')])])
    if init_labels == 'free':
        mbflow.connect([(atlasflow, eval_surf_labels,
                         [('Free_to_VTK.vtk_file','labels_file2')])])
    elif init_labels == 'max':
        mbflow.connect([(atlasflow, eval_surf_labels,
                         [('Label_vote.maxlabel_file','labels_file2')])])
    mbflow.connect([(eval_surf_labels, sink,
                     [('overlaps', 'evaluation.@surface')])])

##############################################################################
#
#   Fill volume prep workflow:
#   Convert labels from VTK to .annot format
#
##############################################################################
if fill_volume:

    annotflow = Workflow(name='Fill_volume_prep')

    #=============================================================================
    #   Convert VTK labels to .annot format
    #=============================================================================
    #-----------------------------------------------------------------------------
    # Write .label files for surface vertices
    #-----------------------------------------------------------------------------
    writelabels = Node(name='Write_label_files',
                       interface = Fn(function = vtk_to_label_files,
                                      input_names = ['hemi',
                                                     'surface_file',
                                                     'label_numbers',
                                                     'label_names',
                                                     'RGBs',
                                                     'scalar_name'],
                                      output_names = ['label_files',
                                                      'colortable']))
    annotflow.add_nodes([writelabels])
    mbflow.connect([(info, annotflow, [('hemi', 'Write_label_files.hemi')])])
    writelabels.inputs.label_numbers = ctx_label_numbers
    writelabels.inputs.label_names = ctx_label_names
    writelabels.inputs.RGBs = RGBs
    if init_labels == 'free':
        writelabels.inputs.scalar_name = 'Labels'
        mbflow.connect([(atlasflow, annotflow,
                         [('Free_to_VTK.vtk_file',
                           'Write_label_files.surface_file')])])
    elif init_labels == 'max':
        writelabels.inputs.scalar_name = 'Max_(majority_labels)'
        mbflow.connect([(atlasflow, annotflow,
                         [('Label_vote.maxlabel_file',
                           'Write_label_files.surface_file')])])
    #-------------------------------------------------------------------------
    # Write .annot file from .label files
    # NOTE:  incorrect labels to be corrected below!
    #-------------------------------------------------------------------------
    writeannot = Node(name='Write_annot_file',
                      interface = Fn(function = labels_to_annot,
                                     input_names = ['hemi',
                                                    'subjects_path',
                                                    'subject',
                                                    'label_files',
                                                    'colortable',
                                                    'annot_name'],
                                     output_names = ['annot_name',
                                                     'annot_file']))
    writeannot.inputs.annot_name = 'temp.labels'
    writeannot.inputs.subjects_path = subjects_path
    annotflow.add_nodes([writeannot])
    mbflow.connect([(info, annotflow,
                     [('hemi', 'Write_annot_file.hemi')])])
    mbflow.connect([(info, annotflow,
                     [('subject', 'Write_annot_file.subject')])])
    annotflow.connect([(writelabels, writeannot,
                      [('label_files','label_files')])])
    annotflow.connect([(writelabels, writeannot,
                      [('colortable','colortable')])])

##############################################################################
#
#   Label volumes workflow:
#   * Fill volume
#   * Label evaluation
#
##############################################################################
mbflow2 = Workflow(name='Label_volumes')
mbflow2.base_dir = temp_path

#-----------------------------------------------------------------------------
# Iterate inputs over subjects
#-----------------------------------------------------------------------------
info2 = info.clone('Inputs2')
info2.iterables = ([('subject', subjects)])
sink2 = sink.clone('Results2')

#-----------------------------------------------------------------------------
# Fill volume mask with surface vertex labels from .annot file
#-----------------------------------------------------------------------------
if fill_volume:

    fillvolume = Node(name='Fill_volume',
                      interface = Fn(function = labels_to_volume,
                                     input_names = ['subject',
                                                    'annot_name'],
                                     output_names = ['output_file']))
    mbflow2.add_nodes([fillvolume])
    fillvolume.inputs.annot_name = 'labels.' + init_labels
    mbflow2.connect([(info2, fillvolume, [('subject', 'subject')])])
    #-------------------------------------------------------------------------
    # Relabel file, replacing colortable labels with real labels
    #-------------------------------------------------------------------------
    relabel = Node(name='Correct_labels',
                   interface = Fn(function = relabel_volume,
                                  input_names = ['input_file',
                                                 'old_labels',
                                                 'new_labels'],
                                  output_names = ['output_file']))
    mbflow2.add_nodes([relabel])
    mbflow2.connect([(fillvolume, relabel,
                          [('output_file', 'input_file')])])
    relabel_file = os.path.join(info_path,
                                'label_volume_errors.' + protocol + '.txt')
    old_labels, new_labels = read_columns(relabel_file, 2)
    relabel.inputs.old_labels = old_labels
    relabel.inputs.new_labels = new_labels
    mbflow2.connect([(relabel, sink2,
                      [('output_file', 'labels.@volume')])])

##############################################################################
#
#   Volume label evaluation workflow
#
##############################################################################
if evaluate_volume_labels:

    #-------------------------------------------------------------------------
    # Evaluation inputs: location and structure of atlas volumes
    #-------------------------------------------------------------------------
    atlas_vol = Node(name = 'Atlas',
                     interface = DataGrabber(infields=['subject'],
                     outfields=['atlas_vol_file']))
    atlas_vol.inputs.base_directory = atlases_path
#    atlas_vol.inputs.template = 'os.path.join(%s, 'mri',
#                                         'labels.' + protocol + '.nii.gz')
    atlas_vol.inputs.template = '%s/mri/aparcNMMjt+aseg.nii.gz'
    atlas_vol.inputs.template_args['atlas_vol_file'] = [['subject']]
    mbflow2.connect([(info2, atlas_vol, [('subject','subject')])])
    #-------------------------------------------------------------------------
    # Evaluate volume labels
    #-------------------------------------------------------------------------
    eval_vol_labels = Node(name='Evaluate_volume_labels',
                           interface = Fn(function = measure_volume_overlap,
                                          input_names = ['labels',
                                                         'atlas_file',
                                                         'input_file'],
                                          output_names = ['overlaps',
                                                          'out_file']))
    labels_file = os.path.join(info_path, 'labels.volume.' + protocol + '.txt')
    labels = read_columns(labels_file, 1)[0]
    eval_vol_labels.inputs.labels = labels
    mbflow2.add_nodes([eval_vol_labels])
    mbflow2.connect([(atlas_vol, eval_vol_labels,
                      [('atlas_vol_file','atlas_file')])])
    mbflow2.connect([(relabel, eval_vol_labels,
                      [('output_file', 'input_file')])])
    mbflow2.connect([(eval_vol_labels, sink2,
                      [('out_file', 'evaluation.@volume')])])

##############################################################################
#
#    Run workflow
#
##############################################################################
if __name__== '__main__':

    args = sys.argv[:]
    if len(args) < 4:
        print("\n\t Please provide the names of the output directory \
                    and one or more subject names corresponding to the names \
                    of directories within FreeSurfer's subjects directory.")
        print("\t Example: python " + args[0] + \
              " output HLN-12-1 HLN-12-2")
        sys.exit()
    else:
        output_path = str(args[1])
        subjects = str(args[2::])

    run_flow1 = True
    run_flow2 = True
    generate_graphs = True
    if generate_graphs:
        if run_flow1:
            mbflow.write_graph(graph2use='flat')
            mbflow.write_graph(graph2use='hierarchical')
        if run_flow2:
            mbflow2.write_graph(graph2use='flat')
            mbflow2.write_graph(graph2use='hierarchical')
    if run_flow1:
        mbflow.run()
    if run_flow2:
        mbflow2.run()
