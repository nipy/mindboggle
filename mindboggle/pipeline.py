#!/usr/bin/python
"""
This is Mindboggle's Nipype pipeline!

For more information about Mindboggle,
see the website: http://www.mindboggle.info
and read the README.

For information on Nipype: http://www.nipy.org/nipype/
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3159964/


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

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
# Setup: import libraries, set file paths, and initialize main workflow
#=============================================================================
from settings import *
#-----------------------------------------------------------------------------
# Import system and nipype Python libraries
#-----------------------------------------------------------------------------
import os
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import Function as Fn
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import DataGrabber, DataSink
#from nipype import config, logging
#config.set('logging', 'interface_level', 'DEBUG')
#config.set('logging', 'workflow_level', 'DEBUG')
#logging.update_logging(config)
#config.enable_debug_mode()
#-----------------------------------------------------------------------------
# Import Mindboggle Python libraries
#-----------------------------------------------------------------------------
from utils.io_vtk import surf_to_vtk, load_scalar, write_scalars, annot_to_vtk
from utils.io_file import read_list_strings, read_list_2strings, np_loadtxt
from label.multiatlas_labeling import register_template,\
    transform_atlas_labels, majority_vote_label
from measure.measure_functions import compute_depth, compute_curvature
from extract.fundi_hmmf.extract_folds import extract_folds
from extract.fundi_hmmf.extract_fundi import extract_fundi
from label.surface_labels_to_volume import write_label_file,\
    label_to_annot_file, fill_label_volume
from label.evaluate_volume_labels import measure_volume_overlap
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
mbflow.connect([(info, surf, [('subject','subject'), ('hemi','hemi')])])
#-----------------------------------------------------------------------------
# Outputs
#-----------------------------------------------------------------------------
sink = Node(DataSink(), name = 'Results')
sink.inputs.base_directory = output_path
sink.inputs.container = 'results'
if not os.path.isdir(output_path):  os.makedirs(output_path)
#-----------------------------------------------------------------------------
#   Convert surfaces to VTK
#-----------------------------------------------------------------------------
if not input_vtk_surfaces:
    convertsurf = Node(name = 'Convert_surfaces',
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
                            protocol + label_method + '.vtk'
    atlas.inputs.template_args['atlas_file'] = [['subject','hemi']]

    mbflow.connect([(info, atlas, [('subject','subject'),('hemi','hemi')])])

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
    freelabels = Node(name = 'Convert_labels',
                    interface = Fn(function = annot_to_vtk,
                                   input_names = ['surface_file',
                                                  'hemi',
                                                  'subject',
                                                  'subjects_path',
                                                  'annot_name'],
                                   output_names = ['vtk_file']))
    freelabels.inputs.annot_name = 'aparc.annot'
    freelabels.inputs.subjects_path = subjects_path
    atlasflow.add_nodes([freelabels])
    mbflow.connect([(info, atlasflow,
                     [('hemi', 'Convert_labels.hemi'),
                      ('subject', 'Convert_labels.subject')])])
    if input_vtk_surfaces:
        mbflow.connect([(surf, atlasflow,
                         [('surface_files',
                           'Convert_labels.surface_file')])])
    else:
        mbflow.connect([(convertsurf, atlasflow,
                         [('vtk_file',
                           'Convert_labels.surface_file')])])

#=============================================================================
#   Initialize labels using multi-atlas registration
#=============================================================================
elif init_labels == 'max':
    #-------------------------------------------------------------------------
    # Register surfaces to average template
    #-------------------------------------------------------------------------
    register = Node(name = 'Register_template',
                    interface = Fn(function = register_template,
                                   input_names = ['hemi',
                                                  'sphere_file',
                                                  'transform',
                                                  'templates_path',
                                                  'template'],
                                   output_names = ['transform']))
    template = 'OASIS-TRT-20'
    register.inputs.template = template + '.tif'
    register.inputs.transform = 'sphere_to_' + template + '_template.reg'
    register.inputs.templates_path = os.path.join(templates_path, 'freesurfer')
    atlasflow.add_nodes([register])
    mbflow.connect([(info, atlasflow, [('hemi', 'Register_template.hemi')]),
                    (surf, atlasflow, [('sphere_files',
                                        'Register_template.sphere_file')])])
    #-------------------------------------------------------------------------
    # Register atlases to subject via template
    #-------------------------------------------------------------------------
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
    # Load atlas list
    atlas_list_file = os.path.join(info_path, 'atlases.txt')
    atlas_list = read_list_strings(atlas_list_file)

    transform.inputs.atlas = atlas_list
    transform.inputs.subjects_path = subjects_path
    transform.inputs.atlas_string = 'labels.' + protocol + '.' + label_method
    atlasflow.add_nodes([transform])
    mbflow.connect([(info, atlasflow,
                     [('hemi', 'Transform_atlas_labels.hemi'),
                      ('subject', 'Transform_atlas_labels.subject')])])
    atlasflow.connect([(register, transform, [('transform', 'transform')])])
    #transform.inputs.transform = 'sphere_to_' + template + '_template.reg'
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
    if input_vtk_surfaces:
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
#   Feature extraction workflow
#
##############################################################################
featureflow = Workflow(name='Feature_extraction')

#=============================================================================
#   Surface measurements
#=============================================================================
#-----------------------------------------------------------------------------
# Measure surface depth
#-----------------------------------------------------------------------------
depth = Node(name='Compute_depth',
             interface = Fn(function = compute_depth,
                            input_names = ['command',
                                           'depth_type',
                                           'surface_file'],
                            output_names = ['depth_file']))
depth_command = os.path.join(code_path,'measure', 'surface_measures',
                             'bin', 'travel_depth', 'TravelDepthMain')
depth.inputs.command = depth_command
depth.inputs.depth_type = '0'  # 0 = travel depth, 1 = Euclidean
#-----------------------------------------------------------------------------
# Measure surface curvature
#-----------------------------------------------------------------------------
curvature = Node(name='Compute_curvature',
                 interface = Fn(function = compute_curvature,
                                input_names = ['command',
                                               'surface_file'],
                                output_names = ['mean_curvature_file',
                                                'gauss_curvature_file',
                                                'max_curvature_file',
                                                'min_curvature_file',
                                                'min_curvature_vector_file']))
curvature_command = os.path.join(code_path, 'measure', 'surface_measures',
                                 'bin', 'curvature', 'CurvatureMain')
curvature.inputs.command = curvature_command
#-----------------------------------------------------------------------------
# Add and connect nodes, save output files
#-----------------------------------------------------------------------------
featureflow.add_nodes([depth, curvature])
if input_vtk_surfaces:
    mbflow.connect([(surf, featureflow,
                     [('surface_files','Compute_depth.surface_file')])])
    mbflow.connect([(surf, featureflow,
                     [('surface_files','Compute_curvature.surface_file')])])
else:
    mbflow.connect([(convertsurf, featureflow,
                     [('vtk_file', 'Compute_depth.surface_file')])])
    mbflow.connect([(convertsurf, featureflow,
                     [('vtk_file', 'Compute_curvature.surface_file')])])
mbflow.connect([(featureflow, sink,
                 [('Compute_depth.depth_file', 'measures.@depth')])])
mbflow.connect([(featureflow, sink,
                 [('Compute_curvature.mean_curvature_file',
                   'measures.@mean_curvature'),
                  ('Compute_curvature.gauss_curvature_file',
                   'measures.@gauss_curvature'),
                  ('Compute_curvature.max_curvature_file',
                   'measures.@max_curvature'),
                  ('Compute_curvature.min_curvature_file',
                   'measures.@min_curvature'),
                  ('Compute_curvature.min_curvature_vector_file',
                   'measures.@min_curvature_vectors')])])

#=============================================================================
#   Feature extraction
#=============================================================================
#-----------------------------------------------------------------------------
# Load depth and curvature files
#-----------------------------------------------------------------------------
load_depth = Node(name='Load_depth',
                  #overwrite=True,
                  interface = Fn(function = load_scalar,
                                 input_names = ['filename'],
                                 output_names = ['Points',
                                                 'Faces',
                                                 'Scalars']))
load_curvature = Node(name='Load_curvature',
                      interface = Fn(function = load_scalar,
                                     input_names = ['filename'],
                                     output_names = ['Points',
                                                     'Faces',
                                                     'Scalars']))
load_directions = Node(name='Load_directions',
                       interface = Fn(function = np_loadtxt,
                                      input_names = ['filename'],
                                      output_names = ['output']))
featureflow.connect([(depth, load_depth,
                      [('depth_file','filename')])])
featureflow.connect([(curvature, load_curvature,
                      [('mean_curvature_file','filename')])])
featureflow.connect([(curvature, load_directions,
                      [('min_curvature_vector_file','filename')])])

#-----------------------------------------------------------------------------
# Extract folds
#-----------------------------------------------------------------------------
folds = Node(name='Extract_folds',
             interface = Fn(function = extract_folds,
                            input_names = ['faces',
                                           'depths',
                                           'fraction_folds',
                                           'min_fold_size'],
                            output_names = ['folds',
                                            'n_folds',
                                            'indices_folds',
                                            'index_lists_folds',
                                            'neighbor_lists',
                                            'faces_folds',
                                            'LUTs',
                                            'LUT_names']))
fraction_folds = 0.5
min_fold_size = 50
folds.inputs.fraction_folds = fraction_folds
folds.inputs.min_fold_size = min_fold_size
featureflow.connect([(load_depth, folds, [('Faces','faces'),
                                          ('Scalars','depths')])])
#-----------------------------------------------------------------------------
# Extract fundi (curves at the bottoms of folds)
#-----------------------------------------------------------------------------
fundi = Node(name='Extract_fundi',
             interface = Fn(function = extract_fundi,
                            input_names = ['index_lists_folds',
                                           'n_folds',
                                           'neighbor_lists',
                                           'vertices',
                                           'faces',
                                           'depths',
                                           'mean_curvatures',
                                           'min_directions',
                                           'min_fold_size',
                                           'thr',
                                           'min_distance'],
                            output_names = ['fundi',
                                            'fundus_lists',
                                            'likelihoods']))
thr = 0.5
min_distance = 5.0
fundi.inputs.thr = thr
fundi.inputs.min_fold_size = min_fold_size
fundi.inputs.min_distance = min_distance
featureflow.connect([(folds, fundi, [('index_lists_folds','index_lists_folds'),
                                     ('n_folds','n_folds'),
                                     ('neighbor_lists','neighbor_lists'),
                                     ('faces_folds','faces')]),
                     (load_depth, fundi, [('Points','vertices'),
                                          ('Scalars','depths')]),
                     (load_curvature, fundi, [('Scalars','mean_curvatures')]),
                     (load_directions, fundi, [('output','min_directions')])])
#-----------------------------------------------------------------------------
# Extract medial surfaces
#-----------------------------------------------------------------------------
#medial = Node(name='Extract_medial',
#                 interface = Fn(function = extract_midaxis,
#                                input_names = ['depth_file',
#                                               'mean_curv_file',
#                                               'gauss_curv_file'],
#                                output_names = ['midaxis']))
#-----------------------------------------------------------------------------
# Write folds, likelihoods, and fundi to VTK files
#-----------------------------------------------------------------------------
save_folds = True
save_likelihoods = False
save_fundi = True
if save_folds:
    save_folds = Node(name='Save_folds',
                      interface = Fn(function = write_scalars,
                                     input_names = ['vtk_file',
                                                    'Points',
                                                    'Vertices',
                                                    'Faces',
                                                    'LUTs',
                                                    'LUT_names'],
                                     output_names = ['vtk_file']))
    save_folds.inputs.vtk_file = 'folds.vtk'
    featureflow.connect([(load_depth, save_folds, [('Points','Points')])])
    featureflow.connect([(folds, save_folds, [('indices_folds','Vertices'),
                                              ('faces_folds','Faces'),
                                              ('LUTs','LUTs'),
                                              ('LUT_names','LUT_names')])])
    mbflow.connect([(featureflow, sink,
                     [('Save_folds.vtk_file','features.@folds')])])
if save_likelihoods:
    save_likelihoods = Node(name='Save_likelihoods',
                            interface = Fn(function = write_scalars,
                                           input_names = ['vtk_file',
                                                          'Points',
                                                          'Vertices',
                                                          'Faces',
                                                          'LUTs',
                                                          'LUT_names'],
                                           output_names = ['vtk_file']))
    save_likelihoods.inputs.vtk_file = 'likelihoods.vtk'
    save_likelihoods.inputs.LUT_names = ['likelihoods']
    featureflow.connect([(load_depth, save_likelihoods, [('Points','Points')])])
    featureflow.connect([(folds, save_likelihoods, [('indices_folds','Vertices'),
                                                    ('faces_folds','Faces')])])
    featureflow.connect([(fundi, save_likelihoods, [('likelihoods','LUTs')])])
    mbflow.connect([(featureflow, sink,
                     [('Save_likelihoods.vtk_file','measures.@likelihoods')])])
if save_fundi:
    save_fundi = Node(name='Save_fundi',
                      interface = Fn(function = write_scalars,
                                     input_names = ['vtk_file',
                                                    'Points',
                                                    'Vertices',
                                                    'Faces',
                                                    'LUTs',
                                                    'LUT_names'],
                                     output_names = ['vtk_file']))
    save_fundi.inputs.vtk_file = 'fundi.vtk'
    save_fundi.inputs.LUT_names = ['fundi']
    featureflow.connect([(load_depth, save_fundi, [('Points','Points')])])
    featureflow.connect([(folds, save_fundi, [('indices_folds','Vertices'),
                                              ('faces_folds','Faces')])])
    featureflow.connect([(fundi, save_fundi, [('fundi','LUTs')])])
    mbflow.connect([(featureflow, sink,
                     [('Save_fundi.vtk_file','features.@fundi')])])

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
                                          input_names = ['labels_file1',
                                                         'labels_file2'],
                                          output_names = ['overlaps']))
    mbflow.add_nodes([eval_surf_labels])
    mbflow.connect([(atlas, eval_surf_labels, [('atlas_file','labels_file1')])])
    if init_labels == 'free':
        mbflow.connect([(atlasflow, eval_surf_labels,
                         [('Convert_labels.vtk_file','labels_file2')])])
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
if fill_volume_labels:

    annotflow = Workflow(name='Fill_volume_prep')

    #=============================================================================
    #   Convert VTK labels to .annot format
    #=============================================================================
    #-------------------------------------------------------------------------
    # Write .label files for surface vertices
    #-------------------------------------------------------------------------
    writelabels = MapNode(name='Write_label_files',
                          iterfield = ['label_number', 'label_name'],
                          interface = Fn(function = write_label_file,
                                         input_names = ['hemi',
                                                        'surface_file',
                                                        'label_number',
                                                        'label_name',
                                                        'scalar_name'],
                                         output_names = ['label_file']))
    # List of cortical labels
    ctx_labels_file = os.path.join(info_path, 'labels.' + protocol + '.txt')
    ctx_label_numbers, ctx_label_names = read_list_2strings(ctx_labels_file)

    writelabels.inputs.label_number = ctx_label_numbers
    writelabels.inputs.label_name = ctx_label_names
    writelabels.inputs.scalar_name = 'Labels'
    annotflow.add_nodes([writelabels])
    mbflow.connect([(info, annotflow, [('hemi', 'Write_label_files.hemi')])])

    if init_labels == 'free':
        mbflow.connect([(atlasflow, annotflow,
                         [('Convert_labels.vtk_file',
                           'Write_label_files.surface_file')])])
    elif init_labels == 'max':
        mbflow.connect([(atlasflow, annotflow,
                         [('Label_vote.maxlabel_file',
                           'Write_label_files.surface_file')])])
    #-------------------------------------------------------------------------
    # Write .annot file from .label files
    #-------------------------------------------------------------------------
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
    writeannot.inputs.annot_name = 'labels.' + init_labels
    writeannot.inputs.subjects_path = subjects_path
    writeannot.inputs.colortable = os.path.join(info_path,
                                           'labels.' + protocol + '.txt')
    annotflow.add_nodes([writeannot])
    mbflow.connect([(info, annotflow,
                     [('hemi', 'Write_annot_file.hemi')])])
    mbflow.connect([(info, annotflow,
                     [('subject', 'Write_annot_file.subject')])])
    annotflow.connect([(writelabels, writeannot,
                      [('label_file','label_files')])])
    mbflow.connect([(annotflow, sink,
                     [('Write_annot_file.annot_file', 'labels.@annot')])])

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
# Evaluation inputs: location and structure of atlas volumes
#-----------------------------------------------------------------------------
if evaluate_volume_labels:

    atlas_vol = Node(name = 'Atlas_volume',
                     interface = DataGrabber(infields=['subject'],
                     outfields=['atlas_vol_file']))
    atlas_vol.inputs.base_directory = atlases_path
#    atlas_vol.inputs.template = 'os.path.join(%s, 'mri',
#                                         'labels.' + protocol + '.nii.gz')
    atlas_vol.inputs.template = '%s/mri/aparcNMMjt+aseg.nii.gz'
    atlas_vol.inputs.template_args['atlas_vol_file'] = [['subject']]
    mbflow2.connect([(info2, atlas_vol, [('subject','subject')])])

#-----------------------------------------------------------------------------
# Fill volume mask with surface vertex labels from .annot file
#-----------------------------------------------------------------------------
if fill_volume_labels:

    fillvolume = Node(name='Fill_volume_labels',
                      interface = Fn(function = fill_label_volume,
                                     input_names = ['subject',
                                                    'annot_name'],
                                     output_names = ['output_file']))
    mbflow2.add_nodes([fillvolume])
    fillvolume.inputs.annot_name = 'labels.' + init_labels
    mbflow2.connect([(info2, fillvolume, [('subject', 'subject')])])
    mbflow2.connect([(fillvolume, sink2,
                      [('output_file', 'labels.@volume')])])

##############################################################################
#
#   Volume label evaluation workflow
#
##############################################################################
if evaluate_volume_labels:

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
    labels_file = os.path.join(info_path, 'labels.' + protocol + '.txt')
    labels = read_list_strings(labels_file)
    eval_vol_labels.inputs.labels = labels
    mbflow2.add_nodes([eval_vol_labels])
    mbflow2.connect([(atlas_vol, eval_vol_labels,
                      [('atlas_vol_file','atlas_file')])])
    mbflow2.connect([(fillvolume, eval_vol_labels,
                      [('output_file', 'input_file')])])
    mbflow2.connect([(eval_vol_labels, sink2,
                      [('out_file', 'evaluation.@volume')])])

##############################################################################
#
#    Run workflow
#
##############################################################################
if __name__== '__main__':

    run_flow1 = True
    run_flow2 = True
    generate_graphs = True
    if run_flow1:
        if generate_graphs:
            mbflow.write_graph(graph2use='flat')
            mbflow.write_graph(graph2use='hierarchical')
        mbflow.run()  #(plugin='Linear') #(updatehash=False)
    if run_flow2:
        if generate_graphs:
            mbflow2.write_graph(graph2use='flat')
            mbflow2.write_graph(graph2use='hierarchical')
        mbflow2.run()  #(plugin='Linear') #(updatehash=False)
