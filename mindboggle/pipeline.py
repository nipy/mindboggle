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
#   Mindboggle workflow combining:
#   * Multi-atlas labeling workflow
#   * Feature extraction workflow
#   * Feature-based labeling workflow
#   * Shape analysis workflow
#   * Label evaluation workflow
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
#from nipype import config
#config.enable_debug_mode()
#-----------------------------------------------------------------------------
# Import Mindboggle Python libraries
#-----------------------------------------------------------------------------
from utils.io_vtk import surf_to_vtk, load_scalar, write_scalars, annot_to_vtk
from utils.io_file import read_list_strings, read_list_2strings, read_surface,\
                          np_loadtxt
from label.multiatlas_labeling import register_template,\
    transform_atlas_labels, majority_vote_label
from label.relabel_surface import relabel_surface
from measure.measure_functions import compute_depth, compute_curvature
from extract.fundi_hmmf.extract_folds import extract_folds
from extract.fundi_hmmf.extract_fundi import extract_fundi
from label.surface_labels_to_volume import write_label_file,\
    label_to_annot_file, fill_label_volume
from label.evaluate_volume_labels import measure_volume_overlap
#-----------------------------------------------------------------------------
# Initialize main workflow
#-----------------------------------------------------------------------------
mbflow = Workflow(name='Mindboggle_workflow')
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
datasink = Node(DataSink(), name = 'Results')
datasink.inputs.base_directory = results_path
datasink.inputs.container = 'output'
if not os.path.isdir(results_path):  os.makedirs(results_path)
#-----------------------------------------------------------------------------
#   Convert surfaces to VTK
#-----------------------------------------------------------------------------
if not do_load_vtk_surfaces:
    convertsurf = Node(name = 'Convert_surface',
                       interface = Fn(function = surf_to_vtk,
                                      input_names = ['surface_file'],
                                      output_names = ['vtk_file']))
    mbflow.connect([(surf, convertsurf, [('surface_files','surface_file')])])
#-------------------------------------------------------------------------
#   Combine atlas labels
#-------------------------------------------------------------------------
if do_combine_atlas_labels:
    combine_atlas = Node(name='do_combine_atlas_labels',
                         interface = Fn(function = relabel_surface,
                                        input_names = ['vtk_file',
                                                       'relabel_list',
                                                       'old_string',
                                                       'new_string'],
                                        output_names = ['relabeled_vtk']))
    combine_atlas.inputs.old_string = 'labels.DKT31.manual.vtk'
    combine_atlas.inputs.new_string = label_string + '.manual.vtk'
    combine_atlas.inputs.relabel_list = os.path.join(info_path,
                                                     'labels.DKT31to25.txt')
    mbflow.connect([(convertsurf, combine_atlas,
                     [('vtk_file', 'vtk_file')])])
    mbflow.connect([(combine_atlas, datasink,
                   [('relabeled_vtk','labels.@do_combine_atlas_labels')])])
#-----------------------------------------------------------------------------
# Evaluation inputs: location and structure of atlas volumes
#-----------------------------------------------------------------------------
if do_evaluate_labels:
    atlas = Node(name = 'Atlas',
                 interface = DataGrabber(infields=['subject'],
                                         outfields=['atlas_file']))
    atlas.inputs.base_directory = atlases_path
    atlas.inputs.template = '%s/mri/aparcNMMjt+aseg.nii.gz'
    atlas.inputs.template_args['atlas_file'] = [['subject']]
    mbflow.connect([(info, atlas, [('subject','subject')])])

##############################################################################
#
#   Multi-atlas labeling workflow
#
##############################################################################
atlasflow = Workflow(name='Atlas_workflow')

#=============================================================================
#   Initialize labels with a classifier atlas (default to FreeSurfer labels)
#=============================================================================
if do_init_fs_labels:
    fslabels = Node(name = 'Convert_FreeSurfer_labels',
                    interface = Fn(function = annot_to_vtk,
                                   input_names = ['surface_file',
                                                  'hemi',
                                                  'subject',
                                                  'subjects_path',
                                                  'annot_name'],
                                   output_names = ['fslabels_file']))
    fslabels.inputs.annot_name = 'aparc.annot'
    fslabels.inputs.subjects_path = subjects_path
    atlasflow.add_nodes([fslabels])
    mbflow.connect([(info, atlasflow,
                     [('hemi', 'Convert_FreeSurfer_labels.hemi'),
                      ('subject', 'Convert_FreeSurfer_labels.subject')])])
    if do_load_vtk_surfaces:
        mbflow.connect([(surf, atlasflow,
                         [('surface_files',
                           'Convert_FreeSurfer_labels.surface_file')])])
    else:
        mbflow.connect([(convertsurf, atlasflow,
                         [('vtk_file',
                           'Convert_FreeSurfer_labels.surface_file')])])
    #-------------------------------------------------------------------------
    #   Combine labels
    #-------------------------------------------------------------------------
    if do_combine_atlas_labels:
        combine = Node(name='Combine_labels',
                         interface = Fn(function = relabel_surface,
                         input_names = ['vtk_file',
                                        'relabel_list',
                                        'old_string',
                                        'new_string'],
                         output_names = ['combine_labels_file']))
        combine.inputs.relabel_list = os.path.join(info_path,
                                                   'labels.DKT31to25.txt')
        combine.inputs.old_string = '.vtk'
        combine.inputs.new_string = label_string + '.vtk'
        atlasflow.connect([(fslabels, combine, [('fslabels_file', 'vtk_file')])])
        mbflow.connect([(atlasflow, datasink,
                         [('Combine_labels.combine_labels_file',
                           'labels.@combineFSlabels')])])

#=============================================================================
#   Initialize labels using multi-atlas registration
#=============================================================================
else:
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
    transform = MapNode(name = 'Transform_atlas_labels',
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
    transform.inputs.atlas_string = label_string + '.manual'
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
    if do_load_vtk_surfaces:
        mbflow.connect([(surf, atlasflow,
                         [('surface_files', 'Label_vote.surface_file')])])
    else:
        mbflow.connect([(convertsurf, atlasflow,
                         [('vtk_file', 'Label_vote.surface_file')])])
    atlasflow.connect([(transform, vote, [('output_file', 'annot_files')])])
    mbflow.connect([(atlasflow, datasink,
                     [('Label_vote.maxlabel_file', 'labels.@max'),
                      ('Label_vote.labelcounts_file', 'labels.@counts'),
                      ('Label_vote.labelvotes_file', 'labels.@votes')])])

##############################################################################
#
#   Feature extraction workflow
#
##############################################################################
featureflow = Workflow(name='Feature_workflow')

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
if do_load_vtk_surfaces:
    mbflow.connect([(surf, featureflow,
                     [('surface_files','Compute_depth.surface_file')])])
    mbflow.connect([(surf, featureflow,
                     [('surface_files','Compute_curvature.surface_file')])])
else:
    mbflow.connect([(convertsurf, featureflow,
                     [('vtk_file', 'Compute_depth.surface_file')])])
    mbflow.connect([(convertsurf, featureflow,
                 [('vtk_file', 'Compute_curvature.surface_file')])])
mbflow.connect([(featureflow, datasink,
                 [('Compute_depth.depth_file', 'surfaces.@depth')])])
mbflow.connect([(featureflow, datasink,
                 [('Compute_curvature.mean_curvature_file',
                   'surfaces.@mean_curvature'),
                  ('Compute_curvature.gauss_curvature_file',
                   'surfaces.@gauss_curvature'),
                  ('Compute_curvature.max_curvature_file',
                   'surfaces.@max_curvature'),
                  ('Compute_curvature.min_curvature_file',
                   'surfaces.@min_curvature'),
                  ('Compute_curvature.min_curvature_vector_file',
                   'surfaces.@min_curvature_vectors')])])

#=============================================================================
#   Feature extraction
#=============================================================================
#-----------------------------------------------------------------------------
# Load depth and curvature files
#-----------------------------------------------------------------------------
load_depth = Node(name='Load_depth',
                  interface = Fn(function = load_scalar,
                                 input_names = ['filename'],
                                 output_names = ['Points',
                                                 'Faces',
                                                 'Scalars']))
featureflow.connect([(depth, load_depth,
                      [('depth_file','filename')])])

load_curvature = Node(name='Load_curvature',
                      interface = Fn(function = load_scalar,
                                     input_names = ['filename'],
                                     output_names = ['Points',
                                                     'Faces',
                                                     'Scalars']))
featureflow.connect([(curvature, load_curvature,
                      [('mean_curvature_file','filename')])])

load_directions = Node(name='Load_directions',
                       interface = Fn(function = np_loadtxt,
                                      input_names = ['filename'],
                                      output_names = ['output']))
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

mbflow.connect([(featureflow, datasink,
                 [('Save_folds.vtk_file','folds'),
                  ('Save_likelihoods.vtk_file','likelihoods'),
                  ('Save_fundi.vtk_file','fundi')])])

##############################################################################
#
#   Label volume workflow
#
##############################################################################
if do_fill_volume_labels:

    volflow = Workflow(name='Fill_volume_workflow')

    #=============================================================================
    #   Filling a volume (e.g., gray matter) mask with labels
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
    ctx_labels_file = os.path.join(info_path, label_string + '.txt')
    ctx_label_numbers, ctx_label_names = read_list_2strings(ctx_labels_file)

    writelabels.inputs.label_number = ctx_label_numbers
    writelabels.inputs.label_name = ctx_label_names
    volflow.add_nodes([writelabels])
    mbflow.connect([(info, volflow, [('hemi', 'Write_label_files.hemi')])])

    if do_init_fs_labels:
        if do_combine_atlas_labels:
            writelabels.inputs.scalar_name = 'Max_(majority_labels)'
            mbflow.connect([(atlasflow, volflow,
                            [('Combine_labels.combine_labels_file',
                              'Write_label_files.surface_file')])])
        else:
            writelabels.inputs.scalar_name = 'Max_(majority_labels)'
            mbflow.connect([(atlasflow, volflow,
                             [('Convert_FreeSurfer_labels.fslabels_file',
                               'Write_label_files.surface_file')])])
    else:
        writelabels.inputs.scalar_name = 'Max_(majority_labels)'
        mbflow.connect([(atlasflow, volflow,
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
    if do_init_fs_labels:
        if do_combine_atlas_labels:
            label_type = label_string
        else:
            label_type = 'labels.fs'
    else:
        label_type = 'labels.max'
    writeannot.inputs.annot_name = label_type
    writeannot.inputs.subjects_path = subjects_path
    writeannot.inputs.colortable = os.path.join(info_path, label_string + '.txt')
    volflow.add_nodes([writeannot])
    mbflow.connect([(info, volflow,
                     [('hemi', 'Write_annot_file.hemi')])])
    mbflow.connect([(info, volflow,
                     [('subject', 'Write_annot_file.subject')])])
    volflow.connect([(writelabels, writeannot,
                      [('label_file','label_files')])])
    mbflow.connect([(volflow, datasink,
                     [('Write_annot_file.annot_file',
                       label_type + '.@annot')])])
    #-------------------------------------------------------------------------
    # Fill volume mask with surface vertex labels from .annot file
    #-------------------------------------------------------------------------
    fillvolume = Node(name='Fill_volume_labels',
                      interface = Fn(function = fill_label_volume,
                                     input_names = ['subject',
                                                    'annot_name'],
                                     output_names = ['output_file']))
    volflow.add_nodes([fillvolume])
    mbflow.connect([(info, volflow,
                     [('subject', 'Fill_volume_labels.subject')])])
    fillvolume.inputs.annot_name = label_type
    mbflow.connect([(volflow, datasink,
                     [('Fill_volume_labels.output_file',
                       label_type + '.@volume')])])

##############################################################################
#
#   Label evaluation
#
##############################################################################
if do_evaluate_labels:
    #-------------------------------------------------------------------------
    # Evaluate volume labels
    #-------------------------------------------------------------------------
    eval_labels = Node(name='Evaluate_volume_labels',
                       interface = Fn(function = measure_volume_overlap,
                                      input_names = ['labels',
                                                     'atlas_file',
                                                     'input_file'],
                                      output_names = ['overlaps']))
    if do_init_fs_labels:
        if do_combine_atlas_labels:
            label_type = label_string
        else:
            label_type = 'labels.fs'
    else:
        label_type = 'labels.max'
    labels_file = os.path.join(info_path, label_string + '.txt')
    labels = read_list_strings(labels_file)
    eval_labels.inputs.labels = labels
    mbflow.connect([(volflow, eval_labels,
                     [('Fill_volume_labels.output_file',
                       'Evaluate_volume_labels.input_file')])])
    mbflow.connect([(atlas, eval_labels,
                     [('atlas_file','Evaluate_volume_labels.atlas_file')])])
    mbflow.connect([(eval_labels, datasink,
                     [('Evaluate_volume_labels.overlaps',
                       label_type + '.@overlaps')])])

##############################################################################
#
#    Run workflow
#
##############################################################################
if __name__== '__main__':

    if do_generate_graph:
        mbflow.write_graph(graph2use='flat')
        mbflow.write_graph(graph2use='hierarchical')
    mbflow.run(updatehash=False)  #mbflow.run(updatehash=True)
