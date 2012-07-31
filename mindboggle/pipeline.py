#!/usr/bin/python
"""
This is Mindboggle's Nipype pipeline!

Example usage:

 # doctest: +SKIP


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

#-----------------------------------------------------------------------------
# Import Python libraries
#-----------------------------------------------------------------------------
from os import path, environ, makedirs
import sys
from re import findall
#import numpy as np
from nipype.pipeline.engine import Workflow as workflow
from nipype.pipeline.engine import Node as node
from nipype.pipeline.engine import MapNode as mapnode
from nipype.interfaces.utility import Function as fn
from nipype.interfaces.utility import IdentityInterface as identity
from nipype.interfaces.io import DataGrabber as datain
from nipype.interfaces.io import DataSink as dataout
#-----------------------------------------------------------------------------
# Options
#-----------------------------------------------------------------------------
use_freesurfer = 1
debug_skip_register = 0
do_label_volume = 1
do_evaluate_labels = 0
do_create_graph = 1
#-----------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------
subjects = ['MMRR-21-1']
subjects_path = environ['SUBJECTS_DIR']  # FreeSurfer subjects directory
basepath = '/projects/Mindboggle/mindboggle'  # mindboggle directory
mbpath = path.join(basepath, 'mindboggle')
utils_path = path.join(mbpath, 'utils')
results_path = '/projects/mindboggle/results'  # Where to save output
fundus_path = path.join(mbpath, 'extract/fundi_hmmf')
temp_path = path.join(results_path, 'workingdir')  # Where to save temp files
sys.path.append(mbpath)
sys.path.append(utils_path)
sys.path.append(fundus_path)
#-----------------------------------------------------------------------------
# Import Mindboggle Python libraries
#-----------------------------------------------------------------------------
from freesurfer2vtk import freesurfer2vtk
from io_vtk import load_vtk_map #, write_vtk_map
from atlas_functions import register_template, transform_atlas_labels,\
    majority_vote_label
from volume_functions import write_label_file, label_to_annot_file,\
    fill_label_volume, measure_volume_overlap
from surface_functions import compute_depth, compute_curvature
from extract_folds import extract_folds
#-----------------------------------------------------------------------------
# Commands (must be compiled)
#-----------------------------------------------------------------------------
depth_command = path.join(mbpath,'measure', 'surface_measures', \
                                 'bin', 'travel_depth', 'TravelDepthMain')
curvature_command = path.join(mbpath, 'measure', 'surface_measures', \
                                      'bin', 'curvature', 'CurvatureMain')

##############################################################################
#
#   Mindboggle workflow combining:
#   * Multi-atlas registration-based labeling workflow
#   * Feature-based labeling and shape analysis workflow
#   * Analytics
#
##############################################################################
mbflow = workflow(name='Mindboggle_workflow')
mbflow.base_dir = temp_path
if not path.isdir(temp_path):  makedirs(temp_path)
# Paths within mindboggle base directory
templates_path = path.join(basepath, 'data', 'templates')
atlases_path = path.join(basepath, 'data', 'atlases')
#label_string = 'labels.DKT26'
label_string = 'labels.DKT32'
atlas_string = label_string + '.manual'

##############################################################################
#   Inputs and outputs
##############################################################################
#-----------------------------------------------------------------------------
# Iterate inputs over subjects, hemispheres
# (surfaces are assumed to take the form: lh.pial or lh.pial.vtk)
#-----------------------------------------------------------------------------
info = node(name = 'Inputs',
            interface = identity(fields=['subject', 'hemi']))
info.iterables = ([('subject', subjects), ('hemi', ['lh','rh'])])
#-----------------------------------------------------------------------------
# Location and structure of the surface inputs
#-----------------------------------------------------------------------------
surf = node(name = 'Surfaces',
            interface = datain(infields=['subject', 'hemi'],
                               outfields=['surface_files', 'sphere_files']))
surf.inputs.base_directory = subjects_path
surf.inputs.template = '%s/surf/%s.%s'
surf.inputs.template_args['surface_files'] = [['subject', 'hemi', 'pial']]
surf.inputs.template_args['sphere_files'] = [['subject', 'hemi', 'sphere']]
mbflow.connect([(info, surf, [('subject','subject'), ('hemi','hemi')])])
#-----------------------------------------------------------------------------
# Location and structure of the volume inputs
#-----------------------------------------------------------------------------
"""
if do_label_volume:
    vol = node(name = 'Volume',
               interface = datain(infields=['subject', 'hemi'],
                                  outfields=['volume_file']))
    vol.inputs.base_directory = subjects_path
    vol.inputs.template = '%s/mri/%s.ribbon.mgz'
    vol.inputs.template_args['volume_file'] = [['subject', 'hemi']]
    mbflow.connect([(info, vol, [('subject','subject'), ('hemi','hemi')])])
"""
#-----------------------------------------------------------------------------
# Outputs
#-----------------------------------------------------------------------------
datasink = node(dataout(), name = 'Results')
datasink.inputs.base_directory = results_path
datasink.inputs.container = 'output'
if not path.isdir(results_path):  makedirs(results_path)

##############################################################################
#   Surface conversion to VTK
##############################################################################
#-----------------------------------------------------------------------------
# Convert FreeSurfer surfaces to VTK format
#-----------------------------------------------------------------------------
if use_freesurfer:

    convertsurf = mapnode(name = 'Convert_surface',
                          iterfield = ['in_file'],
                          interface = fn(function = freesurfer2vtk,
                          input_names = ['in_file'],
                          output_names = ['out_file']))
    mbflow.connect([(surf, convertsurf, [('surface_files','in_file')])])

    """
    if do_label_volume and use_freesurfer:
        convertvol = mapnode(name = 'Convert_volume',
            iterfield = ['in_file'],
            interface = fs.MRIConvert(out_type='niigz'))
        mbflow.connect([(vol, convertvol, [('volume_file','in_file')])])
    """

##############################################################################
#
#   Multi-atlas registration-based labeling workflow
#
##############################################################################
atlasflow = workflow(name='Atlas_workflow')

##############################################################################
#   Multi-atlas registration
##############################################################################
#-----------------------------------------------------------------------------
# Template registration
#-----------------------------------------------------------------------------
template = 'OASIS-TRT-20'
if not debug_skip_register:
    register = node(name = 'Register_template',
                    interface = fn(function = register_template,
                                   input_names = ['hemi',
                                                  'sphere_file',
                                                  'transform',
                                                  'templates_path',
                                                  'template'],
                                   output_names = ['transform']))
    register.inputs.template = template + '.tif'
    register.inputs.transform = 'sphere_to_' + template + '_template.reg'
    register.inputs.templates_path = path.join(templates_path, 'freesurfer')
    atlasflow.add_nodes([register])
    mbflow.connect([(info, atlasflow, [('hemi', 'Register_template.hemi')]),
                    (surf, atlasflow, [('sphere_files',
                                        'Register_template.sphere_file')])])
#-----------------------------------------------------------------------------
# Atlas registration
#-----------------------------------------------------------------------------
transform = mapnode(name = 'Transform_atlas_labels',
                    iterfield = ['atlas'],
                    interface = fn(function = transform_atlas_labels,
                                   input_names = ['hemi',
                                                  'subject',
                                                  'transform',
                                                  'subjects_path',
                                                  'atlas',
                                                  'atlas_string'],
                                   output_names = ['output_file']))
# List of atlas subjects
atlases_file = path.join(atlases_path, 'list_atlases.txt')
f1 = open(atlases_file)
lines1 = f1.readlines()
atlases = []
for line1 in lines1:
    atlases.append(findall(r'\S+', line1)[0])
f1.close()

transform.inputs.atlas = atlases
transform.inputs.subjects_path = subjects_path
transform.inputs.atlas_string = atlas_string
atlasflow.add_nodes([transform])
mbflow.connect([(info, atlasflow,
                 [('hemi', 'Transform_atlas_labels.hemi'),
                  ('subject', 'Transform_atlas_labels.subject')])])
if debug_skip_register:
    transform.inputs.transform = 'sphere_to_' + template + '_template.reg'
else:
    atlasflow.connect([(register, transform, [('transform', 'transform')])])
#-----------------------------------------------------------------------------
# Majority vote labeling
#-----------------------------------------------------------------------------
vote = node(name='Label_vote',
            interface = fn(function = majority_vote_label,
                           input_names = ['surface_file',
                                          'annot_files'],
                           output_names = ['maxlabel_file',
                                           'labelcounts_file',
                                           'labelvotes_file',
                                           'consensus_vertices']))
atlasflow.add_nodes([vote])
if use_freesurfer:
    mbflow.connect([(convertsurf, atlasflow,
                     [('out_file', 'Label_vote.surface_file')])])
else:
    mbflow.connect([(surf, atlasflow,
                     [('surface_files', 'Label_vote.surface_file')])])
atlasflow.connect([(transform, vote, [('output_file', 'annot_files')])])
mbflow.connect([(atlasflow, datasink,
                 [('Label_vote.maxlabel_file', 'labels.@max'),
                  ('Label_vote.labelcounts_file', 'labels.@counts'),
                  ('Label_vote.labelvotes_file', 'labels.@votes')])])

##############################################################################
#   Label propagation through a mask
##############################################################################
#-----------------------------------------------------------------------------
# Filling a volume (e.g., gray matter) mask with majority vote labels
#-----------------------------------------------------------------------------
if do_label_volume:

    #-------------------------------------------------------------------------
    # Write labels for surface vertices in .label and .annot files
    #-------------------------------------------------------------------------
    writelabels = mapnode(name='Write_label_files',
                          iterfield = ['label_number', 'label_name'],
                          interface = fn(function = write_label_file,
                                         input_names = ['hemi',
                                                        'surface_file',
                                                        'label_number',
                                                        'label_name'],
                                         output_names = ['label_file']))
    # List of cortical labels
    ctx_labels_file = path.join(atlases_path, label_string + '.txt')
    f2 = open(ctx_labels_file)
    lines2 = f2.readlines()
    ctx_label_numbers = []
    ctx_label_names = []
    for line2 in lines2:
        ctx_label_numbers.append(findall(r'\S+', line2)[0])
        ctx_label_names.append(findall(r'\S+', line2)[1])
    f2.close()
    writelabels.inputs.label_number = ctx_label_numbers
    writelabels.inputs.label_name = ctx_label_names
    atlasflow.add_nodes([writelabels])
    mbflow.connect([(info, atlasflow, [('hemi', 'Write_label_files.hemi')])])
    atlasflow.connect([(vote, writelabels, [('maxlabel_file','surface_file')])])

    # Write .annot file
    writeannot = node(name='Write_annot_file',
                      interface = fn(function = label_to_annot_file,
                                     input_names = ['hemi',
                                                    'subjects_path',
                                                    'subject',
                                                    'label_files',
                                                    'colortable'],
                                     output_names = ['annot_name',
                                                     'annot_file']))
    writeannot.inputs.subjects_path = subjects_path
    writeannot.inputs.colortable = path.join(atlases_path, label_string + '.txt')
    atlasflow.add_nodes([writeannot])
    mbflow.connect([(info, atlasflow,
                     [('hemi', 'Write_annot_file.hemi')])])
    mbflow.connect([(info, atlasflow,
                     [('subject', 'Write_annot_file.subject')])])
    atlasflow.connect([(writelabels, writeannot, [('label_file','label_files')])])
    mbflow.connect([(atlasflow, datasink,
                     [('Write_annot_file.annot_file', 'labels.@max_annot')])])

    #-------------------------------------------------------------------------
    # Fill volume mask with surface vertex labels
    #-------------------------------------------------------------------------
    fillvolume = node(name='Fill_volume_maxlabels',
                      interface = fn(function = fill_label_volume,
                                     input_names = ['subject', 'annot_name'],
                                     output_names = ['output_file']))
    atlasflow.add_nodes([fillvolume])
    mbflow.connect([(info, atlasflow,
                     [('subject', 'Fill_volume_maxlabels.subject')])])
    atlasflow.connect([(writeannot, fillvolume,
                        [('annot_name','annot_name')])])
    mbflow.connect([(atlasflow, datasink,
                     [('Fill_volume_maxlabels.output_file',
                       'labels.@maxvolume')])])
    """
    ##########################################################################
    #   Evaluation of the volume maxlabels
    ##########################################################################
    #-------------------------------------------------------------------------
    # Evaluate volume labels
    #-------------------------------------------------------------------------
    if do_evaluate_labels:
        eval_maxlabels = node(name='Evaluate_volume_maxlabels',
                              interface = fn(function = measure_volume_overlap,
                                             input_names = ['subject',
                                                            'labels',
                                                            'input_file',
                                                            'atlases_path',
                                                            'atlases',
                                                            'atlases2'],
                                             output_names = ['output_table']))
        #---------------------------------------------------------------------
        # Table with unique, non-zero labels
        #---------------------------------------------------------------------
        # List of labels
        #labels_file = path.join(atlases_path, 'labels.txt')
        #f3 = open(labels_file)
        #lines3 = f3.readlines()
        #label_numbers = []
        #for line3 in lines3:
        #    label_numbers.append(findall(r'\S+', line3)[0])
        #f3.close()
        eval_maxlabels.inputs.labels = label_numbers
        eval_maxlabels.inputs.atlases_path = atlases_path
        eval_maxlabels.inputs.atlases = atlases
        eval_maxlabels.inputs.atlases2 = atlases
        atlasflow.add_nodes([eval_maxlabels])
        mbflow.connect([(info, atlasflow, [('subject',
                                            'Evaluate_volume_maxlabels.subject')])])
        atlasflow.connect([(fill_maxlabels, eval_maxlabels,
                            [('output_file', 'input_file')])])
        mbflow.connect([(atlasflow, datasink,
                         [('Evaluate_volume_maxlabels.output_table',
                           'labels.@eval')])])
    """

##############################################################################
#
#   Feature-based labeling and shape analysis workflow
#
##############################################################################
featureflow = workflow(name='Feature_workflow')

##############################################################################
#   Surface calculations
##############################################################################
#-----------------------------------------------------------------------------
# Measure surface depth
#-----------------------------------------------------------------------------
depth = node(name='Compute_depth',
             interface = fn(function = compute_depth,
                            input_names = ['command',
                                           'surface_file'],
                            output_names = ['depth_file']))
depth.inputs.command = depth_command
#-----------------------------------------------------------------------------
# Measure surface curvature
#-----------------------------------------------------------------------------
curvature = node(name='Compute_curvature',
                 interface = fn(function = compute_curvature,
                                input_names = ['command',
                                               'surface_file'],
                                output_names = ['mean_curvature_file',
                                                'gauss_curvature_file',
                                                'max_curvature_file',
                                                'min_curvature_file',
                                                'min_curvature_vector_file']))
curvature.inputs.command = curvature_command
#-----------------------------------------------------------------------------
# Add and connect nodes
#-----------------------------------------------------------------------------
featureflow.add_nodes([depth, curvature])
if use_freesurfer:
    mbflow.connect([(convertsurf, featureflow,
                   [('out_file', 'Compute_depth.surface_file')])])
    mbflow.connect([(convertsurf, featureflow,
                   [('out_file', 'Compute_curvature.surface_file')])])
else:
    # Connect input to surface depth and curvature nodes
    mbflow.connect([(atlasflow, featureflow,
                     [('Surfaces.surface_files',
                       'Compute_depth.surface_file')])])
    mbflow.connect([(atlasflow, featureflow,
                     [('Surfaces.surface_files',
                       'Compute_curvature.surface_file')])])
#-----------------------------------------------------------------------------
# Save
#-----------------------------------------------------------------------------
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

##############################################################################
#   Feature extraction
##############################################################################
#-----------------------------------------------------------------------------
# Load depth file
#-----------------------------------------------------------------------------
depth_load = node(name='Load_depth',
                  interface = fn(function = load_vtk_map,
                                 input_names = ['filename'],
                                 output_names = ['vertices, faces, scalars']))
#featureflow.add_nodes([depth_load])
featureflow.connect([(depth, depth_load, [('depth_file','filename')])])
#-----------------------------------------------------------------------------
# Extract folds
#-----------------------------------------------------------------------------
folds = node(name='Extract_folds',
             interface = fn(function = extract_folds,
                            input_names = ['faces',
                                           'depths',
                                           'fraction_folds',
                                           'min_fold_size'],
                            output_names = ['folds',
                                            'n_folds',
                                            'index_lists_folds',
                                            'neighbor_lists']))
folds.inputs.fraction_folds = 0.5
folds.inputs.min_fold_size = 50
featureflow.add_nodes([folds])
featureflow.connect([(depth_load, folds, [('faces','faces'),
                                          ('scalars','depths')])])
#-----------------------------------------------------------------------------
# Save folds
#-----------------------------------------------------------------------------
"""
folds_save = node(name='Load_features',
                  interface = fn(function = write_Sulci,
                                 input_names = ['depth_file'],
                                 output_names = ['vertices, faces, depths']))
featureflow.connect([(depth, depth_load, [('depth_file','depth_file')])])



indices_folds = [x for lst in index_lists_folds for x in lst]
# Remove faces that do not contain three fold vertices
fs = frozenset(indices_folds)
faces_folds = [lst for lst in faces if len(fs.intersection(lst)) == 3]
faces_folds = np.reshape(np.ravel(faces_folds), (-1, 3))
print('  Reduced {} to {} faces.'.format(len(faces),
                                         len(faces_folds)))
# Save vtk file
folds_for_vtk = folds.copy()
folds_for_vtk[folds == 0] = -1
LUTs = [[int(x) for x in folds_for_vtk]]
LUT_names = ['fold'+str(i+1) for i in range(n_folds)]
io_vtk.writeSulci(load_path + 'folds.vtk', vertices, indices_folds,
                  faces_folds, LUTs=LUTs, LUTNames=LUT_names)


# Remove faces that do not contain three fold vertices
indices_folds = [x for lst in index_lists_folds for x in lst]
fs = frozenset(indices_folds)
faces_folds = [lst for lst in faces if len(fs.intersection(lst)) == 3]
faces_folds = np.reshape(np.ravel(faces_folds), (-1, 3))

    # Save fold likelihoods
        io_vtk.writeSulci(load_path + 'likelihoods.vtk', vertices,
                          indices_folds, faces_folds,
                          LUTs=[likelihoods],
                          LUTNames=['likelihoods'])
    # Save fundi
        fundi_for_vtk = np.ones(n_vertices)
        for fundus in fundi:
            if len(fundus) > 0:
                fundi_for_vtk += fundus
        io_vtk.writeSulci(load_path + 'fundi.vtk', vertices,
            indices_folds, faces_folds,
            LUTs=[fundi_for_vtk], LUTNames=['fundi'])
"""
"""
#-----------------------------------------------------------------------------
# Load curvature file
#-----------------------------------------------------------------------------
depth_load = node(name='Load_features',
                interface = fn(function = load_VTK_Map,
                            input_names = ['depth_file'],
                            output_names = ['vertices, faces, depths']))
featureflow.connect([(depth, depth_load, [('depth_file','depth_file')])])
#-----------------------------------------------------------------------------
# Extract fundi (curves at the bottoms of folds
#-----------------------------------------------------------------------------
curv_file = load_path + 'lh.pial.curv.avg.vtk'
dir_file = load_path + 'lh.pial.curv.min.dir.csv'
vertices, faces, mean_curvatures = io_vtk.load_VTK_Map(curv_file)
vertices = np.array(vertices)
mean_curvatures = np.array(mean_curvatures)
min_directions = np.loadtxt(dir_file)

fundi = node(name='Extract_fundi',
             interface = fn(function = extract_fundi,
                            input_names = [index_lists_folds, n_folds, neighbor_lists,
                                              vertices, faces, depths, mean_curvatures, min_directions,
                                              min_fold_size, thr, min_distance]
                            output_names = ['fundi']))
fundi.inputs.command = extract_fundi_command

"""
"""
#-----------------------------------------------------------------------------
# Extract medial surfaces
#-----------------------------------------------------------------------------
medial = node(name='Extract_medial',
                 interface = fn(
                                  function = extract_midaxis,
                                  input_names = ['depth_file',
                                                 'mean_curv_file',
                                                 'gauss_curv_file'],
                                  output_names = ['midaxis']))
"""
"""
#-----------------------------------------------------------------------------
# Connect surface depth to feature extraction nodes
#-----------------------------------------------------------------------------
featureflow.connect([(depth, fundi,
               [('depth_file', 'depth_file')])])
featureflow.connect([(fundi, datasink,
               [('fundi', 'fundi')])])

featureflow.connect([(surfaces, sulcus_extraction,
               [('depth_file', 'depth_file'),
                ('mean_curv_file', 'mean_curv_file'),
                ('gauss_curv_file', 'gauss_curv_file')])])
featureflow.connect([(surfaces, midaxis_extraction,
               [('depth_file', 'depth_file'),
                ('mean_curv_file', 'mean_curv_file'),
                ('gauss_curv_file', 'gauss_curv_file')])])

##############################################################################
#   Label propagation
##############################################################################
#-----------------------------------------------------------------------------
# Label propagation node
#-----------------------------------------------------------------------------
propagate = node(name='Propagate_labels',
                    interface = fn(
                                     function = propagate_labels,
                                     input_names=['labels', 'fundi'],
                                     output_names=['labels']))

#-----------------------------------------------------------------------------
# Volume label propagation node
#-----------------------------------------------------------------------------
propagate_volume = node(name='Propagate_volume_labels',
                           interface = fn(
                                            function = propagate_volume_labels,
                                            input_names = ['labels'],
                                            output_names = ['volume_labels']))

#-----------------------------------------------------------------------------
# Labeled surface patch and volume extraction nodes
#-----------------------------------------------------------------------------
extract_patch = node(name='Extract_patch',
                        interface = fn(
                                         function = extract_patches,
                                         input_names=['labels'],
                                         output_names=['patches']))

extract_patch = node(name='Extract_region',
                        interface = fn(
                                         function = extract_patches,
                                         input_names=['labels'],
                                         output_names=['patches']))

#-----------------------------------------------------------------------------
# Connect registration and extraction to label propagation nodes
#-----------------------------------------------------------------------------
featureflow.connect([(transform, propagate, [('labels','labels')]),
                     (extract_fundi, propagate, [('fundi','fundi')])])

# Connect label propagation to labeled surface patch and volume extraction nodes
featureflow.connect([(propagate, propagate_volume, [('labels', 'labels')])])

featureflow.connect([(propagate_volume, extract_region, [('labels', 'labels')])])
featureflow.connect([(propagate, extract_patch, [('labels', 'labels')])])

##############################################################################
#   Feature segmentation / identification
##############################################################################
#-----------------------------------------------------------------------------
# Feature segmentation nodes
#-----------------------------------------------------------------------------
segment_sulci = node(name='Segment_sulci',
                        interface = fn(
                                         function = segment_sulci,
                                         input_names=['sulci','labels'],
                                         output_names=['segmented_sulci']))

segment_fundi = node(name='Segment_fundi',
                        interface = fn(
                                         function = segment_fundi,
                                         input_names=['fundi','labels'],
                                         output_names=['segmented_fundi']))

segment_medial = node(name='Segment_medial',
                         interface = fn(
                                          function = segment_medial,
                                          input_names=['medial','labels'],
                                          output_names=['segmented_medial']))

# Connect feature and feature segmentation nodes
featureflow.connect([(extract_sulci, segment_sulci, [('sulci','sulci')]),
              (extract_fundi, segment_fundi, [('fundi','fundi')]),
              (extract_medial, segment_medial, [('medial','medial')])])

# Connect label propagation and feature segmentation nodes
featureflow.connect([(propagate, segment_sulci, [('labels','labels')]),
              (propagate, segment_fundi, [('labels','labels')]),
              (segment_sulci, segment_medial, [('segmented_sulci','labels')])])

##############################################################################
#   Shape measurement
##############################################################################
#-----------------------------------------------------------------------------
# Shape measurement nodes
#-----------------------------------------------------------------------------
positions = node(interface = fn(input_names = ['segmented_sulci',
                                                 'segmented_fundi',
                                                 'segmented_midaxis',
                                                 'pits',
                                                 'patches',
                                                 'regions'],
                                  output_names=['positions_sulci',
                                                'positions_fundi',
                                                'positions_midaxis',
                                                'positions_pits',
                                                'positions_patches',
                                                'positions_regions'],
                                  function = measure_positions),
                    name='Measure_positions')

extents = node(interface = fn(input_names = ['segmented_sulci',
                                               'segmented_fundi',
                                               'segmented_midaxis',
                                               'pits',
                                               'patches',
                                               'regions'],
                                output_names=['extents_sulci',
                                              'extents_fundi',
                                              'extents_midaxis',
                                              'extents_pits',
                                              'extents_patches',
                                              'extents_regions'],
                                function = measure_extents),
                    name='Measure_extents')

curvatures = node(interface = fn(input_names = ['segmented_sulci',
                                                  'segmented_fundi',
                                                  'segmented_midaxis',
                                                  'pits',
                                                  'patches',
                                                  'regions'],
                                   output_names=['curvatures_sulci',
                                                 'curvatures_fundi',
                                                 'curvatures_midaxis',
                                                 'curvatures_pits',
                                                 'curvatures_patches',
                                                 'curvatures_regions'],
                                   function = measure_curvatures),
                    name='Measure_curvatures')

depths = node(interface = fn(input_names = ['segmented_sulci',
                                              'segmented_fundi',
                                              'segmented_midaxis',
                                              'pits',
                                              'patches',
                                              'regions'],
                               output_names=['depths_sulci',
                                             'depths_fundi',
                                             'depths_midaxis',
                                             'depths_pits',
                                             'depths_patches',
                                             'depths_regions'],
                               function = measure_depths),
                 name='Measure_depths')

spectra = node(interface = fn(input_names = ['segmented_sulci',
                                               'segmented_fundi',
                                               'segmented_midaxis',
                                               'pits',
                                               'patches',
                                               'regions'],
                                output_names=['spectra_sulci',
                                              'spectra_fundi',
                                              'spectra_midaxis',
                                              'spectra_pits',
                                              'spectra_patches',
                                              'spectra_regions'],
                                function = measure_spectra),
                  name='Measure_spectra')

#-----------------------------------------------------------------------------
# Connect labeled surface patches and volumes to shape measurement nodes
#-----------------------------------------------------------------------------
featureflow.connect([(patch_extraction,  positions, [('patches', 'patches')])])
featureflow.connect([(region_extraction, positions, [('regions', 'regions')])])
featureflow.connect([(patch_extraction,  extents, [('patches', 'patches')])])
featureflow.connect([(region_extraction, extents, [('regions', 'regions')])])
featureflow.connect([(patch_extraction,  depths, [('patches', 'patches')])])
featureflow.connect([(region_extraction, depths, [('regions', 'regions')])])
featureflow.connect([(patch_extraction,  curvatures, [('patches', 'patches')])])
featureflow.connect([(region_extraction, curvatures, [('regions', 'regions')])])
featureflow.connect([(patch_extraction,  spectra, [('patches', 'patches')])])
featureflow.connect([(region_extraction, spectra, [('regions', 'regions')])])

#-----------------------------------------------------------------------------
# Connect feature to shape measurement nodes
#-----------------------------------------------------------------------------
featureflow.connect([(sulcus_segmentation, positions, [('segmented_sulci', 'segmented_sulci')])])
featureflow.connect([(fundus_segmentation, positions, [('segmented_fundi', 'segmented_fundi')])])
featureflow.connect([(pit_extraction, positions, [('pits', 'pits')])])
featureflow.connect([(midaxis_segmentation, positions, [('segmented_midaxis', 'segmented_midaxis')])])

featureflow.connect([(sulcus_segmentation, extents, [('segmented_sulci', 'segmented_sulci')])])
featureflow.connect([(fundus_segmentation, extents, [('segmented_fundi', 'segmented_fundi')])])
featureflow.connect([(midaxis_segmentation, extents, [('segmented_midaxis', 'segmented_midaxis')])])

featureflow.connect([(sulcus_segmentation, curvatures, [('segmented_sulci', 'segmented_sulci')])])
featureflow.connect([(fundus_segmentation, curvatures, [('segmented_fundi', 'segmented_fundi')])])
featureflow.connect([(pit_extraction, curvatures, [('pits', 'pits')])])
featureflow.connect([(midaxis_segmentation, curvatures, [('segmented_midaxis', 'segmented_midaxis')])])

featureflow.connect([(sulcus_segmentation, depths, [('segmented_sulci', 'segmented_sulci')])])
featureflow.connect([(fundus_segmentation, depths, [('segmented_fundi', 'segmented_fundi')])])
featureflow.connect([(pit_extraction, depths, [('pits', 'pits')])])
featureflow.connect([(midaxis_segmentation, depths, [('segmented_midaxis', 'segmented_midaxis')])])

featureflow.connect([(sulcus_segmentation, spectra, [('segmented_sulci', 'segmented_sulci')])])
featureflow.connect([(fundus_segmentation, spectra, [('segmented_fundi', 'segmented_fundi')])])
featureflow.connect([(midaxis_segmentation, spectra, [('segmented_midaxis', 'segmented_midaxis')])])

##############################################################################
#    Store surface maps, features, and measures in database
##############################################################################
#-----------------------------------------------------------------------------
# Database nodes
#-----------------------------------------------------------------------------
maps_database = node(interface = fn(input_names = ['depth_file',
                                                     'mean_curv_file',
                                                     'gauss_curv_file'],
                                      output_names=['success'],
                                      function = write_surfaces_to_database),
                        name='Write_surfaces_to_database')

features_database = node(interface = fn(input_names = ['segmented_sulci',
                                                         'segmented_fundi',
                                                         'pits',
                                                         'segmented_midaxis'],
                                          output_names=['success'],
                                          function = write_features_to_database),
                            name='Write_features_to_database')

measures_database = node(interface = fn(input_names = ['positions_sulci',
                                                         'positions_fundi',
                                                         'positions_pits',
                                                         'positions_midaxis',
                                                         'positions_patches',
                                                         'positions_regions',
                                                         'extents_sulci',
                                                         'extents_fundi',
                                                         'extents_midaxis',
                                                         'extents_patches',
                                                         'extents_regions',
                                                         'curvatures_sulci',
                                                         'curvatures_fundi',
                                                         'curvatures_pits',
                                                         'curvatures_midaxis',
                                                         'curvatures_patches',
                                                         'curvatures_regions',
                                                         'depths_sulci',
                                                         'depths_fundi',
                                                         'depths_pits',
                                                         'depths_midaxis',
                                                         'depths_patches',
                                                         'depths_regions',
                                                         'spectra_sulci',
                                                         'spectra_fundi',
                                                         'spectra_midaxis',
                                                         'spectra_patches',
                                                         'spectra_regions'],
                                          output_names=['measures'],
                                          function = write_measures_to_database),
                            name='Write_measures_to_database')

measures_table = node(interface = fn(input_names = ['measures'],
                                       output_names=['success'],
                                       function = write_measures_to_table),
                         name='Write_measures_to_table')

#-----------------------------------------------------------------------------
# Connect surface maps to database nodes
#-----------------------------------------------------------------------------
featureflow.connect([(surfaces, maps_database, [('depth_file','depth_file'),
                                ('mean_curv_file','mean_curv_file'),
                                ('gauss_curv_file','gauss_curv_file')])])

# Connect feature to database nodes
featureflow.connect([(sulcus_segmentation, features_database, [('segmented_sulci', 'segmented_sulci')]),
              (fundus_segmentation, features_database, [('segmented_fundi', 'segmented_fundi')]),
              (pit_extraction, features_database, [('pits', 'pits')]),
              (midaxis_segmentation, features_database,
                          [('segmented_midaxis', 'segmented_midaxis')])])

# Connect feature measures to database nodes
featureflow.connect([(positions, measures_database, [('positions_sulci', 'positions_sulci'),
                                              ('positions_fundi', 'positions_fundi'),
                                              ('positions_pits', 'positions_pits'),
                                              ('positions_midaxis', 'positions_midaxis')]),
              (extents, measures_database, [('extents_sulci', 'extents_sulci'),
                                            ('extents_fundi', 'extents_fundi'),
                                            ('extents_midaxis', 'extents_midaxis')]),
              (curvatures, measures_database, [('curvatures_sulci', 'curvatures_sulci'),
                                               ('curvatures_fundi', 'curvatures_fundi'),
                                               ('curvatures_pits', 'curvatures_pits'),
                                               ('curvatures_midaxis', 'curvatures_midaxis')]),
              (depths, measures_database, [('depths_sulci', 'depths_sulci'),
                                           ('depths_fundi', 'depths_fundi'),
                                           ('depths_pits', 'depths_pits'),
                                           ('depths_midaxis', 'depths_midaxis')]),
              (spectra, measures_database, [('spectra_sulci', 'spectra_sulci'),
                                            ('spectra_fundi', 'spectra_fundi'),
                                            ('spectra_midaxis', 'spectra_midaxis')])])

# Connect label measures to database nodes
featureflow.connect([(positions, measures_database, [('positions_patches', 'positions_patches'),
                                              ('positions_regions', 'positions_regions')]),
              (extents, measures_database, [('extents_patches', 'extents_patches'),
                                            ('extents_regions', 'extents_regions')]),
              (curvatures, measures_database, [('curvatures_patches', 'curvatures_patches'),
                                               ('curvatures_regions', 'curvatures_regions')]),
              (depths, measures_database, [('depths_patches', 'depths_patches'),
                                           ('depths_regions', 'depths_regions')]),
              (spectra, measures_database, [('spectra_patches', 'spectra_patches'),
                                            ('spectra_regions', 'spectra_regions')])])

# Connect measure to table nodes
featureflow.connect([(measures_database, measures_table, [('measures', 'measures')])])
"""
"""
##############################################################################
#
#   Label evaluation workflow
#
##############################################################################

evalflow = workflow(name='Evaluation_workflow')

##############################################################################
#   Surface calculations
##############################################################################

#-----------------------------------------------------------------------------
# Measure surface depth and curvature nodes
#-----------------------------------------------------------------------------
depth = node(name='Compute_depth',
             interface = fn(function = compute_depth,
                            input_names = ['command',
                                           'labels'],
                            output_names = ['overlap']))
depth.inputs.command = depth_command

mbflow.connect([(atlasflow, evalflow,
                 [('Label_vote.output_files',
                   'Evaluation_workflow.maxlabels')])])

featureflow.connect([(propagate_volume, extract_region, [('labels', 'labels')])])
featureflow.connect([(propagate, extract_patch, [('labels', 'labels')])])
"""
##############################################################################
#    Run workflow
##############################################################################
if __name__== '__main__':

    if do_create_graph:
        mbflow.write_graph(graph2use='flat')
        mbflow.write_graph(graph2use='hierarchical')
    mbflow.run(updatehash=False)  #mbflow.run(updatehash=True)



