#!/usr/bin/python

"""
This is Mindboggle's NiPype pipeline!

Example usage:

wf.run() # doctest: +SKIP


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.utility as util     # utility
import nipype.interfaces.io as nio
import numpy as np

from atlases import register_to_template, transform_atlas_labels,\
                    majority_vote_label
from features import *

use_freesurfer_surfaces = 1
hemis = ['lh','rh']
surface_types = ['pial'] #,'inflated']
template_id = 'KKI'
template_name = template_id + '_2.tif'
template_transform = 'sphere_to_' + template_id + '_template.reg'
atlas_annot_name = 'aparcNMMjt.annot'

# Subjects
subjects_list = ['KKI2009-11'] #, 'KKI2009-14']

use_linux_paths = 0
if use_linux_paths:
    subjects_path = '/usr/local/freesurfer/subjects'
else:
    subjects_path = '/Applications/freesurfer/subjects'

# Paths
templates_path = '/projects/mindboggle/data/templates_freesurfer'
atlases_path = subjects_path

# Output directory
results_path = '/projects/mindboggle/results/'
working_path = results_path + 'workingdir'
if not os.path.isdir(results_path):
    os.makedirs(results_path)
if not os.path.isdir(working_path):
    os.makedirs(working_path)

# Commands
depth_command = './measure/surface_measures/bin/travel_depth/TravelDepthMain'
curvature_command = './measure/surface_measures/bin/curvature/CurvatureMain'
extract_fundi_command = './extract/fundi/vtk_extract.py'

# List of atlas subjects
atlas_list_file = os.path.join(atlases_path, 'MMRR.txt')
f = open(atlas_list_file)
atlas_list_lines = f.readlines()
atlas_names = [a.strip("\n") for a in atlas_list_lines if a.strip("\n")]

##############################################################################
#
#   Mindboggle workflow combining:
#   * Multi-atlas registration-based labeling workflow
#   * Feature-based labeling and shape analysis workflow
#   * Analytics
#
##############################################################################
mbflow = pe.Workflow(name='Mindboggle_workflow')
mbflow.base_dir = working_path

# Iterate inputs over subjects, hemispheres, surface types
infosource = pe.Node(name = 'Inputs',
                     interface = util.IdentityInterface(fields=['subject_id',
                                                                'hemi',
                                                                'surface_type']))
infosource.iterables = ([('subject_id', subjects_list),
                         ('hemi', hemis)])
datasource = pe.Node(name = 'Surfaces',
                     interface = nio.DataGrabber(infields=['subject_id',
                                                           'hemi'],
                                                 outfields=['surface_files',
                                                            'inf_surface_files',
                                                            'sph_surface_files']))

# Specify the location and structure of the inputs and outputs
datasource.inputs.base_directory = subjects_path
datasource.inputs.template = '%s/surf/%s.%s'
datasource.inputs.template_args['surface_files'] = [['subject_id', 
                                                     'hemi', 
                                                     'pial']]
datasource.inputs.template_args['inf_surface_files'] = [['subject_id', 
                                                         'hemi', 
                                                         'inflated']]
datasource.inputs.template_args['sph_surface_files'] = [['subject_id', 
                                                         'hemi', 
                                                         'sphere']]
datasink = pe.Node(nio.DataSink(), name = 'Results')
datasink.inputs.base_directory = results_path
datasink.inputs.container = 'output'

# Connect input nodes
mbflow.connect([(infosource, datasource, 
                 [('subject_id','subject_id'),
                  ('hemi','hemi')])])

##############################################################################
#   Surface input and conversion
##############################################################################

# Convert FreeSurfer surfaces to VTK format
if use_freesurfer_surfaces:

    import nipype.interfaces.freesurfer as fs

    convert = pe.MapNode(name = 'Convert_surface',
                         iterfield=['in_file'],
                         interface = fs.MRIsConvert(out_datatype='vtk'))
    mbflow.connect([(datasource, convert, 
                     [('surface_files','in_file')])])

##############################################################################
#
#   Multi-atlas registration-based labeling workflow
#
##############################################################################
atlasflow = pe.Workflow(name='Atlas_workflow')

##############################################################################
#   Multi-atlas registration
##############################################################################

# Template registration
register = pe.Node(name = 'Register_to_template',
                   interface = util.Function(
                                    function = register_to_template,
                                    input_names = ['hemi',
                                                   'sph_surface_file',
                                                   'template_name',
                                                   'templates_path',
                                                   'template_transform'],
                                    output_names = ['template_transform']))
register.inputs.template_name = template_name
register.inputs.templates_path = templates_path
register.inputs.template_transform = template_transform

# Atlas registration
transform = pe.MapNode(name = 'Transform_atlas_labels',
                       iterfield = ['atlas_name'],
                       interface = util.Function(
                                        function = transform_atlas_labels,
                                        input_names = ['hemi',
                                                       'subject_id',
                                                       'template_transform',
                                                       'atlas_name',
                                                       'atlases_path',
                                                       'atlas_annot_name'],
                                        output_names = ['output_file']))
transform.inputs.atlas_name = atlas_names
transform.inputs.atlases_path = atlases_path
transform.inputs.atlas_annot_name = atlas_annot_name

# Majority vote labeling
vote = pe.Node(name='Majority_vote',
               interface = util.Function(
                                function = majority_vote_label,
                                input_names = ['surface_file',
                                               'annot_files'],
                                output_names = ['output_files']))

# Add and connect the above nodes
atlasflow.add_nodes([register, transform, vote])

mbflow.connect([(infosource, atlasflow,
                 [('hemi', 'Register_to_template.hemi')]),
                (datasource, atlasflow, 
                 [('sph_surface_files',
                   'Register_to_template.sph_surface_file')])])

mbflow.connect([(infosource, atlasflow, 
                 [('hemi', 'Transform_atlas_labels.hemi'),
                  ('subject_id', 'Transform_atlas_labels.subject_id')])])
atlasflow.connect([(register, transform, 
                    [('template_transform', 'template_transform')])])

if use_freesurfer_surfaces:
    mbflow.connect([(convert, atlasflow, 
                     [('converted', 'Majority_vote.surface_file')])])
else:
    mbflow.connect([(datasource, atlasflow, 
                     [('surface_files', 'Majority_vote.surface_file')])])
atlasflow.connect([(transform, vote,
                    [('output_file', 'annot_files')])])
#mbflow.connect([(atlasflow, datasink,
#                 [('Majority_vote.output_files', 'maxlabels')])])

##############################################################################
#
#   Feature-based labeling and shape analysis workflow
#
##############################################################################

featureflow = pe.Workflow(name='Feature_workflow')

##############################################################################
#   Surface calculations
##############################################################################

# Measure surface depth and curvature nodes
depth = pe.Node(name='Compute_depth',
                interface = util.Function(
                                 function = compute_depth,
                                 input_names = ['command',
                                                'surface_file'],
                                 output_names = ['depth_file']))
depth.inputs.command = depth_command

curvature = pe.Node(name='Compute_curvature',
                    interface = util.Function(
                                     function = compute_curvature),
                                     input_names = ['command',
                                                    'surface_file'],
                                     output_names = ['mean_curvature_file',
                                                     'gauss_curvature_file',
                                                     'max_curvature_file',
                                                     'min_curvature_file'])
curvature.inputs.command = curvature_command

# Add and connect nodes
featureflow.add_nodes([depth, curvature])

if use_freesurfer_surfaces:
    mbflow.connect([(convert, featureflow, 
                   [('converted', 'Compute_depth.surface_file')])])
    mbflow.connect([(convert, featureflow, 
                   [('converted', 'Compute_curvature.surface_file')])])
else:
    # Connect input to surface depth and curvature nodes
    mbflow.connect([(atlasflow, featureflow, 
                     [('Surfaces.surface_files',
                       'Compute_depth.surface_file')])])
    mbflow.connect([(atlasflow, featureflow, 
                     [('Surfaces.surface_files',
                       'Compute_curvature.surface_file')])])

# Save
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
                   'surfaces.@min_curvature')])])

##############################################################################
#   Feature extraction
##############################################################################

# Extract features
fundi = pe.Node(name='Extract_fundi',
                interface = util.Function(
                                 function = extract_fundi
                                 input_names = ['command',
                                                'depth_file'],
                                 output_names = ['fundi']))
fundi.inputs.command = extract_fundi_command

"""
sulcus_extraction = pe.Node(interface = util.Function(input_names = ['depth_file',
                                                         'mean_curv_file',
                                                         'gauss_curv_file'],
                                          output_names = ['sulci'],
                                          function = extract_sulci),
                            name='Extract_sulci')

midaxis_extraction = pe.Node(interface = util.Function(input_names = ['depth_file',
                                                          'mean_curv_file',
                                                          'gauss_curv_file'],
                                           output_names = ['midaxis'],
                                           function = extract_midaxis),
                             name='Extract_midaxis')

"""
# Connect surface depth to feature extraction nodes
featureflow.connect([(surface_depth, fundus_extraction, 
               [('depth_file', 'depth_file')])])
featureflow.connect([(surface_depth, datasink, 
               [('depth_file', 'surface_depth')])])
"""
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

# Label propagation node
label_propagation = pe.Node(interface = util.Function(input_names=['labels', 'fundi'],
                                          output_names=['labels'],
                                          function = propagate_labels),
                            name='Propagate_labels')

# Volume label propagation node
volume_propagation = pe.Node(interface = util.Function(input_names=['labels'],
                                           output_names=['labels'],
                                           function = propagate_volume_labels),
                             name='Propagate_volume_labels')

# Labeled surface patch and volume extraction nodes
patch_extraction = pe.Node(interface = util.Function(input_names=['labels'],
                                         output_names=['patches'],
                                         function = extract_patches),
                           name='Extract_patches')

region_extraction = pe.Node(interface = util.Function(input_names=['labels'],
                                          output_names=['regions'],
                                          function = extract_regions),
                            name='Extract_regions')

# Connect multiatlas registration(-based labeling) to label propagation nodes
featureflow.connect([(transform, label_propagation, [('labels','labels')]),
              (fundus_extraction, label_propagation, [('fundi','fundi')])])

# Connect label propagation to labeled surface patch and volume extraction nodes
featureflow.connect([(label_propagation, volume_propagation, [('labels', 'labels')])])
featureflow.connect([(volume_propagation, region_extraction, [('labels', 'labels')])])
featureflow.connect([(label_propagation, patch_extraction, [('labels', 'labels')])])

##############################################################################
#   Feature segmentation / identification
##############################################################################

# Feature segmentation nodes
sulcus_segmentation = pe.Node(interface = util.Function(input_names=['sulci','labels'],
                                            output_names=['segmented_sulci'],
                                            function = segment_sulci),
                              name='Segment_sulci')

fundus_segmentation = pe.Node(interface = util.Function(input_names=['fundi','labels'],
                                            output_names=['segmented_fundi'],
                                            function = segment_fundi),
                              name='Segment_fundi')

midaxis_segmentation = pe.Node(interface = util.Function(input_names=['midaxis','labels'],
                                             output_names=['segmented_midaxis'],
                                             function = segment_midaxis),
                               name='Segment_midaxis')

# Connect feature and feature segmentation nodes
featureflow.connect([(sulcus_extraction, sulcus_segmentation, [('sulci','sulci')]),
              (fundus_extraction, fundus_segmentation, [('fundi','fundi')]),
              (midaxis_extraction, midaxis_segmentation, [('midaxis','midaxis')])])

# Connect multiatlas registration(-based labeling) and feature segmentation nodes
featureflow.connect([(label_propagation, sulcus_segmentation, [('labels','labels')]),
              (label_propagation, fundus_segmentation, [('labels','labels')]),
              (sulcus_segmentation, midaxis_segmentation, [('segmented_sulci','labels')])])
              
##############################################################################
#   Shape measurement
##############################################################################

# Shape measurement nodes
positions = pe.Node(interface = util.Function(input_names = ['segmented_sulci', 
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

extents = pe.Node(interface = util.Function(input_names = ['segmented_sulci', 
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

curvatures = pe.Node(interface = util.Function(input_names = ['segmented_sulci', 
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

depths = pe.Node(interface = util.Function(input_names = ['segmented_sulci',
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

spectra = pe.Node(interface = util.Function(input_names = ['segmented_sulci',
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

# Connect labeled surface patches and volumes to shape measurement nodes
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

# Connect feature to shape measurement nodes
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

# Database nodes
maps_database = pe.Node(interface = util.Function(input_names = ['depth_file',
                                                     'mean_curv_file',
                                                     'gauss_curv_file'],
                                      output_names=['success'],
                                      function = write_surfaces_to_database),
                        name='Write_surfaces_to_database')

features_database = pe.Node(interface = util.Function(input_names = ['segmented_sulci',
                                                         'segmented_fundi',
                                                         'pits',
                                                         'segmented_midaxis'],
                                          output_names=['success'],
                                          function = write_features_to_database),
                            name='Write_features_to_database')

measures_database = pe.Node(interface = util.Function(input_names = ['positions_sulci',
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

measures_table = pe.Node(interface = util.Function(input_names = ['measures'],
                                       output_names=['success'],
                                       function = write_measures_to_table),
                         name='Write_measures_to_table')

# Connect surface maps to database nodes
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

##############################################################################
#    Run workflow
##############################################################################
if __name__== '__main__':

    mbflow.write_graph(graph2use='flat')
    mbflow.write_graph(graph2use='hierarchical')
    #atlasflow.write_graph(graph2use='flat')
    #atlasflow.write_graph(graph2use='hierarchical')
    mbflow.run(updatehash=False)  #mbflow.run(updatehash=True)

