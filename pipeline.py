#!/usr/bin/python

"""
This is Mindboggle's NiPype pipeline!

Example usage:

FIX:
from pipeline import create_mindboggle_flow
wf = create_mindboggle_flow()
wf.inputs.feature_extractor.curv_file = \
'/projects/mindboggle/data/ManualSurfandVolLabels/subjects/KKI2009-14/surf/lh.curv'
wf.inputs.feature_extractor.surface_file = \
'/projects/mindboggle/data/ManualSurfandVolLabels/subjects/KKI2009-14/surf/lh.pial'
wf.run() # doctest: +SKIP


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.utility as util     # utility
import nipype.interfaces.io as nio
import numpy as np

from atlas_based import convert_to_vtk, convert_to_vtk2, \
                        register_template, register_atlases, multilabel
from feature_based import *

use_freesurfer_surfaces = 1
use_inflated_surfaces = 1

# Subjects
subjects_list = ['KKI2009-11']
subjects_path = '/usr/local/freesurfer/subjects'

# Paths
base_directory = '/projects/mindboggle'
code_directory = os.path.join(base_directory, 'mindboggle')
templates_path = os.path.join(base_directory, 'data/templates_freesurfer')
atlases_path = subjects_path

working_directory = os.path.join(base_directory, 'results/workingdir')
#os.makedirs(working_directory)

depth_command = os.path.join(code_directory,
                             'measure/bin/travel_depth/TravelDepthMain')
curvature_command = os.path.join(code_directory,
                                 'measure/bin/curvature/CurvatureMain')
extract_fundi_command = os.path.join(code_directory, 
                                     'extract/fundi/extract.py')

##############################################################################
#
#   Mindboggle workflow combining:
#   * Multi-atlas registration-based labeling workflow
#   * Feature-based labeling and shape analysis workflow
#
##############################################################################
mbflow = pe.Workflow(name='pipeline')
mbflow.base_dir = working_directory

##############################################################################
#
#   Multi-atlas registration-based labeling workflow
#
##############################################################################
flo1 = pe.Workflow(name='atlas_based')
flo1.base_dir = working_directory

# Input and output nodes
infosource = pe.Node(interface=util.IdentityInterface(fields=['subject_id']),
                     name = 'Subjects')
if use_freesurfer_surfaces:
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['fs_surface_files']),
                         name = 'Surfaces')
    if use_inflated_surfaces:
        datasource2 = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                        outfields=['fs_inflated_files']),
                              name = 'Inflated_surfaces')
else:
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['surface_files']),
                         name = 'Surfaces')
    if use_inflated_surfaces:
        datasource2 = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                        outfields=['inflated_files']),
                              name = 'Inflated_surfaces')

# Iterate over subjects
infosource.iterables = ('subject_id', subjects_list)

# Specify the location and structure of the inputs and outputs
datasource.inputs.base_directory = subjects_path
datasource.inputs.template = '%s/surf/%s.%s'
if use_inflated_surfaces:
    datasource2.inputs.base_directory = subjects_path
    datasource2.inputs.template = '%s/surf/%s.%s'
if use_freesurfer_surfaces:
    datasources = [['subject_id', ['lh','rh'], 'pial']]
    datasource.inputs.template_args['fs_surface_files'] = datasources
    if use_inflated_surfaces:
        datasources2 = [['subject_id', ['lh','rh'], 'inflated']]
        datasource2.inputs.template_args['fs_inflated_files'] = datasources2
else:
    datasources = [['subject_id', ['lh','rh'], 'pial.vtk']]
    datasource.inputs.template_args['surface_files'] = datasources
    if use_inflated_surfaces:
        datasources2 = [['subject_id', ['lh','rh'], 'inflated.vtk']]   
        datasource2.inputs.template_args['inflated_files'] = datasources2

# Define a template node
template = pe.Node(interface=util.IdentityInterface(fields=[
                   'template_name', 'templates_path', 'reg_name']),
                   name='Template')
template_id = 'KKI'
template.inputs.template_name = template_id + '_2.tif'
template.inputs.templates_path = templates_path
reg_name = 'sphere_to_' + template_id + '_template.reg'
template.inputs.reg_name = reg_name

# Define an atlases node
atlases = pe.Node(interface=util.IdentityInterface(fields=[
                  'atlases_path', 'atlas_list_file', 'annot_name']),
                  name='Atlases')
atlases.inputs.atlases_path = atlases_path
atlases.inputs.atlas_list_file = os.path.join(atlases_path, 'MMRR.txt')
annot_name = 'aparcNMMjt.annot'
atlases.inputs.annot_name = annot_name

##############################################################################
#   Surface input and conversion
##############################################################################

# Connect input nodes
flo1.connect([(infosource, datasource, [('subject_id','subject_id')])])
if use_inflated_surfaces:
    flo1.connect([(infosource, datasource2, [('subject_id','subject_id')])])

# Convert FreeSurfer surfaces to VTK format
if use_freesurfer_surfaces:
    surface_conversion = pe.Node(util.Function(input_names = ['fs_surface_files'],
                                               output_names = ['surface_files'],
                                               function = convert_to_vtk),
                                 name='Convert_surfaces')
    # Connect input to surface node
    flo1.connect([(datasource, surface_conversion,
                   [('fs_surface_files','fs_surface_files')])])

    if use_inflated_surfaces:
        surface_conversion2 = pe.Node(util.Function(input_names = ['fs_inflated_files'],
                                                    output_names = ['inflated_files'],
                                                    function = convert_to_vtk2),
                                      name='Convert_inflated_surfaces')
        flo1.connect([(datasource2, surface_conversion2,
                       [('fs_inflated_files','fs_inflated_files')])])

##############################################################################
#   Multi-atlas registration
##############################################################################

# Registration nodes
"""
Example: ['/Applications/freesurfer/subjects/bert/surf/lh.sphere',
          '/Applications/freesurfer/subjects/bert/surf/rh.sphere']
         ['./templates_freesurfer/lh.KKI_2.tif',
          './templates_freesurfer/rh.KKI_2.tif']
"""
template_registration = pe.Node(util.Function(input_names=['subject_id',
                                              'subjects_path',
                                              'template_name', 
                                              'templates_path', 
                                              'reg_name'],
                                     output_names=['reg_name'],
                                     function = register_template),
                                name='Register_template')
template_registration.inputs.subjects_path = subjects_path

atlas_registration = pe.Node(util.Function(input_names=['subject_id',
                                                        'subjects_path',
                                                        'atlas_list_file',
                                                        'annot_name',
                                                        'reg_name'],
                                           output_names=['annot_name'],
                                           function = register_atlases),
                             name='Register_atlases')
atlas_registration.inputs.subjects_path = subjects_path

# Output majority vote rule labels
majority_vote = pe.Node(util.Function(input_names=['subject_id',
                                                   'subjects_path',
                                                   'annot_name',
                                                   'use_inflated_surfaces'],
                                      output_names=['annot_name'],
                                      function = multilabel),
                        name='Vote_majority')
majority_vote.inputs.subjects_path = subjects_path
majority_vote.inputs.use_inflated_surfaces = use_inflated_surfaces

# Connect input to registration and labeling nodes
flo1.connect([(infosource, template_registration, 
               [('subject_id', 'subject_id')])])
flo1.connect([(infosource, atlas_registration, 
               [('subject_id', 'subject_id')])])
flo1.connect([(infosource, majority_vote,
               [('subject_id', 'subject_id')])])

# Connect template and atlases to registration and labeling nodes
flo1.connect([(template, template_registration, 
               [('template_name', 'template_name'),
                ('templates_path', 'templates_path'),
                ('reg_name', 'reg_name')])])
flo1.connect([(atlases, atlas_registration, 
               [('atlas_list_file', 'atlas_list_file'),
                ('annot_name', 'annot_name')])])

# Connect template registration to labeling nodes
flo1.connect([(template_registration, atlas_registration, 
               [('reg_name', 'reg_name')])])
flo1.connect([(atlas_registration, majority_vote, 
               [('annot_name', 'annot_name')])])


##############################################################################
#
#   Feature-based labeling and shape analysis workflow
#
##############################################################################
flo2 = pe.Workflow(name='feature_based')
flo2.base_dir = working_directory

##############################################################################
#   Surface calculations
##############################################################################

# Measure surface depth and curvature nodes
surface_depth = pe.Node(util.Function(input_names = ['depth_command',
                                                     'surface_files'],
                                 output_names = ['surface_files',
                                                 'depth_files'],
                                 function = measure_surface_depth),
                   name='Measure_surface_depth')
surface_depth.inputs.depth_command = depth_command

surface_curvature = pe.Node(util.Function(input_names = ['curvature_command',
                                                     'surface_files'],
                                 output_names = ['surface_files',
                                                 'curvature_files'],
                                 function = measure_surface_curvature),
                   name='Measure_surface_curvature')
surface_curvature.inputs.curvature_command = curvature_command

# Connect input to surface depth and curvature nodes
flo2.connect([(infosource, datasource, [('subject_id','subject_id')])])

if use_freesurfer_surfaces:
    flo2.connect([(datasource, surface_conversion,
                   [('fs_surface_files','fs_surface_files')])])
    flo2.connect([(surface_conversion, surface_depth, 
                   [('surface_files','surface_files')])])
    flo2.connect([(surface_conversion, surface_curvature, 
                   [('surface_files','surface_files')])])
#    mbflow.connect([(flo1, flo2, 
#                   [('surface_conversion.surface_files','surface_depth.surface_files')])])
else:
    # Connect input to surface surfaces node
    flo2.connect([(datasource, surfaces,
                   [('surface_files','surface_files')])])

##############################################################################
#   Feature extraction
##############################################################################

# Feature extraction nodes
fundus_extraction = pe.Node(util.Function(input_names = ['extract_fundi_command',
                                                         'surface_files', 
                                                         'depth_curv_files'],
                                          output_names = ['fundi'],
                                          function = extract_fundi),
                            name='Extract_fundi')
fundus_extraction.inputs.extract_fundi_command = extract_fundi_command

sulcus_extraction = pe.Node(util.Function(input_names = ['depth_file',
                                                         'mean_curv_file',
                                                         'gauss_curv_file'],
                                          output_names = ['sulci'],
                                          function = extract_sulci),
                            name='Extract_sulci')

midaxis_extraction = pe.Node(util.Function(input_names = ['depth_file',
                                                          'mean_curv_file',
                                                          'gauss_curv_file'],
                                           output_names = ['midaxis'],
                                           function = extract_midaxis),
                             name='Extract_midaxis')

# Connect surface surfaces node to feature extraction nodes
#flo2.connect([(surfaces, fundus_extraction, 
#               [('surface_files', 'surface_files'),
#                ('depth_files', 'depth_files')])])
"""
flo2.connect([(surfaces, sulcus_extraction, 
               [('depth_file', 'depth_file'),
                ('mean_curv_file', 'mean_curv_file'),
                ('gauss_curv_file', 'gauss_curv_file')])])
flo2.connect([(surfaces, midaxis_extraction, 
               [('depth_file', 'depth_file'),
                ('mean_curv_file', 'mean_curv_file'),
                ('gauss_curv_file', 'gauss_curv_file')])])

# Save output
#flo2.connect([(surfaces, datasink, [('depth_file', 'depth_file'),
#                                ('mean_curv_file', 'mean_curv_file'),
#                                ('gauss_curv_file', 'gauss_curv_file')])])
"""
"""
##############################################################################
#   Label propagation
##############################################################################

# Label propagation node
label_propagation = pe.Node(util.Function(input_names=['labels', 'fundi'],
                                          output_names=['labels'],
                                          function = propagate_labels),
                            name='Propagate_labels')

# Volume label propagation node
volume_propagation = pe.Node(util.Function(input_names=['labels'],
                                           output_names=['labels'],
                                           function = propagate_volume_labels),
                             name='Propagate_volume_labels')

# Labeled surface patch and volume extraction nodes
patch_extraction = pe.Node(util.Function(input_names=['labels'],
                                         output_names=['patches'],
                                         function = extract_patches),
                           name='Extract_patches')

region_extraction = pe.Node(util.Function(input_names=['labels'],
                                          output_names=['regions'],
                                          function = extract_regions),
                            name='Extract_regions')

# Connect multiatlas registration(-based labeling) to label propagation nodes
flo2.connect([(atlas_registration, label_propagation, [('labels','labels')]),
              (fundus_extraction, label_propagation, [('fundi','fundi')])])

# Connect label propagation to labeled surface patch and volume extraction nodes
flo2.connect([(label_propagation, volume_propagation, [('labels', 'labels')])])
flo2.connect([(volume_propagation, region_extraction, [('labels', 'labels')])])
flo2.connect([(label_propagation, patch_extraction, [('labels', 'labels')])])

##############################################################################
#   Feature segmentation / identification
##############################################################################

# Feature segmentation nodes
sulcus_segmentation = pe.Node(util.Function(input_names=['sulci','labels'],
                                            output_names=['segmented_sulci'],
                                            function = segment_sulci),
                              name='Segment_sulci')

fundus_segmentation = pe.Node(util.Function(input_names=['fundi','labels'],
                                            output_names=['segmented_fundi'],
                                            function = segment_fundi),
                              name='Segment_fundi')

midaxis_segmentation = pe.Node(util.Function(input_names=['midaxis','labels'],
                                             output_names=['segmented_midaxis'],
                                             function = segment_midaxis),
                               name='Segment_midaxis')

# Connect feature and feature segmentation nodes
flo2.connect([(sulcus_extraction, sulcus_segmentation, [('sulci','sulci')]),
              (fundus_extraction, fundus_segmentation, [('fundi','fundi')]),
              (midaxis_extraction, midaxis_segmentation, [('midaxis','midaxis')])])

# Connect multiatlas registration(-based labeling) and feature segmentation nodes
flo2.connect([(label_propagation, sulcus_segmentation, [('labels','labels')]),
              (label_propagation, fundus_segmentation, [('labels','labels')]),
              (sulcus_segmentation, midaxis_segmentation, [('segmented_sulci','labels')])])
              
##############################################################################
#   Shape measurement
##############################################################################

# Shape measurement nodes
positions = pe.Node(util.Function(input_names = ['segmented_sulci', 
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

extents = pe.Node(util.Function(input_names = ['segmented_sulci', 
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

curvatures = pe.Node(util.Function(input_names = ['segmented_sulci', 
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

depths = pe.Node(util.Function(input_names = ['segmented_sulci',
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

spectra = pe.Node(util.Function(input_names = ['segmented_sulci',
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
flo2.connect([(patch_extraction,  positions, [('patches', 'patches')])])
flo2.connect([(region_extraction, positions, [('regions', 'regions')])])
flo2.connect([(patch_extraction,  extents, [('patches', 'patches')])])
flo2.connect([(region_extraction, extents, [('regions', 'regions')])])
flo2.connect([(patch_extraction,  depths, [('patches', 'patches')])])
flo2.connect([(region_extraction, depths, [('regions', 'regions')])])
flo2.connect([(patch_extraction,  curvatures, [('patches', 'patches')])])
flo2.connect([(region_extraction, curvatures, [('regions', 'regions')])])
flo2.connect([(patch_extraction,  spectra, [('patches', 'patches')])])
flo2.connect([(region_extraction, spectra, [('regions', 'regions')])])

# Connect feature to shape measurement nodes
flo2.connect([(sulcus_segmentation, positions, [('segmented_sulci', 'segmented_sulci')])])
flo2.connect([(fundus_segmentation, positions, [('segmented_fundi', 'segmented_fundi')])])
flo2.connect([(pit_extraction, positions, [('pits', 'pits')])])
flo2.connect([(midaxis_segmentation, positions, [('segmented_midaxis', 'segmented_midaxis')])])

flo2.connect([(sulcus_segmentation, extents, [('segmented_sulci', 'segmented_sulci')])])
flo2.connect([(fundus_segmentation, extents, [('segmented_fundi', 'segmented_fundi')])])
flo2.connect([(midaxis_segmentation, extents, [('segmented_midaxis', 'segmented_midaxis')])])

flo2.connect([(sulcus_segmentation, curvatures, [('segmented_sulci', 'segmented_sulci')])])
flo2.connect([(fundus_segmentation, curvatures, [('segmented_fundi', 'segmented_fundi')])])
flo2.connect([(pit_extraction, curvatures, [('pits', 'pits')])])
flo2.connect([(midaxis_segmentation, curvatures, [('segmented_midaxis', 'segmented_midaxis')])])

flo2.connect([(sulcus_segmentation, depths, [('segmented_sulci', 'segmented_sulci')])])
flo2.connect([(fundus_segmentation, depths, [('segmented_fundi', 'segmented_fundi')])])
flo2.connect([(pit_extraction, depths, [('pits', 'pits')])])
flo2.connect([(midaxis_segmentation, depths, [('segmented_midaxis', 'segmented_midaxis')])])

flo2.connect([(sulcus_segmentation, spectra, [('segmented_sulci', 'segmented_sulci')])])
flo2.connect([(fundus_segmentation, spectra, [('segmented_fundi', 'segmented_fundi')])])
flo2.connect([(midaxis_segmentation, spectra, [('segmented_midaxis', 'segmented_midaxis')])])

##############################################################################
#    Store surface maps, features, and measures in database
##############################################################################

# Database nodes
maps_database = pe.Node(util.Function(input_names = ['depth_file',
                                                     'mean_curv_file',
                                                     'gauss_curv_file'],
                                      output_names=['success'],
                                      function = write_surfaces_to_database),
                        name='Write_surfaces_to_database')

features_database = pe.Node(util.Function(input_names = ['segmented_sulci',
                                                         'segmented_fundi',
                                                         'pits',
                                                         'segmented_midaxis'],
                                          output_names=['success'],
                                          function = write_features_to_database),
                            name='Write_features_to_database')

measures_database = pe.Node(util.Function(input_names = ['positions_sulci',
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

measures_table = pe.Node(util.Function(input_names = ['measures'],
                                       output_names=['success'],
                                       function = write_measures_to_table),
                         name='Write_measures_to_table')

# Connect surface maps to database nodes
flo2.connect([(surfaces, maps_database, [('depth_file','depth_file'),
                                ('mean_curv_file','mean_curv_file'),
                                ('gauss_curv_file','gauss_curv_file')])])

# Connect feature to database nodes
flo2.connect([(sulcus_segmentation, features_database, [('segmented_sulci', 'segmented_sulci')]),
              (fundus_segmentation, features_database, [('segmented_fundi', 'segmented_fundi')]),
              (pit_extraction, features_database, [('pits', 'pits')]),
              (midaxis_segmentation, features_database, 
                          [('segmented_midaxis', 'segmented_midaxis')])])

# Connect feature measures to database nodes
flo2.connect([(positions, measures_database, [('positions_sulci', 'positions_sulci'),
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
flo2.connect([(positions, measures_database, [('positions_patches', 'positions_patches'),
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
flo2.connect([(measures_database, measures_table, [('measures', 'measures')])])

"""

##############################################################################
#    Run workflow
##############################################################################
if __name__== '__main__':

    flo1.run()
    #flo2.run()
"""
    mbflow.connect([(flo1, flo2,
                     [('', '')])])
    flo2.connect([(surface_conversion, surfaces, 
                   [('surface_files','surface_files')])])

    mbflow.write_graph(graph2use='flat')
    mbflow.write_graph(graph2use='hierarchical')

    mbflow.run()
"""
