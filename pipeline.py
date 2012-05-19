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

import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.utility as util     # utility
import nipype.interfaces.io as nio
import numpy as np
from pipeline_functions import *

use_freesurfer_surfaces = 1

##############################################################################
#   Workflow setup
##############################################################################
flow = pe.Workflow(name='pipeline')
flow.base_dir = '.'

# Input and output nodes
infosource = pe.Node(interface=util.IdentityInterface(fields=['subject_id']),
                     name = 'Input')
if use_freesurfer_surfaces:
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['fs_surface_files']),
                         name = 'DataSource')
else:
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['surface_files']),
                         name = 'DataSource')
templates = pe.Node(interface=nio.DataGrabber(infields=['template_id'],
                                             outfields=['template_files']),
                   name = 'Template')
atlases = pe.Node(interface=nio.DataGrabber(infields=['atlas_list_file'],
                                            outfields=['atlas_list_file']),
                  name = 'Atlases')
datasink = pe.Node(interface=nio.DataSink(),
                   name = 'DataSink')

# Iterate over subjects
infosource.iterables = ('subject_id', ['KKI2009-11'])

# Specify the location and structure of the inputs and outputs
datasource.inputs.base_directory = '/projects/mindboggle/data'
datasource.inputs.template = '%s/surf/%s.%s'
if use_freesurfer_surfaces:
    datasources = [['subject_id', ['lh','rh'], 'pial']]
    datasource.inputs.template_args['fs_surface_files'] = datasources
else:
    datasources = [['subject_id', ['lh','rh'], 'pial.vtk']]
    datasource.inputs.template_args['surface_files'] = datasources

datasink.inputs.base_directory = '/projects/mindboggle'
datasink.inputs.container = 'output'

# Specify the location and structure of the templates
templates.inputs.base_directory = '../data/templates_freesurfer'
templates.inputs.template = '%s.%s_2.tif'
templates.inputs.template_args['template_files'] = [[['lh','rh'], 'template_id']]

# Specify the location and structure of the atlases
atlases.inputs.base_directory = '../data/atlases'
atlases.inputs.template = '%s/%s'
atlases.inputs.template_args['atlas_list_file'] = [['MMRR_labeled', 'KKI.txt']]

##############################################################################
#   Surface input and conversion
##############################################################################

# Connect input nodes
flow.connect([(infosource, datasource, [('subject_id','subject_id')])])

# Convert FreeSurfer surfaces to VTK format
if use_freesurfer_surfaces:
    surface_conversion = pe.Node(util.Function(input_names = ['fs_surface_files'],
                                               output_names = ['surface_files'],
                                               function = convert_to_vtk),
                                 name='Convert_surfaces')

    # Connect input to surface surface_maps node
    flow.connect([(datasource, surface_conversion,
                   [('fs_surface_files','fs_surface_files')])])

##############################################################################
#   Surface map calculation
##############################################################################

# Measure surface surface_maps node
surface_maps = pe.Node(util.Function(input_names = ['surface_files'],
                                     output_names = ['depth_curv_map_files'],
                                     function = measure_surface_maps),
                       name='Measure_surface_maps')

# Connect input to surface maps nodes
if use_freesurfer_surfaces:
    flow.connect([(surface_conversion, surface_maps, 
                   [('surface_files','surface_files')])])
else:
    # Connect input to surface surface_maps node
    flow.connect([(datasource, surface_maps,
                   [('surface_files','surface_files')])])

##############################################################################
#   Feature extraction
##############################################################################

# Feature extraction nodes
fundus_extraction = pe.Node(util.Function(input_names = ['depth_curv_map_files'],
                                          output_names = ['fundi'],
                                          function = extract_fundi),
                            name='Extract_fundi')

"""
sulcus_extraction = pe.Node(util.Function(input_names = ['depth_map',
                                                         'mean_curv_map',
                                                         'gauss_curv_map'],
                                          output_names = ['sulci'],
                                          function = extract_sulci),
                            name='Extract_sulci')

midaxis_extraction = pe.Node(util.Function(input_names = ['depth_map',
                                                             'mean_curv_map',
                                                             'gauss_curv_map'],
                                              output_names = ['midaxis'],
                                              function = extract_midaxis),
                                name='Extract_midaxis')
"""

# Connect surface surface_maps node to feature extraction nodes
#flow.connect([(surface_maps, fundus_extraction, 
#               [('depth_map', 'depth_map'),
#                ('mean_curv_map', 'mean_curv_map'),
#                ('gauss_curv_map', 'gauss_curv_map')])])

"""
flow.connect([(surface_maps, sulcus_extraction, 
               [('depth_map', 'depth_map'),
                ('mean_curv_map', 'mean_curv_map'),
                ('gauss_curv_map', 'gauss_curv_map')])])
flow.connect([(surface_maps, midaxis_extraction, 
               [('depth_map', 'depth_map'),
                ('mean_curv_map', 'mean_curv_map'),
                ('gauss_curv_map', 'gauss_curv_map')])])

# Save output
#flow.connect([(surface_maps, datasink, [('depth_map', 'depth_map'),
#                                ('mean_curv_map', 'mean_curv_map'),
#                                ('gauss_curv_map', 'gauss_curv_map')])])

##############################################################################
#   Multi-atlas registration
##############################################################################

# Registration nodes
"""
"""
Example: ['/Applications/freesurfer/subjects/bert/surf/lh.sphere',
          '/Applications/freesurfer/subjects/bert/surf/rh.sphere']
         ['./templates_freesurfer/lh.KKI_2.tif',
          './templates_freesurfer/rh.KKI_2.tif']
"""
"""
template_registration = pe.Node(util.Function(input_names=['subject_id',
                                              'subjects_path',
                                              'template_id', 
                                              'template_path', 
                                              'registration_name'],
                                     output_names=['registration_name'],
                                     function = register_template),
                                name='Register_template')

atlas_registration = pe.Node(util.Function(input_names=['subject_id',
                                                        'atlas_list_file',
                                                        'registration_name',
                                                        'output_path'],
                                           output_names=['labels'],
                                           function = register_atlases),
                             name='Register_atlases')

# Connect template and atlases to registration nodes
flow.connect([(template, template_registration, 
               [('template_id','template_id')])])
flow.connect([(atlases, atlas_registration, 
               [('atlas_list_file','atlas_list_file')])])

# Connect input to registration nodes
flow.connect([(datasource, template_registration, 
               [('subject_id','subject_id')])])

# Connect template registration to multiatlas registration-based labeling nodes
flow.connect([(template_registration, atlas_registration, 
               [('registration_name','registration_name')])])

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
flow.connect([(atlas_registration, label_propagation, [('labels','labels')]),
              (fundus_extraction, label_propagation, [('fundi','fundi')])])

# Connect label propagation to labeled surface patch and volume extraction nodes
flow.connect([(label_propagation, volume_propagation, [('labels', 'labels')])])
flow.connect([(volume_propagation, region_extraction, [('labels', 'labels')])])
flow.connect([(label_propagation, patch_extraction, [('labels', 'labels')])])

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
flow.connect([(sulcus_extraction, sulcus_segmentation, [('sulci','sulci')]),
              (fundus_extraction, fundus_segmentation, [('fundi','fundi')]),
              (midaxis_extraction, midaxis_segmentation, [('midaxis','midaxis')])])

# Connect multiatlas registration(-based labeling) and feature segmentation nodes
flow.connect([(label_propagation, sulcus_segmentation, [('labels','labels')]),
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
flow.connect([(patch_extraction,  positions, [('patches', 'patches')])])
flow.connect([(region_extraction, positions, [('regions', 'regions')])])
flow.connect([(patch_extraction,  extents, [('patches', 'patches')])])
flow.connect([(region_extraction, extents, [('regions', 'regions')])])
flow.connect([(patch_extraction,  depths, [('patches', 'patches')])])
flow.connect([(region_extraction, depths, [('regions', 'regions')])])
flow.connect([(patch_extraction,  curvatures, [('patches', 'patches')])])
flow.connect([(region_extraction, curvatures, [('regions', 'regions')])])
flow.connect([(patch_extraction,  spectra, [('patches', 'patches')])])
flow.connect([(region_extraction, spectra, [('regions', 'regions')])])

# Connect feature to shape measurement nodes
flow.connect([(sulcus_segmentation, positions, [('segmented_sulci', 'segmented_sulci')])])
flow.connect([(fundus_segmentation, positions, [('segmented_fundi', 'segmented_fundi')])])
flow.connect([(pit_extraction, positions, [('pits', 'pits')])])
flow.connect([(midaxis_segmentation, positions, [('segmented_midaxis', 'segmented_midaxis')])])

flow.connect([(sulcus_segmentation, extents, [('segmented_sulci', 'segmented_sulci')])])
flow.connect([(fundus_segmentation, extents, [('segmented_fundi', 'segmented_fundi')])])
flow.connect([(midaxis_segmentation, extents, [('segmented_midaxis', 'segmented_midaxis')])])

flow.connect([(sulcus_segmentation, curvatures, [('segmented_sulci', 'segmented_sulci')])])
flow.connect([(fundus_segmentation, curvatures, [('segmented_fundi', 'segmented_fundi')])])
flow.connect([(pit_extraction, curvatures, [('pits', 'pits')])])
flow.connect([(midaxis_segmentation, curvatures, [('segmented_midaxis', 'segmented_midaxis')])])

flow.connect([(sulcus_segmentation, depths, [('segmented_sulci', 'segmented_sulci')])])
flow.connect([(fundus_segmentation, depths, [('segmented_fundi', 'segmented_fundi')])])
flow.connect([(pit_extraction, depths, [('pits', 'pits')])])
flow.connect([(midaxis_segmentation, depths, [('segmented_midaxis', 'segmented_midaxis')])])

flow.connect([(sulcus_segmentation, spectra, [('segmented_sulci', 'segmented_sulci')])])
flow.connect([(fundus_segmentation, spectra, [('segmented_fundi', 'segmented_fundi')])])
flow.connect([(midaxis_segmentation, spectra, [('segmented_midaxis', 'segmented_midaxis')])])

##############################################################################
#    Store surface maps, features, and measures in database
##############################################################################

# Database nodes
maps_database = pe.Node(util.Function(input_names = ['depth_map',
                                                     'mean_curv_map',
                                                     'gauss_curv_map'],
                                      output_names=['success'],
                                      function = write_maps_to_database),
                        name='Write_maps_to_database')

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
flow.connect([(surface_maps, maps_database, [('depth_map','depth_map'),
                                ('mean_curv_map','mean_curv_map'),
                                ('gauss_curv_map','gauss_curv_map')])])

# Connect feature to database nodes
flow.connect([(sulcus_segmentation, features_database, [('segmented_sulci', 'segmented_sulci')]),
              (fundus_segmentation, features_database, [('segmented_fundi', 'segmented_fundi')]),
              (pit_extraction, features_database, [('pits', 'pits')]),
              (midaxis_segmentation, features_database, 
                          [('segmented_midaxis', 'segmented_midaxis')])])

# Connect feature measures to database nodes
flow.connect([(positions, measures_database, [('positions_sulci', 'positions_sulci'),
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
flow.connect([(positions, measures_database, [('positions_patches', 'positions_patches'),
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
flow.connect([(measures_database, measures_table, [('measures', 'measures')])])

"""

##############################################################################
#    Run workflow
##############################################################################
if __name__== '__main__':

    flow.write_graph(graph2use='flat')
    flow.write_graph(graph2use='hierarchical')
    flow.run()
