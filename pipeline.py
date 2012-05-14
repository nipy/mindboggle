#!/usr/bin/python

"""
This is Mindboggle's NiPype pipeline!

Example usage:

FIX:
from pipeline import create_mindboggle_flow
wf = create_mindboggle_flow()
wf.inputs.feature_extractor.curvature_file = '/projects/mindboggle/data/ManualSurfandVolLabels/subjects/KKI2009-14/surf/lh.curv'
wf.inputs.feature_extractor.surface_file = '/projects/mindboggle/data/ManualSurfandVolLabels/subjects/KKI2009-14/surf/lh.pial'
wf.run() # doctest: +SKIP


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.utility as util     # utility
import nipype.interfaces.io as nio
import numpy as np
from pipeline_functions import *

##############################################################################
#   Workflow setup
##############################################################################
flow = pe.Workflow(name='pipeline')
flow.base_dir = '.'

# Input and output nodes
infosource = pe.Node(interface=util.IdentityInterface(fields=['subject_id']),
                     name = 'Input')
datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                               outfields=['surface_file']),
                     name = 'DataSource')
template_source = pe.Node(interface=nio.DataGrabber(infields=['template_id'],
                                                    outfields=['template_surface_file']),
                          name = 'Template')
atlas_source = pe.Node(interface=nio.DataGrabber(infields=['atlas_list_file'],
                                                 outfields=['atlas_list_file']),
                       name = 'Atlases')
datasink = pe.Node(interface=nio.DataSink(),
                   name = 'DataSink')

# Iterate over subjects
infosource.iterables = ('subject_id', ['KKI2009-14'])

# Specify the location and structure of the inputs and outputs
datasource.inputs.base_directory = '/projects/mindboggle/data/ManualSurfandVolLabels/subjects'
datasource.inputs.template = '%s/surf/lh.%s'
datasource.inputs.template_args['surface_file'] = [['subject_id', 'pial']]

datasink.inputs.base_directory = '/projects/mindboggle'
datasink.inputs.container = 'output'

subjects_directory = '/Applications/freesurfer/subjects/'
subject_id = 'bert'
multilabel_directory = 'multilabels'
atlas_list = '../../data/KKI_OASIS.txt' 

# Specify the location and structure of the template
template_source.inputs.template = '%s/surf/lh.%s'
template_source.inputs.template_args['template_surface_file'] = [['template_id', 'pial']]

# Specify the location and structure of the atlases
atlas_source.inputs.template = '%s/%s'
atlas_source.inputs.template_args['atlas_list'] = [['atlases', 'atlas_list.txt']]

##############################################################################
#   Surface map calculation
##############################################################################

# Measure surface surface_maps node
surface_maps = pe.Node(util.Function(input_names = ['surface_file'],
                                     output_names = ['depth_map',
                                                     'mean_curvature_map',
                                                     'gauss_curvature_map'],
                                     function = measure_surface_maps),
                       name='Measure_surface_maps')

##############################################################################
#   Feature extraction
##############################################################################

# Feature extraction nodes
sulci = pe.Node(util.Function(input_names = ['depth_map',
                                             'mean_curvature_map',
                                             'gauss_curvature_map'],
                              output_names = ['sulci'],
                              function = extract_sulci),
                name='Extract_sulci')
fundi = pe.Node(util.Function(input_names = ['depth_map',
                                             'mean_curvature_map',
                                             'gauss_curvature_map'],
                              output_names = ['fundi'],
                              function = extract_fundi),
                name='Extract_fundi')
pits = pe.Node(util.Function(input_names = ['depth_map',
                                            'mean_curvature_map',
                                            'gauss_curvature_map'],
                             output_names = ['pits'],
                             function = extract_pits),
               name='Extract_pits')
medialaxis = pe.Node(util.Function(input_names = ['depth_map',
                                                  'mean_curvature_map',
                                                  'gauss_curvature_map'],
                                   output_names = ['medialaxis'],
                                   function = extract_medialaxis),
                     name='Extract_medialaxis')

# Connect input to surface surface_maps node
flow.connect([(infosource, datasource, [('subject_id','subject_id')])])
flow.connect([(datasource, surface_maps, [('surface_file','surface_file')])])

# Connect surface surface_maps node to feature extraction nodes
flow.connect([(surface_maps, sulci, [('depth_map','depth_map'),
                                     ('mean_curvature_map','mean_curvature_map'),
                                     ('gauss_curvature_map','gauss_curvature_map')])])
flow.connect([(surface_maps, fundi, [('depth_map','depth_map'),
                                     ('mean_curvature_map','mean_curvature_map'),
                                     ('gauss_curvature_map','gauss_curvature_map')])])
flow.connect([(surface_maps, pits, [('depth_map','depth_map'),
                                    ('mean_curvature_map','mean_curvature_map'),
                                    ('gauss_curvature_map','gauss_curvature_map')])])
flow.connect([(surface_maps, medialaxis, [('depth_map','depth_map'),
                                          ('mean_curvature_map','mean_curvature_map'),
                                          ('gauss_curvature_map','gauss_curvature_map')])])

"""
# Save output
flow.connect([(surface_maps, datasink, [('depth_map','depth_map'),
                                ('mean_curvature_map','mean_curvature_map'),
                                ('gauss_curvature_map','gauss_curvature_map')])])
"""

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
flow.connect([(template_source, template_registration, [('template_id','template_id')])])
flow.connect([(atlas_source, atlas_registration, [('atlas_list_file','atlas_list_file')])])

# Connect input to registration nodes
flow.connect([(datasource, template_registration, [('subject_id','subject_id')])])

# Connect template registration and multiatlas registration(-based labeling) nodes
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

# Connect multiatlas registration(-based labeling) and label propagation nodes
flow.connect([(atlas_registration, label_propagation,
               [('labels','labels')]),
              (fundi, label_propagation,
               [('fundi','fundi')])])

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

medialaxis_segmentation = pe.Node(util.Function(input_names=['medialaxis','labels'],
                                                output_names=['segmented_medialaxis'],
                                                function = segment_medialaxis),
                                  name='Segment_medialaxis')

# Connect feature and feature segmentation nodes
flow.connect([(sulci, sulcus_segmentation, [('sulci','sulci')]),
              (fundi, fundus_segmentation, [('fundi','fundi')]),
              (medialaxis, medialaxis_segmentation, [('medialaxis','medialaxis')])])

# Connect multiatlas registration(-based labeling) and feature segmentation nodes
flow.connect([(label_propagation, sulcus_segmentation, [('labels','labels')]),
              (label_propagation, fundus_segmentation, [('labels','labels')]),
              (sulcus_segmentation, medialaxis_segmentation, [('segmented_sulci','labels')])])
              
##############################################################################
#   Shape measurement
##############################################################################

# Shape measurement nodes
positions = pe.Node(util.Function(input_names = ['segmented_sulci', 'segmented_fundi',
                                                 'segmented_medialaxis', 'pits'],
                                  output_names=['positions_sulci', 'positions_fundi',
                                                'positions_medialaxis', 'positions_pits'],
                                  function = measure_positions),
                    name='Measure_positions')

extents = pe.Node(util.Function(input_names = ['segmented_sulci', 'segmented_fundi',
                                               'segmented_medialaxis', 'pits'],
                                output_names=['extents_sulci', 'extents_fundi',
                                              'extents_medialaxis', 'extents_pits'],
                                function = measure_extents),
                    name='Measure_extents')

curvatures = pe.Node(util.Function(input_names = ['segmented_sulci', 'segmented_fundi',
                                                  'segmented_medialaxis', 'pits'],
                                   output_names=['curvatures_sulci', 'curvatures_fundi',
                                                 'curvatures_medialaxis', 'curvatures_pits'],
                                   function = measure_curvatures),
                    name='Measure_curvatures')

depths = pe.Node(util.Function(input_names = ['segmented_sulci', 'segmented_fundi',
                                              'segmented_medialaxis', 'pits'],
                               output_names=['depths_sulci', 'depths_fundi',
                                             'depths_medialaxis', 'depths_pits'],
                               function = measure_depths),
                 name='Measure_depths')

spectra = pe.Node(util.Function(input_names = ['segmented_sulci', 'segmented_fundi',
                                               'segmented_medialaxis', 'pits'],
                                output_names=['spectra_sulci', 'spectra_fundi',
                                              'spectra_medialaxis', 'spectra_pits'],
                                function = measure_spectra),
                  name='Measure_spectra')

# Connect feature to shape measurement nodes
flow.connect([(sulcus_segmentation, positions, [('segmented_sulci', 'segmented_sulci')])])
flow.connect([(fundus_segmentation, positions, [('segmented_fundi', 'segmented_fundi')])])
flow.connect([(pits, positions, [('pits', 'pits')])])
flow.connect([(medialaxis_segmentation, positions, [('segmented_medialaxis', 'segmented_medialaxis')])])

flow.connect([(sulcus_segmentation, extents, [('segmented_sulci', 'segmented_sulci')])])
flow.connect([(fundus_segmentation, extents, [('segmented_fundi', 'segmented_fundi')])])
flow.connect([(medialaxis_segmentation, extents, [('segmented_medialaxis', 'segmented_medialaxis')])])

flow.connect([(sulcus_segmentation, curvatures, [('segmented_sulci', 'segmented_sulci')])])
flow.connect([(fundus_segmentation, curvatures, [('segmented_fundi', 'segmented_fundi')])])
flow.connect([(pits, curvatures, [('pits', 'pits')])])
flow.connect([(medialaxis_segmentation, curvatures, [('segmented_medialaxis', 'segmented_medialaxis')])])

flow.connect([(sulcus_segmentation, depths, [('segmented_sulci', 'segmented_sulci')])])
flow.connect([(fundus_segmentation, depths, [('segmented_fundi', 'segmented_fundi')])])
flow.connect([(pits, depths, [('pits', 'pits')])])
flow.connect([(medialaxis_segmentation, depths, [('segmented_medialaxis', 'segmented_medialaxis')])])

flow.connect([(sulcus_segmentation, spectra, [('segmented_sulci', 'segmented_sulci')])])
flow.connect([(fundus_segmentation, spectra, [('segmented_fundi', 'segmented_fundi')])])
flow.connect([(medialaxis_segmentation, spectra, [('segmented_medialaxis', 'segmented_medialaxis')])])

##############################################################################
#    Store surface maps, features, and measures in database
##############################################################################

# Database nodes
maps_database = pe.Node(util.Function(input_names = ['depth_map',
                                                     'mean_curvature_map',
                                                     'gauss_curvature_map'],
                                      output_names=['success'],
                                      function = write_maps_to_database),
                        name='Write_maps_to_database')

features_database = pe.Node(util.Function(input_names = ['segmented_sulci',
                                                         'segmented_fundi',
                                                         'pits',
                                                         'segmented_medialaxis'],
                                          output_names=['success'],
                                          function = write_features_to_database),
                            name='Write_features_to_database')

measures_database = pe.Node(util.Function(input_names = ['positions_sulci',
                                                         'positions_fundi',
                                                         'positions_pits',
                                                         'positions_medialaxis',
                                                         'extents_sulci',
                                                         'extents_fundi',
                                                         'extents_medialaxis',
                                                         'curvatures_sulci',
                                                         'curvatures_fundi',
                                                         'curvatures_pits',
                                                         'curvatures_medialaxis',
                                                         'depths_sulci',
                                                         'depths_fundi',
                                                         'depths_pits',
                                                         'depths_medialaxis',
                                                         'spectra_sulci',
                                                         'spectra_fundi',
                                                         'spectra_medialaxis'],
                                          output_names=['measures'],
                                          function = write_measures_to_database),
                            name='Write_measures_to_database')

measures_table = pe.Node(util.Function(input_names = ['measures'],
                                       output_names=['success'],
                                       function = write_measures_to_table),
                         name='Write_measures_to_table')

# Connect surface maps to database nodes
flow.connect([(surface_maps, maps_database, [('depth_map','depth_map'),
                                ('mean_curvature_map','mean_curvature_map'),
                                ('gauss_curvature_map','gauss_curvature_map')])])

# Connect feature to database nodes
flow.connect([(sulcus_segmentation, features_database, [('segmented_sulci', 'segmented_sulci')]),
              (fundus_segmentation, features_database, [('segmented_fundi', 'segmented_fundi')]),
              (pits, features_database, [('pits', 'pits')]),
              (medialaxis_segmentation, features_database, 
                          [('segmented_medialaxis', 'segmented_medialaxis')])])

# Connect measure to database nodes
flow.connect([(positions, measures_database, [('positions_sulci', 'positions_sulci'),
                                              ('positions_fundi', 'positions_fundi'),
                                              ('positions_pits', 'positions_pits'),
                                              ('positions_medialaxis', 'positions_medialaxis')]),
              (extents, measures_database, [('extents_sulci', 'extents_sulci'),
                                            ('extents_fundi', 'extents_fundi'),
                                            ('extents_medialaxis', 'extents_medialaxis')]),
              (curvatures, measures_database, [('curvatures_sulci', 'curvatures_sulci'),
                                               ('curvatures_fundi', 'curvatures_fundi'),
                                               ('curvatures_pits', 'curvatures_pits'),
                                               ('curvatures_medialaxis', 'curvatures_medialaxis')]),
              (depths, measures_database, [('depths_sulci', 'depths_sulci'),
                                           ('depths_fundi', 'depths_fundi'),
                                           ('depths_pits', 'depths_pits'),
                                           ('depths_medialaxis', 'depths_medialaxis')]),
              (spectra, measures_database, [('spectra_sulci', 'spectra_sulci'),
                                            ('spectra_fundi', 'spectra_fundi'),
                                            ('spectra_medialaxis', 'spectra_medialaxis')])])

# Connect measure to table nodes
flow.connect([(measures_database, measures_table, [('measures', 'measures')])])

##############################################################################
#    Run workflow
##############################################################################
if __name__== '__main__':

    flow.write_graph(graph2use='flat')
    flow.write_graph(graph2use='hierarchical')
    #flow.run()

