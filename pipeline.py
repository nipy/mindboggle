#!/usr/bin/python

"""
This is Mindboggle's NiPype pipeline!

Example usage:

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

##############################################################################
#   Feature extraction
##############################################################################

# Measure surface maps node
maps = pe.Node(util.Function(input_names = ['surface_file'],
                             output_names = ['depth_map',
                                             'mean_curvature_map',
                                             'gauss_curvature_map'],
                             function = measure_depth),
               name='Measure_surface_maps')

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
medial = pe.Node(util.Function(input_names = ['depth_map',
                                              'mean_curvature_map',
                                              'gauss_curvature_map'],
                               output_names = ['medial'],
                               function = extract_medial),
                 name='Extract_medial')

# Connect input to surface maps node
flow.connect([(infosource, datasource, [('subject_id','subject_id')])])
flow.connect([(datasource, maps, [('surface_file','surface_file')])])

# Connect surface maps node to feature extraction nodes
flow.connect([(maps, sulci, [('depth_map','depth_map'),
                             ('mean_curvature_map','mean_curvature_map'),
                             ('gauss_curvature_map','gauss_curvature_map')])])
flow.connect([(maps, fundi, [('depth_map','depth_map'),
                             ('mean_curvature_map','mean_curvature_map'),
                             ('gauss_curvature_map','gauss_curvature_map')])])
flow.connect([(maps, pits, [('depth_map','depth_map'),
                            ('mean_curvature_map','mean_curvature_map'),
                            ('gauss_curvature_map','gauss_curvature_map')])])
flow.connect([(maps, medial, [('depth_map','depth_map'),
                              ('mean_curvature_map','mean_curvature_map'),
                              ('gauss_curvature_map','gauss_curvature_map')])])

# Save output
flow.connect([(maps, datasink, [('depth_map','depth_map'),
                                ('mean_curvature_map','mean_curvature_map'),
                                ('gauss_curvature_map','gauss_curvature_map')])])
flow.connect([(sulci,  datasink, [('sulci', 'sulci')]),
              (fundi,  datasink, [('fundi', 'fundi')]),
              (pits,   datasink, [('pits',  'pits')]),
              (medial, datasink, [('medial','medial')])])

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
                                                  'subject_surf_path', 
                                                  'template_path', 
                                                  'template_name', 
                                                  'registration_name'],
                                     output_names=['registration_name'],
                                     function = register_template),
                       name='Register_template')

atlas_registration = pe.Node(util.Function(input_names=['subject_id',
                                                         'atlas_list_file',
                                                         'registration_name',
                                                         'output_path'],
                                            output_names=['atlas_list'],
                                            function = register_atlases),
                             name='Register_atlases')

# Connect input to registration nodes
flow.connect([(datasource, template_registration, [('subject_id','subject_id')])])

# Connect registration and multiatlas labeling nodes
flow.connect([(template_registration, atlas_registration, 
               [('registration_name','registration_name')])])

##############################################################################
#   Feature identification
##############################################################################

# Feature identification node
"""
"""
fundus_identification = pe.Node(util.Function(input_names=['fundi','atlas_list'],
                                     output_names=['labeled_fundi'],
                                     function = identify_fundi),
                       name='Identify_fundi')

# Connect multi-atlas registration-based labeling and feature identification nodes
flow.connect([(atlas_registration, fundus_identification,
               [('atlas_list','atlas_list')]),
              (fundi, fundus_identification,
               [('fundi','fundi')])])

#flow.connect([(labeled_fundi, fundus_identification, [('fundi', 'labeled_fundi')])])

##############################################################################
#   Shape measurement
##############################################################################

# Shape measurement nodes
position = pe.MapNode(util.Function(input_names = ['labeled_sulci','labeled_fundi',
                                                   'labeled_pits','labeled_medial'],
                                    output_names=['position'],
                                    function = measure_position),
                      iterfield=['sulci'],
                      name='Measure_position')
extent = pe.MapNode(util.Function(input_names = ['labeled_sulci','labeled_fundi',
                                                 'labeled_pits','labeled_medial'],
                                  output_names=['extent'],
                                  function = measure_extent),
                    iterfield=['feature'],
                    name='Measure_extent')
curvature = pe.MapNode(util.Function(input_names = ['labeled_sulci','labeled_fundi',
                                                    'labeled_pits','labeled_medial'],
                                     output_names=['curvature'],
                                     function = measure_curvature),
                       iterfield=['feature'],
                       name='Measure_curvature')
depth = pe.MapNode(util.Function(input_names = ['labeled_sulci','labeled_fundi',
                                                'labeled_pits','labeled_medial'],
                                 output_names=['depth'],
                                 function = measure_depth),
                       iterfield=['feature'],
                       name='Measure_depth')
spectra = pe.MapNode(util.Function(input_names = ['labeled_sulci','labeled_fundi',
                                                  'labeled_pits','labeled_medial'],
                                   output_names=['spectra'],
                                   function = measure_LaplaceBeltrami),
                     iterfield=['feature'],
                     name='Measure_spectra')

# Connect feature to shape measurement nodes
#flow.connect([(sulcus_identification,  position,  [('labeled_sulci', 'labeled_sulci')])])
flow.connect([(fundus_identification,  position,  [('labeled_fundi', 'labeled_fundi')])])
#flow.connect([(pit_identification,     position,  [('labeled_pits',  'labeled_pits')])])
#flow.connect([(medial_identification,  position,  [('labeled_medial','labeled_medial')])])
#flow.connect([(sulcus_identification,  extent,    [('labeled_sulci', 'labeled_sulci')])])
flow.connect([(fundus_identification,  extent,    [('labeled_fundi', 'labeled_fundi')])])
#flow.connect([(pit_identification,     extent,    [('labeled_pits',  'labeled_pits')])])
#flow.connect([(medial_identification,  extent,    [('labeled_medial','labeled_medial')])])
#flow.connect([(sulcus_identification,  curvature, [('labeled_sulci', 'labeled_sulci')])])
flow.connect([(fundus_identification,  curvature, [('labeled_fundi', 'labeled_fundi')])])
#flow.connect([(pit_identification,     curvature, [('labeled_pits',  'labeled_pits')])])
#flow.connect([(medial_identification,  curvature, [('labeled_medial','labeled_medial')])])
#flow.connect([(sulcus_identification,  depth,     [('labeled_sulci', 'labeled_sulci')])])
flow.connect([(fundus_identification,  depth,     [('labeled_fundi', 'labeled_fundi')])])
#flow.connect([(pit_identification,     depth,     [('labeled_pits',  'labeled_pits')])])
#flow.connect([(medial_identification,  depth,     [('labeled_medial','labeled_medial')])])
#flow.connect([(sulcus_identification,  spectra,   [('labeled_sulci', 'labeled_sulci')])])
flow.connect([(fundus_identification,  spectra,   [('labeled_fundi', 'labeled_fundi')])])
#flow.connect([(pit_identification,     spectra,   [('labeled_pits',  'labeled_pits')])])
#flow.connect([(medial_identification,  spectra,   [('labeled_medial','labeled_medial')])])

##############################################################################
#    Store features in database
##############################################################################

"""
# Database nodes
database = pe.MapNode(util.Function(input_names = ['depth_map',
                                                   'mean_curvature_map',
                                                   'gauss_curvature_map',
                                                   'labeled_sulci',
                                                   'labeled_fundi',
                                                   'labeled_pits',
                                                   'labeled_medial',
                                                   'position',
                                                   'extent',
                                                   'curvature',
                                                   'depth',
                                                   'spectra'],
                                    output_names=['success'],
                                    function = write_to_database),
                      iterfield=['storees'],
                      name='Write_to_database')

# Connect maps to database nodes
flow.connect([(maps, database, [('depth_map','depth_map'),
                                ('mean_curvature_map','mean_curvature_map'),
                                ('gauss_curvature_map','gauss_curvature_map'),
                                ('labeled_sulci','labeled_sulci'),
                                ('labeled_fundi','labeled_fundi'),
                                ('labeled_pits','labeled_pits'),
                                ('labeled_medial','labeled_medial'),
                                ('position','position'),
                                ('extent','extent'),
                                ('curvature','curvature'),
                                ('depth','depth'),
                                ('spectra','spectra')])])

# Connect feature to database nodes
flow.connect(sulci,  'sulci',  featureDB, 'sulci')
flow.connect(fundi,  'fundi',  featureDB, 'fundi')
flow.connect(pits,   'pits',   featureDB, 'pits')
flow.connect(medial, 'medial', featureDB, 'medial')

# Connect measure to database nodes
flow.connect(position,  'position',  measureDB, 'position')
flow.connect(extent,    'extent',    measureDB, 'extent')
flow.connect(curvature, 'curvature', measureDB, 'curvature')
flow.connect(depth,     'depth',     measureDB, 'depth')
flow.connect(spectra,   'spectra',   measureDB, 'spectra')

"""

##############################################################################
#    Run workflow
##############################################################################
if __name__== '__main__':

    flow.write_graph(graph2use='flat')
    flow.write_graph(graph2use='hierarchical')
    #flow.run()

