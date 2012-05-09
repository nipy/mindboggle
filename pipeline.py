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

"""
def create_extraction_workflow():
    Create a feature extraction workflow
    # Create an instance of a workflow
    flow = pe.Workflow(name='Feature_extraction_pipeline')
    flow.base_dir = '.'

    # Define nodes

    # Feature extraction node
    features = pe.MapNode(util.Function(input_names=['surface_file',
                                                     'curvature_file'],
                                        output_names=['feature_files'],
                                        function = extract_features),
                          iterfield=['surface_file','curvature_file'],
                          name='Extract_features')

    #def create_extraction_flow():
    #    Create a feature extraction workflow
    #    # Create an instance of a workflow
    #    extraction_flow = pe.Workflow(name='FeatureExtraction')
    #    extraction_flow.base_dir = '.'
    #def create_measurement_workflow():
    #    Create a shape measurement workflow
    #    # Create an instance of a workflow
    #    measurement_flow = pe.Workflow(name='Measurement')
    #    measurement_flow.base_dir = '.'
    #def create_database_workflow():
    #    Create a write-to-database workflow
    #    # Create an instance of a workflow
    #    database_flow = pe.Workflow(name='Database')
    #    database_flow.base_dir = '.'
"""

def create_workflow():
    """Create a Mindboggle workflow"""

    # Create an instance of a workflow and 
    # indicate that it should operate in the current directory:
    flow = pe.Workflow(name='pipeline')
    flow.base_dir = '.'

    # Define nodes

    # DataGrabber node
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['surface_file']),
                         name = 'Data')

    # Measure surface maps node
    maps = pe.MapNode(util.Function(input_names = ['surface_file'],
                                    output_names=['depth_map',
                                                  'mean_curvature_map',
                                                  'gauss_curvature_map'],
                                    function = measure_depth),
                      iterfield=['surface_file'],
                      name='Measure_surface_maps')

    # Feature extraction nodes
    sulci = pe.MapNode(util.Function(input_names=['surface_file',
                                                  'depth_map',
                                                  'mean_curvature_map',
                                                  'gauss_curvature_map'],
                                     output_names=['sulci'],
                                     function = extract_sulci),
                       iterfield=['surface_file'],
                       name='Extract_sulci')
    fundi = pe.MapNode(util.Function(input_names=['surface_file',
                                                  'depth_map',
                                                  'mean_curvature_map',
                                                  'gauss_curvature_map'],
                                     output_names=['fundi'],
                                     function = extract_fundi),
                       iterfield=['surface_file'],
                       name='Extract_fundi')
    pits = pe.MapNode(util.Function(input_names=['surface_file',
                                                 'depth_map',
                                                 'mean_curvature_map',
                                                 'gauss_curvature_map'],
                                    output_names=['pits'],
                                    function = extract_pits),
                      iterfield=['surface_file'],
                      name='Extract_pits')
    medial = pe.MapNode(util.Function(input_names=['surface_file',
                                                   'depth_map',
                                                   'mean_curvature_map',
                                                   'gauss_curvature_map'],
                                      output_names=['medial'],
                                      function = extract_medial),
                        iterfield=['surface_file'],
                        name='Extract_medial')

    # Multi-atlas registration
    """
    registration = pe.MapNode(util.Function(input_names=['surface_file',
                                                         'curvature_file'],
                                            output_names=['feature_files'],
                                            function = register),
                              iterfield=['surface_file','curvature_file'],
                              name='Register_atlases')
    """
    # Shape measurement nodes
    position = pe.MapNode(util.Function(input_names = ['sulci','fundi',
                                                       'pits','medial'],
                                        output_names=['position'],
                                        function = measure_position),
                          iterfield=['sulci'],
                          name='Measure_position')
    extent = pe.MapNode(util.Function(input_names = ['sulci','fundi',
                                                     'pits','medial'],
                                      output_names=['extent'],
                                      function = measure_extent),
                        iterfield=['feature'],
                        name='Measure_extent')
    curvature = pe.MapNode(util.Function(input_names = ['sulci','fundi',
                                                        'pits','medial'],
                                         output_names=['curvature'],
                                         function = measure_curvature),
                           iterfield=['feature'],
                           name='Measure_curvature')
    depth = pe.MapNode(util.Function(input_names = ['sulci','fundi',
                                                    'pits','medial'],
                                     output_names=['depth'],
                                     function = measure_depth),
                           iterfield=['feature'],
                           name='Measure_depth')
    spectra = pe.MapNode(util.Function(input_names = ['sulci','fundi',
                                                      'pits','medial'],
                                       output_names=['spectra'],
                                       function = measure_LaplaceBeltrami),
                         iterfield=['feature'],
                         name='Measure_spectra')

    # Database nodes
    featureDB = pe.MapNode(util.Function(input_names = ['sulci','fundi',
                                                        'pits','medial'],
                                         output_names=['success'],
                                         function = write_features_to_database),
                           iterfield=['feature'],
                           name='Store_features')
    measureDB = pe.MapNode(util.Function(input_names = ['position',
                                                        'extent',
                                                        'curvature',
                                                        'depth',
                                                        'spectra'],
                                         output_names=['success'],
                                         function = write_measures_to_database),
                           iterfield=['measurement'],
                           name='Store_measures')

    # Define connections
    flow.connect(datasource, 'surface_file', maps, 'surface_file')
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
    flow.connect([(sulci,  position,  [('sulci', 'sulci')])])
    flow.connect([(fundi,  position,  [('fundi', 'fundi')])])
    flow.connect([(pits,   position,  [('pits',  'pits')])])
    flow.connect([(medial, position,  [('medial','medial')])])
    flow.connect([(sulci,  extent,    [('sulci', 'sulci')])])
    flow.connect([(fundi,  extent,    [('fundi', 'fundi')])])
    flow.connect([(pits,   extent,    [('pits',  'pits')])])
    flow.connect([(medial, extent,    [('medial','medial')])])
    flow.connect([(sulci,  curvature, [('sulci', 'sulci')])])
    flow.connect([(fundi,  curvature, [('fundi', 'fundi')])])
    flow.connect([(pits,   curvature, [('pits',  'pits')])])
    flow.connect([(medial, curvature, [('medial','medial')])])
    flow.connect([(sulci,  depth,     [('sulci', 'sulci')])])
    flow.connect([(fundi,  depth,     [('fundi', 'fundi')])])
    flow.connect([(pits,   depth,     [('pits',  'pits')])])
    flow.connect([(medial, depth,     [('medial','medial')])])
    flow.connect([(sulci,  spectra,   [('sulci', 'sulci')])])
    flow.connect([(fundi,  spectra,   [('fundi', 'fundi')])])
    flow.connect([(pits,   spectra,   [('pits',  'pits')])])
    flow.connect([(medial, spectra,   [('medial','medial')])])

    flow.connect(sulci,  'sulci',  featureDB, 'sulci')
    flow.connect(fundi,  'fundi',  featureDB, 'fundi')
    flow.connect(pits,   'pits',   featureDB, 'pits')
    flow.connect(medial, 'medial', featureDB, 'medial')

    flow.connect(position,  'position',  measureDB, 'position')
    flow.connect(extent,    'extent',    measureDB, 'extent')
    flow.connect(curvature, 'curvature', measureDB, 'curvature')
    flow.connect(depth,     'depth',     measureDB, 'depth')
    flow.connect(spectra,   'spectra',   measureDB, 'spectra')
    return flow


if __name__ == '__main__':
    flow = create_workflow()

    # Initiate the DataGrabber node with the infield: 'subject_id'
    # and the outfield: 'curv' and 'pial'
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['surface_file']),
                              name = 'DataSource')
    # Specify the location of the data folder
    datasource.inputs.base_directory = '/projects/mindboggle/data/ManualSurfandVolLabels/subjects'

    # Define the structure of the data folders and files.
    # Each '%s' will later be filled by a template argument.
    datasource.inputs.template = '%s/surf/lh.%s'

    # Define the arguments for the template
    datasource.inputs.template_args['surface_file'] = [['subject_id', 'pial']]

    # Define the list of subjects
    datasource.iterables = ('subject_id', ['KKI2009-14'])

    # Write graphs
    flow.write_graph(graph2use='flat')
    flow.write_graph(graph2use='hierarchical')
    #flow.run()


