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

def create_workflow():
    """Create a Mindboggle workflow"""

    # Create an instance of a workflow and 
    # indicate that it should operate in the current directory:
    flow = pe.Workflow(name='pipeline')
    flow.base_dir = '.'

    #def create_extraction_flow():
    #    """Create a feature extraction workflow"""
    #    # Create an instance of a workflow
    #    extraction_flow = pe.Workflow(name='FeatureExtraction')
    #    extraction_flow.base_dir = '.'
    #def create_measurement_workflow():
    #    """Create a shape measurement workflow"""
    #    # Create an instance of a workflow
    #    measurement_flow = pe.Workflow(name='Measurement')
    #    measurement_flow.base_dir = '.'
    #def create_database_workflow():
    #    """Create a write-to-database workflow"""
    #    # Create an instance of a workflow
    #    database_flow = pe.Workflow(name='Database')
    #    database_flow.base_dir = '.'

    # Define nodes

    # DataGrabber node
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['surface_file',
                                                              'curvature_file']),
                         name = 'Data')
    # Feature extraction node
    features = pe.MapNode(util.Function(input_names=['surface_file',
                                                     'curvature_file'],
                                        output_names=['feature_files'],
                                        function = extract_features),
                          iterfield=['surface_file','curvature_file'],
                          name='Extract_features')
    # Shape measurement nodes
    position = pe.MapNode(util.Function(input_names = ['feature'],
                                        output_names=['position'],
                                        function = measure_position),
                          iterfield=['feature'],
                          name='Measure_position')
    extent = pe.MapNode(util.Function(input_names = ['feature'],
                                      output_names=['extent'],
                                      function = measure_extent),
                        iterfield=['feature'],
                        name='Measure_extent')
    curvature = pe.MapNode(util.Function(input_names = ['feature'],
                                         output_names=['curvature'],
                                         function = measure_curvature),
                           iterfield=['feature'],
                           name='Measure_curvature')
    depth = pe.MapNode(util.Function(input_names = ['feature'],
                                     output_names=['depth'],
                                     function = measure_depth),
                           iterfield=['feature'],
                           name='Measure_depth')
    spectra = pe.MapNode(util.Function(input_names = ['feature'],
                                       output_names=['spectra'],
                                       function = measure_LaplaceBeltrami),
                         iterfield=['feature'],
                         name='Measure_spectra')
    # Database nodes
    featuresDB = pe.MapNode(util.Function(input_names = ['feature'],
                                          output_names=['success'],
                                          function = write_features_to_database),
                            iterfield=['feature'],
                            name='Store_features')
    measuresDB = pe.MapNode(util.Function(input_names = ['position',
                                                         'extent',
                                                         'curvature',
                                                         'depth',
                                                         'spectra'],
                                          output_names=['success'],
                                          function = write_measures_to_database),
                            iterfield=['measurement'],
                            name='Store_measures')

    # Define connections
    flow.connect([(datasource,features,[('surface_file','surface_file'),
                                            ('curvature_file','curvature_file')]
                    )])
    flow.connect(features, 'feature_files', position, 'feature')
    flow.connect(features, 'feature_files', extent, 'feature')
    flow.connect(features, 'feature_files', curvature, 'feature')
    flow.connect(features, 'feature_files', depth, 'feature')
    flow.connect(features, 'feature_files', spectra, 'feature')

    flow.connect(features, 'feature_files', featuresDB, 'feature')

    flow.connect(position, 'position', measuresDB, 'position')
    flow.connect(extent, 'extent', measuresDB, 'extent')
    flow.connect(curvature, 'curvature', measuresDB, 'curvature')
    flow.connect(depth, 'depth', measuresDB, 'depth')
    flow.connect(spectra, 'spectra', measuresDB, 'spectra')
    return flow

if __name__ == '__main__':
    flow = create_workflow()

    # Initiate the DataGrabber node with the infield: 'subject_id'
    # and the outfield: 'curv' and 'pial'
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['surface_file',
                                                              'curvature_file']),
                              name = 'DataSource')
    # Specify the location of the data folder
    datasource.inputs.base_directory = '/projects/mindboggle/data/ManualSurfandVolLabels/subjects'

    # Define the structure of the data folders and files.
    # Each '%s' will later be filled by a template argument.
    datasource.inputs.template = '%s/surf/lh.%s'

    # Define the arguments for the template
    datasource.inputs.template_args['surface_file'] = [['subject_id', 'pial']]
    datasource.inputs.template_args['curvature_file'] = [['subject_id', 'curv']]

    # Define the list of subjects
    datasource.iterables = ('subject_id', ['KKI2009-14'])

    # Write graphs
    flow.write_graph(graph2use='flat')
    flow.write_graph(graph2use='hierarchical')
    #flow.run()


