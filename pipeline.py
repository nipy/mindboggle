#!/usr/bin/python

"""
This is Mindboggle's NiPype pipeline!

Example usage:

<<<<<<< HEAD
from pipeline import create_mindboggle_flow
wf = create_mindboggle_flow()
wf.inputs.feature_extractor.curvature_file = '/projects/mindboggle/data/ManualSurfandVolLabels/subjects/KKI2009-14/surf/lh.curv'
wf.inputs.feature_extractor.surface_file = '/projects/mindboggle/data/ManualSurfandVolLabels/subjects/KKI2009-14/surf/lh.pial'
=======
from pipeline import create_mindboggle_workflow
wf = create_mindboggle_workflow()
subject_path = '../data/ManualSurfandVolLabels/subjects/KKI2009-14/'
wf.inputs.feature_extractor.curvature_file = subjects_path + 'surf/lh.curv'
wf.inputs.feature_extractor.surface_file = subjects_path + 'surf/lh.pial'
>>>>>>> 19781b6f52136eced17165422f48fb329e1e19ef
wf.run() # doctest: +SKIP


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.utility as util     # utility
import nipype.interfaces.io as nio
import numpy as np
from pipeline_functions import *

<<<<<<< HEAD
def create_workflow():
    """Create a Mindboggle workflow"""

    # Create an instance of a workflow and 
    # indicate that it should operate in the current directory:
    flow = pe.Workflow(name='pipeline')
    flow.base_dir = '.'
=======
def extract_features(surface_file, curvature_file):
    """Extract features: sulcus, medial surface, fundi, and pits

    extract_features('../data/lh.curv','../data/lh.pial')
    """
    from glob import glob
    import subprocess as sp
    cmd = 'extract/fundi/extract.py'
    cmd = ['python', cmd, '%s'%curvature_file, '%s'%surface_file]
    proc = sp.Popen(cmd)
    o, e = proc.communicate()
    if proc.returncode > 0 :
        raise Exception('\n'.join(['extract.py failed', o, e]))
    #output_file = glob('file1.vtk').pop()
    feature_files = glob('*.vtk')
    return feature_files

def measure_position(object):
    """Measure

    measure_()
    """
    from measure.py import measure_position
    if type(object) is np.ndarray:
        measurement = measure_position(object)
    return measurement

def measure_extent(object):
    """Measure

    measure_()
    """
    from measure.py import measure_extent
    if type(object) is np.ndarray:
        measurement = measure_extent(object)
    return measurement

def measure_depth(object):
    """Measure

    measure_()
    """
    from measure.py import measure_depth
    if type(object) is np.ndarray:
        measurement = measure_(object)
    return measurement

def measure_curvature(object):
    """Measure

    measure_()
    """
    from measure.py import measure_curvature
    if type(object) is np.ndarray:
        measurement = measure_curvature(object)
    return measurement

def measure_LaplaceBeltrami(object):
    """Measure

    measure_()
    """
    from measure.py import measure_LaplaceBeltrami
    if type(object) is np.ndarray:
        measurement = measure_LaplaceBeltrami(object)
    return measurement

def write_features_to_database(object):
    """Write to database

    write_to_database()
    """
    from write_to_database.py import features_to_database
    if type(object) is np.ndarray:
        features_to_database(object)
    success = 'True'
    return success

def write_measures_to_database(measurement):
    """Write to database

    write_to_database()
    """
    from write_to_database.py import measures_to_database
    if type(measurement) is np.ndarray:
        measures_to_database(measurement)
    success = 'True'
    return success

def create_mindboggle_workflow():
    """Create a Mindboggle workflow"""

    # Create an instance of a workflow and indicate that it should operate in the current directory:
    workflow = pe.Workflow(name='pipeline')
    workflow.base_dir = '.'
>>>>>>> 19781b6f52136eced17165422f48fb329e1e19ef

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
<<<<<<< HEAD

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
=======
    extract_features_node = pe.Node(util.Function(input_names = ['surface_file', 'curvature_file'],
                                              output_names=['feature_files'],
                                              function = extract_features),
                                name='ExtractFeatures')
    measure_position_node = pe.Node(util.Function(input_names = ['object'],
                                              output_names=['position'],
                                              function = measure_position),
                                name='MeasurePosition')
    measure_extent_node = pe.Node(util.Function(input_names = ['object'],
                                              output_names=['extent'],
                                              function = measure_extent),
                                name='MeasureExtent')
    measure_curvature_node = pe.Node(util.Function(input_names = ['object'],
                                              output_names=['curvature'],
                                              function = measure_curvature),
                                name='MeasureCurvature')
    measure_LaplaceBeltrami_node = pe.Node(util.Function(input_names = ['object'],
                                              output_names=['LaplaceBeltrami'],
                                              function = measure_LaplaceBeltrami),
                                name='MeasureLaplaceBeltrami')
    write_features_to_database_node = pe.Node(util.Function(input_names = ['object'],
                                              output_names=['success'],
                                              function = write_features_to_database),
                                name='WriteFeaturesToDatabase')
    write_measures_to_database_node = pe.Node(util.Function(input_names = ['measurement'],
                                              output_names=['success'],
                                              function = write_measures_to_database),
                                name='WriteMeasuresToDatabase')

# Define connections
    #workflow.add_nodes([extract_features_node])
    workflow.connect(extract_features_node, 'feature_files', measure_position_node, 'object')
    workflow.connect(extract_features_node, 'feature_files', measure_extent_node, 'object')
    workflow.connect(extract_features_node, 'feature_files', measure_curvature_node, 'object')
    workflow.connect(extract_features_node, 'feature_files', measure_LaplaceBeltrami_node, 'object')

    workflow.connect(extract_features_node, 'feature_files', write_features_to_database_node, 'object')

    workflow.connect(measure_position_node, 'position', write_measures_to_database_node, 'measurement')
    #workflow.connect(measure_extent_node, 'extent', write_measures_to_database_node, 'measurement')
    """
    workflow.connect(measure_depth_node, 'measurement', write_measures_to_database_node, 'measurement')
    workflow.connect(measure_curvature_node, 'measurement', write_measures_to_database_node, 'measurement')
    workflow.connect(measure_LaplaceBeltrami_node, 'measurement', write_measures_to_database_node, 'measurement')
    #workflow.connect(measure__node, '', write_measures_to_database_node, 'measurement')
    """
    return workflow

if __name__ == '__main__':
    subject_path = '../data/ManualSurfandVolLabels/subjects/KKI2009-14/'
    workflow = create_mindboggle_workflow()
    workflow.inputs.ExtractFeatures.curvature_file = subject_path + 'surf/lh.curv'
    workflow.inputs.ExtractFeatures.surface_file = subject_path + 'surf/lh.pial'
    workflow.write_graph()
    #workflow.run()
>>>>>>> 19781b6f52136eced17165422f48fb329e1e19ef


