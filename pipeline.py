#!/usr/bin/python

"""
This is Mindboggle's NiPype pipeline!

Example usage:

from pipeline import create_mindboggle_workflow
wf = create_mindboggle_workflow()
subject_path = '../data/ManualSurfandVolLabels/subjects/KKI2009-14/'
wf.inputs.feature_extractor.curvature_file = subjects_path + 'surf/lh.curv'
wf.inputs.feature_extractor.surface_file = subjects_path + 'surf/lh.pial'
wf.run() # doctest: +SKIP


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
import numpy as np

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

    # Define nodes
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

