#!/usr/bin/python

"""
This is Mindboggle's NiPype pipeline!

Example usage:

from pipeline import create_mindboggle_flow
wf = create_mindboggle_flow()
wf.inputs.feature_extractor.curvature_file = '/home/arno/Documents/Projects/mindboggle/data/ManualSurfandVolLabels/subjects/KKI2009-14/surf/lh.curv'
wf.inputs.feature_extractor.surface_file = '/home/arno/Documents/Projects/mindboggle/data/ManualSurfandVolLabels/subjects/KKI2009-14/surf/lh.pial'
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
    cmd = '/home/arno/Documents/Projects/mindboggle/mindboggle/extract/fundi/extract.py'
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

"""
def measure_(object):
    ""Measure

    measure_()
    ""
    from measure.py import measure_
    if type(object) is np.ndarray:
        measurement = measure_(object)
    return measurement
"""

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

def create_mindboggle_flow():
    """Create a Mindboggle workflow"""

    # Create an instance of a workflow and indicate that it should operate in the current directory:
    flow = pe.Workflow(name='pipeline')
    flow.base_dir = '.'

    # Define nodes
    extract_features_node = pe.Node(util.Function(input_names = ['surface_file', 'curvature_file'],
                                              output_names=['feature_files'],
                                              function = extract_features),
                                name='extract_features_node')
    measure_position_node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=['position'],
                                              function = measure_position),
                                       iterfield=['object'],
                                name='measure_position_node')
    measure_extent_node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=['extent'],
                                              function = measure_extent),
                                       iterfield=['object'],
                                name='measure_extent_node')
    measure_curvature_node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=['curvature'],
                                              function = measure_curvature),
                                       iterfield=['object'],
                                name='measure_curvature_node')
    measure_LaplaceBeltrami_node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=['LaplaceBeltrami'],
                                              function = measure_LaplaceBeltrami),
                                       iterfield=['object'],
                                name='measure_LaplaceBeltrami_node')
    """
    measure__node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=[''],
                                              function = measure_),
                                       iterfield=['object'],
                                name='position_')
    """
    write_features_to_database_node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=['success'],
                                              function = write_features_to_database),
                                       iterfield=['object'],
                                name='write_features_to_database_node')
    write_measures_to_database_node = pe.MapNode(util.Function(input_names = ['measurement'],
                                              output_names=['success'],
                                              function = write_measures_to_database),
                                       iterfield=['measurement'],
                                name='write_measures_to_database_node')

# Define connections
    #flow.add_nodes([extract_features_node])

    flow.connect(extract_features_node, 'feature_files', write_features_to_database_node, 'object')

    flow.connect(extract_features_node, 'feature_files', measure_position_node, 'object')
    flow.connect(extract_features_node, 'feature_files', measure_extent_node, 'object')
    flow.connect(extract_features_node, 'feature_files', measure_curvature_node, 'object')
    flow.connect(extract_features_node, 'feature_files', measure_LaplaceBeltrami_node, 'object')
    #flow.connect(extract_features_node, 'feature_files', measure__node, 'object')

    """
    flow.connect(measure_position_node, 'position', write_measures_to_database_node, 'measurement')
    flow.connect(measure_extent_node, 'extent', write_measures_to_database_node, 'measurement')
    flow.connect(measure_depth_node, 'depth', write_measures_to_database_node, 'measurement')
    flow.connect(measure_curvature_node, 'curvature', write_measures_to_database_node, 'measurement')
    flow.connect(measure_LaplaceBeltrami_node, 'LaplaceBeltrami', write_measures_to_database_node, 'measurement')
    #flow.connect(measure__node, '', write_measures_to_database_node, 'measurement')
    """

    return flow

if __name__ == '__main__':
    flow = create_mindboggle_flow()
    flow.inputs.extract_features_node.curvature_file = '/home/arno/Documents/Projects/mindboggle/data/ManualSurfandVolLabels/subjects/KKI2009-14/surf/lh.curv'
    flow.inputs.extract_features_node.surface_file = '/home/arno/Documents/Projects/mindboggle/data/ManualSurfandVolLabels/subjects/KKI2009-14/surf/lh.pial'
    flow.write_graph()
    #flow.run()

