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

def extract_features(surface_file, curvature_file):
    """Extract features: sulcus, medial surface, fundi, and pits

    extract_features('../data/lh.curv','../data/lh.pial')
    """
    from glob import glob
    import subprocess as sp
    cmd = '/projects/mindboggle/mindboggle/extract/fundi/extract.py'
    cmd = ['python', cmd, '%s'%curvature_file, '%s'%surface_file]
    proc = sp.Popen(cmd)
    o, e = proc.communicate()
    if proc.returncode > 0 :
        raise Exception('\n'.join(['extract.py failed', o, e]))
    #output_file = glob('file1.vtk').pop()
    feature_files = glob('*.vtk')
    return feature_files

def measure_position(feature):
    """Measure

    measure_()
    """
    from measure.py import measure_position
    if type(feature) is np.ndarray:
        measurement = measure_position(feature)
    return measurement

def measure_extent(feature):
    """Measure

    measure_()
    """
    from measure.py import measure_extent
    if type(feature) is np.ndarray:
        measurement = measure_extent(feature)
    return measurement

def measure_depth(feature):
    """Measure

    measure_()
    """
    from measure.py import measure_depth
    if type(feature) is np.ndarray:
        measurement = measure_(feature)
    return measurement

def measure_curvature(feature):
    """Measure

    measure_()
    """
    from measure.py import measure_curvature
    if type(feature) is np.ndarray:
        measurement = measure_curvature(feature)
    return measurement

def measure_depth(feature):
    """Measure

    measure_()
    """
    from measure.py import measure_depth
    if type(feature) is np.ndarray:
        measurement = measure_depth(feature)
    return measurement

def measure_LaplaceBeltrami(feature):
    """Measure

    measure_()
    """
    from measure.py import measure_LaplaceBeltrami
    if type(feature) is np.ndarray:
        measurement = measure_LaplaceBeltrami(feature)
    return measurement

def write_features_to_database(feature):
    """Write to database

    write_to_database()
    """
    from write_to_database.py import features_to_database
    if type(feature) is np.ndarray:
        features_to_database(feature)
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


