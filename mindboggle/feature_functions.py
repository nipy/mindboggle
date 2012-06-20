#!/usr/bin/python

"""
Feature-based functions for morphometry and labeling:

2. Feature extraction

3. Label propagation

4. Feature segmentation / identification

5. Shape measurement


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

##############################################################################
#   Feature extraction
##############################################################################

def extract_fundi(command, depth_file):
    """
    Extract fundi

    extract_fundi
    """
    from nipype.interfaces.base import CommandLine

    args = [command, depth_file]
    cli = CommandLine(command = 'python')
    cli.inputs.args = ' '.join(args)
    cli.cmdline
    cli.run()

    return fundi

def extract_sulci(surface_file, depth_file, mean_curvature_file, gauss_curvature_file):
    """
    Extract sulci

    extract_sulci    
    """
    from nipype.interfaces.base import CommandLine

    return sulci

def extract_midaxis(surface_file, depth_file, mean_curvature_file, gauss_curvature_file):
    """
    Extract midaxis

    extract_midaxis    
    """
    from nipype.interfaces.base import CommandLine

    return midaxis

# Labeled surface patch and volume extraction nodes
def extract_patches(labels):
    """
    Extract labeled surface patches
    
    extract_patches    
    """
    from nipype.interfaces.base import CommandLine

    return patches

def extract_regions(labels):
    """
    Extract labeled region volumes
    
    extract_regions
    """
    from nipype.interfaces.base import CommandLine

    return regions

##############################################################################
#   Label propagation
##############################################################################

def propagate_labels(labels, fundi):
    """
    Propagate labels
    """

    return labels

# Volume label propagation node
def propagate_volume_labels(labels):
    """
    Propagate labels through volume
    """

    return labels


##############################################################################
#   Feature segmentation / identification
##############################################################################

def segment_sulci(sulci):
    """
    Segment and identify sulci
    """
    from nipype.interfaces.base import CommandLine

    return feature_type, segmented_sulci

def segment_fundi(fundi):
    """
    Segment and identify fundi
    """
    from nipype.interfaces.base import CommandLine

    return feature_type, segmented_fundi

def segment_midaxis(midaxis):
    """
    Segment and identify medial axis surfaces
    """
    from nipype.interfaces.base import CommandLine

    return feature_type, segmented_midaxis

##############################################################################
#   Shape measurement functions
##############################################################################

def measure_position(feature):
    """
    Measure

    measure_()
    """
    from measure.py import measure_position

    measurement = measure_position(feature)

    return measurement

def measure_extent(feature):
    """
    Measure

    measure_()
    """
    from measure.py import measure_extent

    measurement = measure_position(feature)

    return measurement

def measure_curvature(feature):
    """
    Measure

    measure_()
    """
    from measure.py import measure_curvature

    measurement = measure_position(feature)

    return measurement

def measure_depth(feature):
    """
    Measure

    measure_()
    """
    from measure.py import measure_depth

    measurement = measure_position(feature)

    return measurement

def measure_spectra(feature):
    """
    Measure

    measure_()
    """
    from measure.py import measure_spectra

    measurement = measure_position(feature)

    return measurement

##############################################################################
#    Store features in database
##############################################################################

def write_surfaces_to_database(surfaces):
    """
    Write to database

    write_to_database()
    """
    from write_to_database.py import surfaces_to_database
    if type(feature) is np.ndarray:
        surfaces_to_database(surfaces)
    success = 'True'
    return success

def write_features_to_database(feature):
    """
    Write to database

    write_to_database()
    """
    from write_to_database.py import features_to_database
    if type(feature) is np.ndarray:
        features_to_database(feature)
    success = 'True'
    return success

def write_measures_to_database(measurements):
    """
    Write to database

    write_to_database()
    """
    from write_to_database.py import measure_to_database
    for measurement in measurements:
        if type(measurement) is np.ndarray:
            measure_to_database(measurement)
    success = 'True'
    return success
    
def write_measures_to_table(measurement):
    """
    Write to database

    write_to_database()
    """
    from write_to_table.py import measure_to_table
    for measurement in measurements:
        if type(measurement) is np.ndarray:
            measure_to_table(measurement)
    success = 'True'
    return success

