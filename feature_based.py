#!/usr/bin/python

""" test
Feature-based functions for morphometry and labeling:

1. Surface calculations

2. Feature extraction

3. Label propagation

4. Feature segmentation / identification

5. Shape measurement


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

##############################################################################
#   Surface calculations
##############################################################################

def measure_surface_depth(command, surface_file):
    """
    Measure

    measure_()
    """
    from nipype.interfaces.base import CommandLine

    depth_file = surface_file.strip('.vtk') + '.depth.vtk'
    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join([surface_file, depth_file])
    cli.cmdline
    return depth_file
    
def measure_surface_curvature(command, surface_file):
    """
    Measure
    CurvatureMain input MeanCurvatureOutput [GaussianCurvatureOutput]
                        [MaximalCurvatureOutput] [MinimalCurvatureOutput]
    measure_()
    """
    from nipype.interfaces.base import CommandLine

    mean_curvature_file = surface_file.strip('.vtk') + '.curvature.mean.vtk'
    gauss_curvature_file = surface_file.strip('.vtk') + '.curvature.gauss.vtk'
    max_curvature_file = surface_file.strip('.vtk') + '.curvature.max.vtk'
    min_curvature_file = surface_file.strip('.vtk') + '.curvature.min.vtk'
    args = [surface_file, 
            mean_curvature_file, gauss_curvature_file,
            max_curvature_file, min_curvature_file]
    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join(args)
    cli.cmdline
    return mean_curvature_file, gauss_curvature_file,\
           max_curvature_file, min_curvature_file  
    
##############################################################################
#   Feature extraction
##############################################################################

def extract_fundi(extract_fundi_command, depth_files):
    """
    Extract fundi

    extract_fundi
    """
    from nipype.interfaces.base import CommandLine

    fundi = []
    for depth_file in depth_files:
        args = ['python', extract_fundi_command, '%s'%depth_file]
        cli = CommandLine(command = curvature_command)
        cli.inputs.args = ' '.join(args)
        cli.cmdline
    return fundi

def extract_sulci(surface_file, depth_file, mean_curvature_file, gauss_curvature_file):
    """
    Extract sulci

    extract_sulci    
    """
    from nipype.interfaces.base import CommandLine

    args = 'feature = extract/sulci/extract.py'
    args = ['python', args, '%s'%surface_file, '%s'%depth_file]
    print(' '.join(args)); os.system(' '.join(args))
    #output_file = glob('file1.vtk').pop()
    #feature_files = glob('*.vtk')
    #return feature_files
    return sulci

def extract_midaxis(surface_file, depth_file, mean_curvature_file, gauss_curvature_file):
    """
    Extract midaxis

    extract_midaxis    
    """
    from nipype.interfaces.base import CommandLine

    args = 'feature = extract/midaxis/extract.py'
    args = ['python', args, '%s'%surface_file, '%s'%depth_file]
    print(' '.join(args)); os.system(' '.join(args))
    return midaxis

# Labeled surface patch and volume extraction nodes
def extract_patches(labels):
    """
    Extract labeled surface patches
    
    extract_patches    
    """
    from nipype.interfaces.base import CommandLine

    args = ['python', 'patch = extract/labels/extract.py']
    print(' '.join(args)); os.system(' '.join(args))
    return patches

def extract_regions(labels):
    """
    Extract labeled region volumes
    
    extract_regions
    """
    from nipype.interfaces.base import CommandLine

    args = ['python', 'region = extract/labels/extract.py']
    print(' '.join(args)); os.system(' '.join(args))
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

    for hemi in ['lh','rh']:
        input_file = os.path.join(subject_surf_path, hemi + '.sphere')
        output_file = os.path.join(subject_surf_path, hemi + '.' + reg_name)
        template_file = os.path.join(templates_path, hemi + '.' + template_name)
        args = ['mris_register -curv', input_file, template_file, output_file]
        print(' '.join(args)); os.system(' '.join(args))
        return feature_type, segmented_sulci

def segment_fundi(fundi):
    """
    Segment and identify fundi
    """
    from nipype.interfaces.base import CommandLine

    for hemi in ['lh','rh']:
        input_file = os.path.join(subject_surf_path, hemi + '.sphere')
        output_file = os.path.join(subject_surf_path, hemi + '.' + reg_name)
        template_file = os.path.join(templates_path, hemi + '.' + template_name)
        args = ['mris_register -curv', input_file, template_file, output_file]
        print(' '.join(args)); os.system(' '.join(args))
        return feature_type, segmented_fundi

def segment_midaxis(midaxis):
    """
    Segment and identify medial axis surfaces
    """
    from nipype.interfaces.base import CommandLine

    for hemi in ['lh','rh']:
        input_file = os.path.join(subject_surf_path, hemi + '.sphere')
        output_file = os.path.join(subject_surf_path, hemi + '.' + reg_name)
        template_file = os.path.join(templates_path, hemi + '.' + template_name)
        args = ['mris_register -curv', input_file, template_file, output_file]
        print(' '.join(args)); os.system(' '.join(args))
        return feature_type, segmented_midaxis

##############################################################################
#   High-level shape measurement functions
##############################################################################

def measure_positions(features):
    """
    Measure

    measure_()
    """
    for feature in features:
        from measure.py import measure_position
        if type(feature) is np.ndarray:
            measurement = measure_position(feature)

def measure_extents(features):
    """
    Measure

    measure_()
    """
    for feature in features:
        from measure.py import measure_extent
        if type(feature) is np.ndarray:
            measurement = measure_extent(feature)

def measure_curvatures(features):
    """
    Measure

    measure_()
    """
    for feature in features:
        from measure.py import measure_curvature
        if type(feature) is np.ndarray:
            measurement = measure_curvature(feature)

def measure_depths(features):
    """
    Measure

    measure_()
    """
    for feature in features:
        from measure.py import measure_depth
        if type(feature) is np.ndarray:
            measurement = measure_depth(feature)

def measure_spectra(features):
    """
    Measure

    measure_()
    """
    for feature in features:
        from measure.py import measure_spectrum
        if type(feature) is np.ndarray:
            measurement = measure_spectrum(feature)

##############################################################################
#   Individual shape measurement functions
##############################################################################

def measure_position(feature):
    """
    Measure

    measure_()
    """
    from measure.py import measure_position
    if type(feature) is np.ndarray:
        measurement = measure_position(feature)
    return measurement

def measure_extent(feature):
    """
    Measure

    measure_()
    """
    from measure.py import measure_extent
    if type(feature) is np.ndarray:
        measurement = measure_extent(feature)
    return measurement

def measure_depth(feature):
    """
    Measure

    measure_()
    """
    from measure.py import measure_depth
    if type(feature) is np.ndarray:
        measurement = measure_(feature)
    return measurement

def measure_curvature(feature):
    """
    Measure

    measure_()
    """
    from measure.py import measure_curvature
    if type(feature) is np.ndarray:
        measurement = measure_curvature(feature)
    return measurement

def measure_spectrum(feature):
    """
    Measure

    measure_()
    """
    from measure.py import measure_spectrum
    if type(feature) is np.ndarray:
        measurement = measure_spectrum(feature)
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

