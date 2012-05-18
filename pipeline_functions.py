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

##############################################################################
#   Surface map calculation
##############################################################################

def convert_to_vtk(fs_surface_file):
    """Measure

    measure_()
    """
    import subprocess as sp
    surface_file = fs_surface_file + '.vtk'
    cmd = ['mris_convert', fs_surface_file, vtk_file]
    proc = sp.Popen(' '.join(cmd))
    o, e = proc.communicate()
    if proc.returncode > 0 :
        raise Exception('\n'.join([cmd + ' failed', o, e]))
    return surface_file

def measure_surface_maps(surface_file):
    """Measure

    measure_()
    """
    from measure.py import measure_surface_maps
    if type(surface_file) is np.ndarray:
        depth_map, mean_curvature_map, gauss_curvature_map = measure_surface_maps(surface_file)
    return depth_map, mean_curvature_map, gauss_curvature_map

##############################################################################
#   Feature extraction
##############################################################################

def extract_sulci(surface_file, depth_map, mean_curvature_map, gauss_curvature_map):
    """Extract sulci

    extract_sulci
    """
    import subprocess as sp
    cmd = 'feature = extract/sulci/extract.py'
    cmd = ['python', cmd, '%s'%surface_file, '%s'%depth_map]
    proc = sp.Popen(cmd)
    o, e = proc.communicate()
    if proc.returncode > 0 :
        raise Exception('\n'.join(['extract.py failed', o, e]))
    #output_file = glob('file1.vtk').pop()
    #feature_files = glob('*.vtk')
    #return feature_files
    return sulci

def extract_fundi(surface_file, depth_map, mean_curvature_map, gauss_curvature_map):
    """Extract fundi

    extract_fundi
    """
    import subprocess as sp
    cmd = 'feature = extract/fundi/extract.py'
    cmd = ['python', cmd, '%s'%surface_file, '%s'%depth_map]
    proc = sp.Popen(cmd)
    o, e = proc.communicate()
    if proc.returncode > 0 :
        raise Exception('\n'.join(['extract.py failed', o, e]))
    return fundi

def extract_pits(surface_file, depth_map, mean_curvature_map, gauss_curvature_map):
    """Extract pits

    extract_pits
    """
    from glob import glob
    import subprocess as sp
    cmd = 'feature = extract/pits/extract.py'
    cmd = ['python', cmd, '%s'%surface_file, '%s'%depth_map]
    proc = sp.Popen(cmd)
    o, e = proc.communicate()
    if proc.returncode > 0 :
        raise Exception('\n'.join(['extract.py failed', o, e]))
    return pits

def extract_midaxis(surface_file, depth_map, mean_curvature_map, gauss_curvature_map):
    """Extract midaxis

    extract_midaxis
    """
    from glob import glob
    import subprocess as sp
    cmd = 'feature = extract/midaxis/extract.py'
    cmd = ['python', cmd, '%s'%surface_file, '%s'%depth_map]
    proc = sp.Popen(cmd)
    o, e = proc.communicate()
    if proc.returncode > 0 :
        raise Exception('\n'.join(['extract.py failed', o, e]))
    return midaxis

# Labeled surface patch and volume extraction nodes
def extract_patches(labels):
    """Extract labeled surface patches
    
    extract_patches
    """
    from glob import glob
    import subprocess as sp
    cmd = 'patch = extract/labels/extract.py'
    cmd = ['python', cmd]
    proc = sp.Popen(cmd)
    o, e = proc.communicate()
    if proc.returncode > 0 :
        raise Exception('\n'.join(['extract.py failed', o, e]))
    return patches

def extract_regions(labels):
    """Extract labeled region volumes
    
    extract_regions
    """
    from glob import glob
    import subprocess as sp
    cmd = 'region = extract/labels/extract.py'
    cmd = ['python', cmd]
    proc = sp.Popen(cmd)
    o, e = proc.communicate()
    if proc.returncode > 0 :
        raise Exception('\n'.join(['extract.py failed', o, e]))
    return regions

##############################################################################
#   Multi-atlas registration
##############################################################################

def register_template(subject_id, subjects_path, 
                      template_id, template_path, registration_name):
    """Register surface to template with FreeSurfer's mris_register

    Example: bert
             /Applications/freesurfer/subjects
             ./templates_freesurfer
             KKI_2.tif
             sphere_to_template.reg
    """
    import os

    for hemi in ['lh','rh']:
        input_file = os.path.join(subjects_path, subject_id, 'surf', hemi + '.sphere')
        output_file = os.path.join(subjects_path, subject_id, 'surf', hemi + '.' + registration_name)
        template_file = os.path.join(template_path, hemi + '.' + template_id)
        args = ['mris_register -curv', input_file, template_file, output_file]
        print(' '.join(args));
        #os.system(' '.join(args)); # p = Popen(args);
        proc = sp.Popen(args)
        o, e = proc.communicate()
        if proc.returncode > 0 :
            raise Exception('\n'.join([cmd + ' failed', o, e]))
        return output_name

def register_atlases(subject_id, atlas_list_file, 
                     registration_name, output_path):
    """Transform the labels from multiple atlases via a template
    using FreeSurfer's mri_surf2surf (wrapped in NiPype)

    nipype.workflows.smri.freesurfer.utils.fs.SurfaceTransform
    wraps command **mri_surf2surf**:
    "Transform a surface file from one subject to another via a spherical registration.
    Both the source and target subject must reside in your Subjects Directory,
    and they must have been processed with recon-all, unless you are transforming
    to one of the icosahedron meshes."
    """
    import os
    from nipype.interfaces.freesurfer import SurfaceTransform

    annot_file_name = 'aparcNMMjt.annot'

    sxfm = SurfaceTransform()
    sxfm.inputs.target_subject = subject_id

    # Get list of atlas subjects from a file
    f = open(atlas_list_file)
    atlas_list = f.readlines()
    for atlas_line in atlas_list:
        # For each atlas
        atlas_name = atlas_line.strip("\n")
        sxfm.inputs.source_subject = atlas_name
        # For each hemisphere
        for hemi in ['lh','rh']:        
            sxfm.inputs.hemi = hemi
            # Specify annotation file
            output_annot = os.path.join(output_path, hemi + '.' + atlas_name + '_to_' + \
                                        subject_id + '_' + annot_file_name)
            args = ['--sval-annot', annot_file_name,
                    '--tval', output_annot,
                    '--srcsurfreg', registration_name,
                    '--trgsurfreg', registration_name]
            sxfm.inputs.args = ' '.join(args)
            sxfm.run()
    return atlas_list

##############################################################################
#   Label propagation
##############################################################################

def propagate_labels(labels, fundi):
    """Propagate labels
    """
    return labels

# Volume label propagation node
def propagate_volume_labels(labels):
    """Propagate labels through volume
    """
    return labels


##############################################################################
#   Feature segmentation / identification
##############################################################################

def segment_sulci(sulci):
    """Segment and identify sulci
    """
    import os

    for hemi in ['lh','rh']:
        input_file = os.path.join(subject_surf_path, hemi + '.sphere')
        output_file = os.path.join(subject_surf_path, hemi + '.' + registration_name)
        template_file = os.path.join(template_path, hemi + '.' + template_id)
        args = ['mris_register -curv', input_file, template_file, output_file]
        print(' '.join(args));
        #os.system(' '.join(args)); # p = Popen(args);
        proc = sp.Popen(args)
        o, e = proc.communicate()
        if proc.returncode > 0 :
            raise Exception('\n'.join([cmd + ' failed', o, e]))
        return feature_type, segmented_sulci

def segment_fundi(fundi):
    """Segment and identify fundi
    """
    import os

    for hemi in ['lh','rh']:
        input_file = os.path.join(subject_surf_path, hemi + '.sphere')
        output_file = os.path.join(subject_surf_path, hemi + '.' + registration_name)
        template_file = os.path.join(template_path, hemi + '.' + template_id)
        args = ['mris_register -curv', input_file, template_file, output_file]
        print(' '.join(args));
        #os.system(' '.join(args)); # p = Popen(args);
        proc = sp.Popen(args)
        o, e = proc.communicate()
        if proc.returncode > 0 :
            raise Exception('\n'.join([cmd + ' failed', o, e]))
        return feature_type, segmented_fundi

def segment_midaxis(midaxis):
    """Segment and identify medial axis surfaces
    """
    import os

    for hemi in ['lh','rh']:
        input_file = os.path.join(subject_surf_path, hemi + '.sphere')
        output_file = os.path.join(subject_surf_path, hemi + '.' + registration_name)
        template_file = os.path.join(template_path, hemi + '.' + template_id)
        args = ['mris_register -curv', input_file, template_file, output_file]
        print(' '.join(args));
        #os.system(' '.join(args)); # p = Popen(args);
        proc = sp.Popen(args)
        o, e = proc.communicate()
        if proc.returncode > 0 :
            raise Exception('\n'.join([cmd + ' failed', o, e]))
        return feature_type, segmented_midaxis

##############################################################################
#   High-level shape measurement functions
##############################################################################

def measure_positions(features):
    """Measure

    measure_()
    """
    for feature in features:
        from measure.py import measure_position
        if type(feature) is np.ndarray:
            measurement = measure_position(feature)

def measure_extents(features):
    """Measure

    measure_()
    """
    for feature in features:
        from measure.py import measure_extent
        if type(feature) is np.ndarray:
            measurement = measure_extent(feature)

def measure_curvatures(features):
    """Measure

    measure_()
    """
    for feature in features:
        from measure.py import measure_curvature
        if type(feature) is np.ndarray:
            measurement = measure_curvature(feature)

def measure_depths(features):
    """Measure

    measure_()
    """
    for feature in features:
        from measure.py import measure_depth
        if type(feature) is np.ndarray:
            measurement = measure_depth(feature)

def measure_spectra(features):
    """Measure

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

def measure_spectrum(feature):
    """Measure

    measure_()
    """
    from measure.py import measure_spectrum
    if type(feature) is np.ndarray:
        measurement = measure_spectrum(feature)
    return measurement

##############################################################################
#    Store features in database
##############################################################################

def write_maps_to_database(surface_map):
    """Write to database

    write_to_database()
    """
    from write_to_database.py import maps_to_database
    if type(feature) is np.ndarray:
        maps_to_database(surface_map)
    success = 'True'
    return success

def write_features_to_database(feature):
    """Write to database

    write_to_database()
    """
    from write_to_database.py import features_to_database
    if type(feature) is np.ndarray:
        features_to_database(feature)
    success = 'True'
    return success

def write_measures_to_database(measurements):
    """Write to database

    write_to_database()
    """
    from write_to_database.py import measure_to_database
    for measurement in measurements:
        if type(measurement) is np.ndarray:
            measure_to_database(measurement)
    success = 'True'
    return success
    
def write_measures_to_table(measurement):
    """Write to database

    write_to_database()
    """
    from write_to_table.py import measure_to_table
    for measurement in measurements:
        if type(measurement) is np.ndarray:
            measure_to_table(measurement)
    success = 'True'
    return success

