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

def extract_features(surface_file, depth_map, mean_curvature_map, gauss_curvature_map, feature_type):
    """Extract features

    extract_features
    """

    if feature_type == 'sulci':
        def extract_sulci(surface_file, depth_map, mean_curvature_map, gauss_curvature_map):
            """Extract sulci
        
            extract_sulci
            """
            from glob import glob
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

    elif feature_type == 'fundi':
        def extract_fundi(surface_file, depth_map, mean_curvature_map, gauss_curvature_map):
            """Extract fundi
        
            extract_fundi
            """
            from glob import glob
            import subprocess as sp
            cmd = 'feature = extract/fundi/extract.py'
            cmd = ['python', cmd, '%s'%surface_file, '%s'%depth_map]
            proc = sp.Popen(cmd)
            o, e = proc.communicate()
            if proc.returncode > 0 :
                raise Exception('\n'.join(['extract.py failed', o, e]))
            return fundi

    elif feature_type == 'pits':
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

    elif feature_type == 'medial':
        def extract_medial(surface_file, depth_map, mean_curvature_map, gauss_curvature_map):
            """Extract medial
        
            extract_medial
            """
            from glob import glob
            import subprocess as sp
            cmd = 'feature = extract/medial/extract.py'
            cmd = ['python', cmd, '%s'%surface_file, '%s'%depth_map]
            proc = sp.Popen(cmd)
            o, e = proc.communicate()
            if proc.returncode > 0 :
                raise Exception('\n'.join(['extract.py failed', o, e]))
            return medial

def register_template(subject_id, subject_surf_path, template_path, 
                         template_name, registration_name):
    """Register surface to template with FreeSurfer's mris_register

    Example: /Applications/freesurfer/subjects/bert/surf
             ./templates_freesurfer
             KKI_2.tif
             sphere_to_template.reg
    """
    import os

    for hemi in ['lh','rh']:
        input_file = os.path.join(subject_surf_path, hemi + '.sphere')
        output_file = os.path.join(subject_surf_path, hemi + '.' + registration_name)
        template_file = os.path.join(template_path, hemi + '.' + template_name)
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

def segment_features(feature_type, atlas_list):
    """Extract features

    extract_features
    """

    if feature_type == 'sulci':
        def segment_sulci(sulci):
            """Segment and identify sulci
            """
            import os
        
            for hemi in ['lh','rh']:
                input_file = os.path.join(subject_surf_path, hemi + '.sphere')
                output_file = os.path.join(subject_surf_path, hemi + '.' + registration_name)
                template_file = os.path.join(template_path, hemi + '.' + template_name)
                args = ['mris_register -curv', input_file, template_file, output_file]
                print(' '.join(args));
                #os.system(' '.join(args)); # p = Popen(args);
                proc = sp.Popen(args)
                o, e = proc.communicate()
                if proc.returncode > 0 :
                    raise Exception('\n'.join([cmd + ' failed', o, e]))
                return segmented_sulci
        
    elif feature_type == 'fundi':
        def segment_fundi(fundi):
            """Segment and identify fundi
            """
            import os
        
            for hemi in ['lh','rh']:
                input_file = os.path.join(subject_surf_path, hemi + '.sphere')
                output_file = os.path.join(subject_surf_path, hemi + '.' + registration_name)
                template_file = os.path.join(template_path, hemi + '.' + template_name)
                args = ['mris_register -curv', input_file, template_file, output_file]
                print(' '.join(args));
                #os.system(' '.join(args)); # p = Popen(args);
                proc = sp.Popen(args)
                o, e = proc.communicate()
                if proc.returncode > 0 :
                    raise Exception('\n'.join([cmd + ' failed', o, e]))
                return segmented_fundi

    elif feature_type == 'medial':
        def segment_medial(medial):
            """Segment and identify medial surfaces
            """
            import os
        
            for hemi in ['lh','rh']:
                input_file = os.path.join(subject_surf_path, hemi + '.sphere')
                output_file = os.path.join(subject_surf_path, hemi + '.' + registration_name)
                template_file = os.path.join(template_path, hemi + '.' + template_name)
                args = ['mris_register -curv', input_file, template_file, output_file]
                print(' '.join(args));
                #os.system(' '.join(args)); # p = Popen(args);
                proc = sp.Popen(args)
                o, e = proc.communicate()
                if proc.returncode > 0 :
                    raise Exception('\n'.join([cmd + ' failed', o, e]))
                return segmented_medial

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

def measure_spectra(feature):
    """Measure

    measure_()
    """
    from measure.py import measure_spectra
    if type(feature) is np.ndarray:
        measurement = measure_spectra(feature)
    return measurement

def write_to_database(feature):
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
