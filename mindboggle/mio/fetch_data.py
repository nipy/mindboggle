#!/usr/bin/env python
"""
Functions for fetching data from a URL or from third party software.

Authors:
    - Arno Klein, 2013-2015  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2015,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def hashes_url():
    """
    Hashes and URL to verify retrieved data used by the Mindboggle software.

    Data include atlases and templates required for registration and labeling.
    See fetch_hash() to create new hashes.

    Returns
    -------
    hashes : dictionary
        dictionary of data file names and hashes for those files
    url : string
        URL to Mindboggle data
    cache_env : string
        cache environment variable
    cache : string
        default path to cache directory (in case not supplied elsewhere)

    """
    import os

    hashes = {}
    url = 'http://media.mindboggle.info/data/cache/atlases/'
    cache_env = 'MINDBOGGLE_CACHE'
    cache = os.environ.get(cache_env, os.path.join(os.environ['HOME'],
                                                   'mindboggle_cache'))

    #-------------------------------------------------------------------------
    # Atlases and templates:
    #-------------------------------------------------------------------------
    hashes['OASIS-30_Atropos_template.nii.gz'] = 'f95dbe37ab40e8ad59c1b1eabc7f230c'
    hashes['OASIS-30_Atropos_template_in_MNI152.nii.gz'] = '2efbe9c0755f4aa4506ad78e60ec319f'
    hashes['OASIS-30_Atropos_template_to_MNI152_affine.txt'] = 'f36e3d5d99f7c4a9bb70e2494ed7340b'
    hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_OASIS-30_v2.nii.gz'] = 'abd9a4dec1017e0d3a69fd61d02e88bf'
    hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152_v2.nii.gz'] = 'ef2b6f7116e05b46746c044bfa995546'
    hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_label_probabilities_in_OASIS-30_v2.nii.gz'] = '70301f193e3cbb1386738c02064adc51'
    hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_label_probabilities_in_MNI152_v2.nii.gz'] = '2756e8b98f3600553bd515af2d3da808'
    hashes['MNI152_T1_1mm_brain.nii.gz'] = '7c47f203858f57c50ee93bcf7b0662fc'
    hashes['lh.DKTatlas40.gcs'] = 'b28a600a713c269f4c561f66f64337b2'
    hashes['rh.DKTatlas40.gcs'] = '50150021c84b0e829aed3a68709404d7'
    hashes['depth_curv_border_nonborder_parameters.pkl'] = 'a943a0f46c2c2b3bfdfdf5a64176c6c3'

    return hashes, url, cache_env, cache


def test_urls():
    """
    URLs corresponding to Mindboggle test (example output) data.

    Returns
    -------
    urls : dictionary
        dictionary of names and urls for data files

    Examples
    --------
    >>> from mindboggle.mio.fetch_data import test_urls
    >>> urls = test_urls()
    >>> urls['left_mean_curvature']
    'http://media.mindboggle.info/data/cache/ex/shapes/left_cortical_surface/mean_curvature.vtk'

    """
    url = 'http://media.mindboggle.info/data/cache/'
    ANTS = url + 'ants/'
    MBW = url + 'mindboggle_working/Mindboggle/'
    F = url + 'mindboggled/features/'
    L = url + 'mindboggled/labels/'
    S = url + 'mindboggled/shapes/'
    T = url + 'mindboggled/tables/'
    left = 'left_cortical_surface/'
    right = 'right_cortical_surface/'
    Fleft = F + left
    Fright = F + right
    Lleft = L + left
    Lright = L + right
    Sleft = S + left
    Sright = S + right
    Tleft = T + left
    Tright = T + right

    urls = {}
    #-------------------------------------------------------------------------
    # ants (antsCorticalThickness.sh) output:
    #-------------------------------------------------------------------------
    urls['ants_segmentation'] = ANTS + 'antsBrainSegmentation.nii.gz'
    #-------------------------------------------------------------------------
    # Mindboggle working directory (including converted FreeSurfer output):
    #-------------------------------------------------------------------------
    urls['freesurfer_segmentation'] = MBW + 'mgh_to_nifti/001.mgz.nii.gz'
    #-------------------------------------------------------------------------
    # Mindboggle features:
    #-------------------------------------------------------------------------
    urls['right_cortex_in_mni'] = Fleft + 'cortex_in_MNI152_space.vtk'
    urls['left_folds'] = Fleft + 'folds.vtk'
    urls['left_fundus_per_sulcus'] = Fleft + 'fundus_per_sulcus.vtk'
    urls['left_sulci'] = Fleft + 'sulci.vtk'
    urls['right_cortex_in_mni'] = Fright + 'cortex_in_MNI152_space.vtk'
    urls['right_folds'] = Fright + 'folds.vtk'
    urls['right_fundus_per_sulcus'] = Fright + 'fundus_per_sulcus.vtk'
    urls['right_sulci'] = Fright + 'sulci.vtk'
    #-------------------------------------------------------------------------
    # Mindboggle labels:
    #-------------------------------------------------------------------------
    urls['ants_labels'] = L + 'ants_filled_labels.nii.gz'
    urls['freesurfer_labels'] = L + 'freesurfer_wmparc_filled_labels.nii.gz'
    urls['left_freesurfer_labels'] = Lleft + 'freesurfer_cortex_labels.vtk'
    urls['right_freesurfer_labels'] = Lright + 'freesurfer_cortex_labels.vtk'
    #-------------------------------------------------------------------------
    # Mindboggle shapes:
    #-------------------------------------------------------------------------
    urls['left_area'] = Sleft + 'area.vtk'
    urls['left_freesurfer_curvature'] = Sleft + 'freesurfer_curvature.vtk'
    urls['left_freesurfer_sulc'] = Sleft + 'freesurfer_sulc.vtk'
    urls['left_freesurfer_thickness'] = Sleft + 'freesurfer_thickness.vtk'
    urls['left_geodesic_depth'] = Sleft + 'geodesic_depth.vtk'
    urls['left_mean_curvature'] = Sleft + 'mean_curvature.vtk'
    urls['left_travel_depth'] = Sleft + 'travel_depth.vtk'
    urls['right_area'] = Sright + 'area.vtk'
    urls['right_freesurfer_curvature'] = Sright + 'freesurfer_curvature.vtk'
    urls['right_freesurfer_sulc'] = Sright + 'freesurfer_sulc.vtk'
    urls['right_freesurfer_thickness'] = Sright + 'freesurfer_thickness.vtk'
    urls['right_geodesic_depth'] = Sright + 'geodesic_depth.vtk'
    urls['right_mean_curvature'] = Sright + 'mean_curvature.vtk'
    urls['right_travel_depth'] = Sright + 'travel_depth.vtk'
    #-------------------------------------------------------------------------
    # Mindboggle tables:
    #-------------------------------------------------------------------------
    urls['thickinthehead_ants_labels_table'] = \
        T + 'thickinthehead_per_ants_cortex_label.csv'
    urls['thickinthehead_freesurfer_labels_table'] = \
        T + 'thickinthehead_per_freesurfer_cortex_label.csv'
    urls['volume_ants_labels_table'] = \
        T + 'volume_for_each_ants_label.csv'
    urls['volume_freesurfer_labels_table'] = \
        T + 'volume_for_each_freesurfer_label.csv'
    urls['left_fundus_shapes_table'] = Tleft + 'fundus_shapes.csv'
    urls['left_label_shapes_table'] = Tleft + 'label_shapes.csv'
    urls['left_sulcus_shapes_table'] = Tleft + 'sulcus_shapes.csv'
    urls['left_vertices_table'] = Tleft + 'vertices.csv'
    urls['right_fundus_shapes_table'] = Tright + 'fundus_shapes.csv'
    urls['right_label_shapes_table'] = Tright + 'label_shapes.csv'
    urls['right_sulcus_shapes_table'] = Tright + 'sulcus_shapes.csv'
    urls['right_vertices_table'] = Tright + 'vertices.csv'

    # Prepend with url:
    #for key in urls:
    #    urls[key] = url + urls[key]

    return urls


def prep_tests():
    """
    Prepare to fetch data in docstring tests.

    Returns
    -------
    urls : dictionary
        dictionary of names and urls for data files
    fetch_data : function
        fetch_data() function

    Examples
    --------
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> urls['left_mean_curvature']
    'http://media.mindboggle.info/data/cache/ex/shapes/left_cortical_surface/mean_curvature.vtk'
    >>> fetch_data
    <function mindboggle.mio.fetch_data.fetch_data>

    """
    from mindboggle.mio.fetch_data import fetch_data
    from mindboggle.mio.fetch_data import test_urls

    urls = test_urls()

    return urls, fetch_data


def fetch_hash(data_file):
    """
    Get hash of data file.

    Parameters
    ----------
    data_file : string
        data file name

    Returns
    -------
    hash : string
        hash of data file

    Examples
    --------
    >>> import os
    >>> from mindboggle.mio.fetch_data import hashes_url
    >>> from mindboggle.mio.fetch_data import fetch_hash
    >>> cache_env = 'MINDBOGGLE_CACHE'
    >>> hashes, url, cache_env, cache = hashes_url()
    >>> data_file = hashes.keys()[0]
    >>> data_path = os.path.join(os.environ[cache_env], data_file)
    >>> fetch_hash(data_path)

    """
    import hashlib

    # Compute the file's hash:
    hash = hashlib.md5(open(data_file, 'rb').read()).hexdigest()

    return hash


def fetch_data(url, output_file=''):
    """
    Get data from URL.

    Parameters
    ----------
    url : string
        URL for data file
    output_file : string
        name of output file (full path)

    Returns
    -------
    output_file : string
        name of output file (full path)

    Examples
    --------
    >>> from mindboggle.mio.fetch_data import fetch_data
    >>> from mindboggle.mio.fetch_data import hashes_url
    >>> hashes, url = hashes_url()
    >>> data_file = hashes.keys()[0]
    >>> output_file = 'test_output.nii.gz'
    >>> fetch_data(url, output_file)

    """
    import urllib

    # Download file to specified output:
    if output_file:
        output_file, foo = urllib.urlretrieve(url, output_file)
    # Download file as a temporary file:
    else:
        output_file, foo = urllib.urlretrieve(url)

    return output_file


def fetch_check_data(data_file, url='', hashes={}, cache_env='', cache='',
                  return_missing=False, lookup=True):
    """
    Get data file through a URL call and check its hash.

    Steps ::
        If hashes provided:
            1. Check hash table for data file.
            2. Check hash subdirectory within cache directory for data file.
            3. If data file not in cache, download file, compute hash,
               and verify hash.
            4. If hash correct, save file.
        Otherwise, simply download file or return file path as a string.

    Parameters
    ----------
    data_file : string
        data file name
    url : string
        URL for data file
    hashes : dictionary
        file names and md5 hashes (if empty, simply download file from url)
    cache_env : string
        environment variable name for cache path
    cache : string
        in case cache_env is not set, use as cache directory
    return_missing : Boolean
        if data_file not in hash, simply download data_file and return path
    lookup : Boolean
        Simply return data_file path

    Returns
    -------
    data_path : string
        data file name (full path)

    Examples
    --------
    >>> from mindboggle.mio.fetch_data import hashes_url
    >>> from mindboggle.mio.fetch_data import fetch_check_data
    >>> hashes, url, cache_env, cache = hashes_url()
    >>> data_file = hashes.keys()[0]
    >>> fetch_check_data(data_file, url, hashes, cache_env, cache)

    """
    import os
    import sys
    import shutil

    from mindboggle.mio.fetch_data import fetch_data, fetch_hash

    if lookup:

        #---------------------------------------------------------------------
        # If hashes provided, go through steps to check/download file:
        #---------------------------------------------------------------------
        if hashes:
    
            if not cache_env:
                cache_env = 'MINDBOGGLE_CACHE'
            if not cache:
                cache = os.path.join(os.environ['HOME'], 'hash_temp')
    
            #-----------------------------------------------------------------
            # Check hash table for file:
            #-----------------------------------------------------------------
            if data_file not in hashes.keys():
                if return_missing:
                    data_path = data_file
                    print("Retrieved file not in hashes: {0}".
                          format(data_path))
                    return data_path
                else:
                    sys.exit("Data file '{0}' not in hash table.".
                    format(data_file))
            else:
                stored_hash = hashes[data_file]
    
                #-------------------------------------------------------------
                # Create missing cache and hash directories:
                #-------------------------------------------------------------
                if cache_env in os.environ.keys():
                    cache = os.environ[cache_env]
                if not os.path.exists(cache):
                    print("Create missing cache directory: {0}".format(cache))
                    os.mkdir(cache)
                hash_dir = os.path.join(cache, stored_hash)
                if not os.path.exists(hash_dir):
                    print("Create missing hash directory: {0}".
                          format(hash_dir))
                    os.mkdir(os.path.join(hash_dir))
        
                #-------------------------------------------------------------
                # Check hash subdirectory for file:
                #-------------------------------------------------------------
                data_path = os.path.join(hash_dir, data_file)
                if os.path.exists(data_path):
                    return data_path
        
                #-------------------------------------------------------------
                # If file not in cache, download, compute hash, and verify:
                #-------------------------------------------------------------
                else:
                    print("Retrieve file from the Mindboggle website: {0}".
                          format(url+data_file))
        
                    # Download file as a temporary file:
                    temp_file = fetch_data(url+data_file)
        
                    # Compute the file's hash:
                    data_hash = fetch_hash(temp_file)
        
                    # If hash matches name of the hash directory, save file:
                    if os.path.join(cache, data_hash) == hash_dir:
                        print("Copy file to cache: {0}".format(data_path))
                        shutil.copyfile(temp_file, data_path)
                        return data_path
                    else:
                        print("Retrieved hash does not match stored hash.")
    
        #---------------------------------------------------------------------
        # If hashes not provided, simply download file:
        #---------------------------------------------------------------------
        elif url:
            # Download file as a temporary file:
            data_path = fetch_data(url+data_file)
            print("Hashes not provided. Retrieved file: {0}".format(data_path))
            return data_path
    
        #---------------------------------------------------------------------
        # If URL also not provided, simply return file path:
        #---------------------------------------------------------------------
        else:
            data_path = data_file
            print("Neither hashes nor URL provided. "
                  "Returning file path: {0}".format(data_path))
            return data_path

    #-------------------------------------------------------------------------
    # Simply return file path:
    #-------------------------------------------------------------------------
    else:
        data_path = data_file
        print("Returning file path: {0}".format(data_path))
        return data_path


def fetch_ants_data(segmented_file, use_ants_transforms=True):
    """
    Fetch antsCorticalThickness.sh output.

    The input argument "segmented_file" is one of the relevant
    antsCorticalThickness.sh output files called by Mindboggle
    (assume path and PREFIX="ants"):
        ants_subjects/subject1/antsBrainExtractionMask.nii.gz
        ants_subjects/subject1/antsBrainSegmentation.nii.gz
        ants_subjects/subject1/antsSubjectToTemplate0GenericAffine.mat
        ants_subjects/subject1/antsSubjectToTemplate1Warp.nii.gz
        ants_subjects/subject1/antsTemplateToSubject0Warp.nii.gz
        ants_subjects/subject1/antsTemplateToSubject1GenericAffine.mat
    The existence of the transform files are checked only if
    use_ants_transforms == True. Transforms can only be included
    if they have been generated by antsCorticalThickness.sh
    with the -k argument.

    Parameters
    ----------
    segmented_file : string
        full path to a subject's antsCorticalThickness.sh segmented file
    use_ants_transforms : Boolean
        include antsCorticalThickness.sh-generated transforms?

    Returns
    -------
    mask : string
        antsBrainExtraction.sh brain volume mask for extracting brain volume
    segments : string
        Atropos-segmented brain volume
    affine_subject2template : string
        subject to template affine transform (antsRegistration)
    warp_subject2template : string
        subject to template nonlinear transform (antsRegistration)
    affine_template2subject : string
        template to subject affine transform (antsRegistration)
    warp_template2subject : string
        template to subject nonlinear transform (antsRegistration)

    Examples
    --------
    >>> from mindboggle.mio.fetch_data import fetch_ants_data
    >>> segmented_file = 'ants_subjects/OASIS-TRT-20-1/tmpBrainSegmentation.nii.gz'
    >>> fetch_ants_data(segmented_file)

    """
    import os

    prefix = segmented_file.strip('BrainSegmentation.nii.gz')

    mask = prefix + 'BrainExtractionMask.nii.gz'
    segments = segmented_file

    if use_ants_transforms:
        affine_subject2template = prefix + 'SubjectToTemplate0GenericAffine.mat'
        warp_subject2template = prefix + 'SubjectToTemplate1Warp.nii.gz'
        affine_template2subject = prefix + 'TemplateToSubject1GenericAffine.mat'
        warp_template2subject = prefix + 'TemplateToSubject0Warp.nii.gz'
        files = [mask, segments,
                 affine_subject2template, warp_subject2template,
                 affine_template2subject, warp_template2subject]
    else:
        affine_subject2template = ''
        warp_subject2template = ''
        affine_template2subject = ''
        warp_template2subject = ''
        files = [mask, segments]

    # The existence of the transform files are checked only if
    # use_ants_transforms == True. Transforms are generated by
    # antsCorticalThickness.sh when the -k argument is used.
    for s in files:
        if not os.path.exists(s):
            str1 = 'antsCorticalThickness.sh output ' + s + ' does not exist.'
            raise(IOError(str1))

    return mask, segments, affine_subject2template, warp_subject2template, \
                           affine_template2subject, warp_template2subject


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod()