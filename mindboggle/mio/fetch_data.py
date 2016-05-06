#!/usr/bin/env python
"""
Functions for fetching data from a URL or from third party software.

Authors:
    - Arno Klein, 2013-2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def cache_hashes():
    """
    Hashes to verify retrieved data cached by the Mindboggle software.

    Data include atlases and templates required for registration and labeling.
    See fetch_hash() to create new hashes.

    Returns
    -------
    hashes : dictionary
        dictionary of data file names and hashes for those files

    """

    hashes = {}

    # ------------------------------------------------------------------------
    # Atlases and templates:
    # ------------------------------------------------------------------------
    ## atlas_volume and atropos_to_MNI152_affine:
    hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152_v2.nii.gz'] = 'ef2b6f7116e05b46746c044bfa995546'
    hashes['OASIS-30_Atropos_template_to_MNI152_affine.txt'] = 'f36e3d5d99f7c4a9bb70e2494ed7340b'
    ## Likelihood parameters:
    hashes['depth_curv_border_nonborder_parameters.pkl'] = 'a943a0f46c2c2b3bfdfdf5a64176c6c3'

    # hashes['OASIS-30_Atropos_template.nii.gz'] = 'f95dbe37ab40e8ad59c1b1eabc7f230c'
    # hashes['OASIS-30_Atropos_template_in_MNI152.nii.gz'] = '2efbe9c0755f4aa4506ad78e60ec319f'
    # hashes['OASIS-30_Atropos_template_to_MNI152_affine.txt'] = 'f36e3d5d99f7c4a9bb70e2494ed7340b'
    # hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_OASIS-30_v2.nii.gz'] = 'abd9a4dec1017e0d3a69fd61d02e88bf'
    # hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152_v2.nii.gz'] = 'ef2b6f7116e05b46746c044bfa995546'
    # hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_label_probabilities_in_OASIS-30_v2.nii.gz'] = '70301f193e3cbb1386738c02064adc51'
    # hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_label_probabilities_in_MNI152_v2.nii.gz'] = '2756e8b98f3600553bd515af2d3da808'
    # hashes['MNI152_T1_1mm_brain.nii.gz'] = '7c47f203858f57c50ee93bcf7b0662fc'
    # hashes['lh.DKTatlas40.gcs'] = 'b28a600a713c269f4c561f66f64337b2'
    # hashes['rh.DKTatlas40.gcs'] = '50150021c84b0e829aed3a68709404d7'
    # hashes['depth_curv_border_nonborder_parameters.pkl'] = 'a943a0f46c2c2b3bfdfdf5a64176c6c3'

    return hashes


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
    >>> urls['left_mean_curvature']  # shapes/left_cortical_surface/mean_curvature.vtk
    'https://osf.io/2q7hb/?action=download&version=1'

    """

    urls = {}
    # ------------------------------------------------------------------------
    # Atlases and templates:
    # ------------------------------------------------------------------------
    # OASIS-30_Atropos_template.nii.gz
    urls['OASIS-30_Atropos_template'] = \
        'https://osf.io/ktxqy/?action=download&version=1'
    # OASIS-30_Atropos_template_in_MNI152.nii.gz
    urls['OASIS-30_Atropos_template_in_MNI152'] = \
        'https://osf.io/mgdyw/?action=download&version=1'
    # OASIS-30_Atropos_template_to_MNI152_affine.txt
    urls['OASIS-30_Atropos_template_to_MNI152_affine'] = \
        'https://osf.io/eyvzm/?action=download&version=1'
    # OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152_v2.nii.gz
    urls['OASIS-TRT-20_jointfusion_labels_in_MNI152'] = \
        'https://osf.io/5q46g/?action=download&version=1'
    # OASIS-TRT-20_jointfusion_DKT31_CMA_label_probabilities_in_MNI152_v2.nii.gz
    urls['OASIS-TRT-20_jointfusion_probabilities_in_MNI152'] = \
        'https://osf.io/dj6px/?action=download&version=1'
    # MNI152_T1_1mm_brain.nii.gz
    urls['MNI152_T1_brain'] = \
        'https://osf.io/tuvx7/?action=download&version=1'
    # depth_curv_border_nonborder_parameters.pkl
    urls['depth_curv_border_nonborder_parameters'] = \
        'https://osf.io/v9rqe/?action=download&version=1'
    # ------------------------------------------------------------------------
    # Manual labels:
    # ------------------------------------------------------------------------
    # lh.labels.DKT31.manual.vtk
    urls['left_manual_labels'] = \
        'https://osf.io/3rmzs/?action=download&version=1'
    # rh.labels.DKT31.manual.vtk
    urls['right_manual_labels'] = \
        'https://osf.io/y2pzr/?action=download&version=1'
    # ------------------------------------------------------------------------
    # FreeSurfer output:
    # ------------------------------------------------------------------------
    # label/lh.aparc.annot
    urls['left_freesurfer_aparc_annot'] = \
        'https://osf.io/awcbk/?action=download&version=1'
    # label/rh.aparc.annot
    urls['right_freesurfer_aparc_annot'] = \
        'https://osf.io/zhwge/?action=download&version=1'
    # mri/orig.mgz
    urls['freesurfer_orig_mgz'] = \
        'https://osf.io/4jcq3/?action=download&version=1'
    # mri/orig/001.mgz
    urls['freesurfer_001_mgz'] = \
        'https://osf.io/d26ja/?action=download&version=1'
    # surf/lh.pial
    urls['left_freesurfer_pial'] = \
        'https://osf.io/pgqms/?action=download&version=1'
    # surf/rh.pial
    urls['right_freesurfer_pial'] = \
        'https://osf.io/xn8ae/?action=download&version=1'
    # surf/lh.thickness
    urls['left_freesurfer_thickness'] = \
        'https://osf.io/rw6nv/?action=download&version=1'
    # surf/rh.thickness
    urls['right_freesurfer_thickness'] = \
        'https://osf.io/fb6nm/?action=download&version=1'
    # mri/t1weighted_brain.MNI152Affine.txt
    urls['affine_mni_transform'] = \
        'https://osf.io/ez7qx/?action=download&version=1'
    # ------------------------------------------------------------------------
    # ants (antsCorticalThickness.sh) output:
    # ------------------------------------------------------------------------
    # antsBrainSegmentation.nii.gz
    urls['ants_segmentation'] = \
        'https://osf.io/tbw95/?action=download&version=1'
    # antsBrainExtractionMask.nii.gz
    urls['ants_mask'] = \
        'https://osf.io/ar35t/?action=download&version=1'
    # antsSubjectToTemplate0GenericAffine.mat
    urls['ants_affine_subject2template'] = \
        'https://osf.io/t6pkv/?action=download&version=1'
    # antsSubjectToTemplate1Warp.nii.gz
    urls['ants_warp_subject2template'] = \
        'https://osf.io/qrhbz/?action=download&version=1'
    # antsTemplateToSubject0Warp.nii.gz
    urls['ants_warp_template2subject'] = \
        'https://osf.io/w5fzx/?action=download&version=1'
    # antsTemplateToSubject1GenericAffine.mat
    urls['ants_affine_template2subject'] = \
        'https://osf.io/n4puk/?action=download&version=1'
    # ------------------------------------------------------------------------
    # Mindboggle working directory (including converted FreeSurfer output):
    # ------------------------------------------------------------------------
    # Volume_labels/Freesurfer_cerebrum_labels_to_graywhite/wmparc.mgz.nii.gz
    urls['freesurfer_segmentation'] = \
        'https://osf.io/gs45h/?action=download&version=1'
    # Convert_MRI_to_nifti_format/001.mgz.nii.gz
    urls['T1_001'] = \
        'https://osf.io/m8vx6/?action=download&version=1'
    # _hemi_lh/Surface_to_vtk/lh.pial.vtk
    urls['left_pial'] = \
        'https://osf.io/fhs53/?action=download&version=1'
    # _hemi_rh/Surface_to_vtk/rh.pial.vtk
    urls['right_pial'] = \
        'https://osf.io/c9j7z/?action=download&version=1'
    # ------------------------------------------------------------------------
    # Mindboggle features:
    # ------------------------------------------------------------------------
    # cortex_in_MNI152_space.vtk
    urls['left_cortex_in_mni'] = \
        'https://osf.io/gv24u/?action=download&version=1'
    # folds.vtk
    urls['left_folds'] = \
        'https://osf.io/z793f/?action=download&version=1'
    # fundus_per_fold.vtk
    urls['left_fundus_per_fold'] = \
        'https://osf.io/syk4t/?action=download&version=1'
    # fundus_per_sulcus.vtk
    urls['left_fundus_per_sulcus'] = \
        'https://osf.io/te5aw/?action=download&version=1'
    # sulci.vtk
    urls['left_sulci'] = \
        'https://osf.io/2839a/?action=download&version=1'
    # cortex_in_MNI152_space.vtk
    urls['right_cortex_in_mni'] = \
        'https://osf.io/48ncx/?action=download&version=1'
    # folds.vtk
    urls['right_folds'] = \
        'https://osf.io/jy7bh/?action=download&version=1'
    # fundus_per_fold.vtk
    urls['right_fundus_per_fold'] = \
        'https://osf.io/7p3m8/?action=download&version=1'
    # fundus_per_sulcus.vtk
    urls['right_fundi'] = \
        'https://osf.io/2ajub/?action=download&version=1'
    # sulci.vtk
    urls['right_sulci'] = \
        'https://osf.io/sy6r8/?action=download&version=1'
    # ------------------------------------------------------------------------
    # Mindboggle labels:
    # ------------------------------------------------------------------------
    # ants_labels_in_hybrid_graywhite.nii.gz
    urls['ants_labels'] = \
        'https://osf.io/9jqg4/?action=download&version=1'
    # freesurfer_wmparc_labels_in_hybrid_graywhite.nii.gz
    urls['freesurfer_labels'] = \
        'https://osf.io/j478b/?action=download&version=1'
    # freesurfer_cortex_labels.vtk
    urls['left_freesurfer_labels'] = \
        'https://osf.io/fj47k/?action=download&version=1'
    # freesurfer_cortex_labels.vtk
    urls['right_freesurfer_labels'] = \
        'https://osf.io/ktcdg/?action=download&version=1'
    # ------------------------------------------------------------------------
    # Mindboggle shapes:
    # ------------------------------------------------------------------------
    # area.vtk
    urls['left_area'] = \
        'https://osf.io/e56rg/?action=download&version=1'
    # freesurfer_curvature.vtk
    urls['left_freesurfer_curvature'] = \
        'https://osf.io/37dxj/?action=download&version=1'
    # freesurfer_sulc.vtk
    urls['left_freesurfer_sulc'] = \
        'https://osf.io/fhabt/?action=download&version=1'
    # freesurfer_thickness.vtk
    urls['left_freesurfer_thickness'] = \
        'https://osf.io/pwj9u/?action=download&version=1'
    # geodesic_depth.vtk
    urls['left_geodesic_depth'] = \
        'https://osf.io/ru6hp/?action=download&version=1'
    # mean_curvature.vtk
    urls['left_mean_curvature'] = \
        'https://osf.io/2q7hb/?action=download&version=1'
    # travel_depth.vtk
    urls['left_travel_depth'] = \
        'https://osf.io/uzghr/?action=download&version=1'
    # area.vtk
    urls['right_area'] = \
        'https://osf.io/tpxj8/?action=download&version=1'
    # freesurfer_curvature.vtk
    urls['right_freesurfer_curvature'] = \
        'https://osf.io/2g95u/?action=download&version=1'
    # freesurfer_sulc.vtk
    urls['right_freesurfer_sulc'] = \
        'https://osf.io/vn9em/?action=download&version=1'
    # freesurfer_thickness.vtk
    urls['right_freesurfer_thickness'] = \
        'https://osf.io/4me2n/?action=download&version=1'
    # geodesic_depth.vtk
    urls['right_geodesic_depth'] = \
        'https://osf.io/4vc3d/?action=download&version=1'
    # mean_curvature.vtk
    urls['right_mean_curvature'] = \
        'https://osf.io/haqj3/?action=download&version=1'
    # travel_depth.vtk
    urls['right_travel_depth'] = \
        'https://osf.io/da5e2/?action=download&version=1'
    # ------------------------------------------------------------------------
    # Mindboggle tables:
    # ------------------------------------------------------------------------
    # thickinthehead_per_ants_cortex_label.csv
    urls['thickinthehead_ants_labels_table'] = \
        'https://osf.io/ujvrd/?action=download&version=1'
    # thickinthehead_per_freesurfer_cortex_label.csv
    urls['thickinthehead_freesurfer_labels_table'] = \
        'https://osf.io/976ew/?action=download&version=1'
    # volume_per_ants_label.csv
    urls['volume_ants_labels_table'] = \
        'https://osf.io/y6d2c/?action=download&version=1'
    # volume_per_freesurfer_label.csv
    urls['volume_freesurfer_labels_table'] = \
        'https://osf.io/ukqdy/?action=download&version=1'
    # fundus_shapes.csv
    urls['left_fundus_shapes_table'] = \
        'https://osf.io/23rqv/?action=download&version=1'
    # label_shapes.csv
    urls['left_label_shapes_table'] = \
        'https://osf.io/jecq3/?action=download&version=1'
    # sulcus_shapes.csv
    urls['left_sulcus_shapes_table'] = \
        'https://osf.io/tmhrp/?action=download&version=1'
    # vertices.csv
    urls['left_vertices_table'] = \
        'https://osf.io/6dcx7/?action=download&version=1'
    # fundus_shapes.csv
    urls['right_fundus_shapes_table'] = \
        'https://osf.io/4trm9/?action=download&version=1'
    # label_shapes.csv
    urls['right_label_shapes_table'] = \
        'https://osf.io/5sgn7/?action=download&version=1'
    # sulcus_shapes.csv
    urls['right_sulcus_shapes_table'] = \
        'https://osf.io/mndya/?action=download&version=1'
    # vertices.csv
    urls['right_vertices_table'] = \
        'https://osf.io/f6h72/?action=download&version=1'

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
    >>> urls['left_mean_curvature']  # shapes/left_cortical_surface/mean_curvature.vtk
    'https://osf.io/2q7hb/?action=download&version=1'

    """
    from mindboggle.mio.fetch_data import fetch_data, test_urls

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
    >>> from mindboggle.mio.fetch_data import fetch_hash, fetch_data
    >>> # osf.io URL for OASIS-30_Atropos_template_to_MNI152_affine.txt
    >>> url = 'https://osf.io/ufydw/?action=download&version=1'
    >>> data_file = fetch_data(url)
    >>> fetch_hash(data_file)
    'f36e3d5d99f7c4a9bb70e2494ed7340b'

    """
    from io import open

    import hashlib

    # Compute the file's hash:
    hash = hashlib.md5(open(data_file, 'rb').read()).hexdigest()

    return hash


def fetch_data(url, output_file='', append=''):
    """
    Download file from a URL to a specified or a temporary file.

    Optionally append to file name.

    Parameters
    ----------
    url : string
        URL for data file
    output_file : string
        name of output file (full path)
    append : string
        append to output file (ex: '.nii.gz')

    Returns
    -------
    output_file : string
        name of output file (full path)

    Examples
    --------
    >>> from mindboggle.mio.fetch_data import fetch_data, fetch_hash
    >>> output_file = ''
    >>> append = ''
    >>> # osf.io URL for OASIS-30_Atropos_template_to_MNI152_affine.txt
    >>> url = 'https://osf.io/ufydw/?action=download&version=1'
    >>> output_file = fetch_data(url, output_file, append)
    >>> fetch_hash(output_file)
    'f36e3d5d99f7c4a9bb70e2494ed7340b'

    """
    import os
    import urllib.request

    output_file, foo = urllib.request.urlretrieve(url, output_file)

    # Add append if assigned:
    if append:
        os.rename(output_file, output_file + append)
        output_file += append

    return output_file


def fetch_check_data(data_file, url, hashes, cache_directory='', append='',
                     verbose=False):
    """
    Get data file through a URL call and check its hash:

        1. Check hash table for data file name.
        2. Check hash subdirectory within cache directory for data file.
        3. If data file not in cache, download, compute hash, and verify hash.
        4. If hash correct, save file (+ append); otherwise, raise an error.

    Parameters
    ----------
    data_file : string
        name of file (not the full path)
    url : string
        URL for data file
    hashes : dictionary
        file names and md5 hashes (if empty, simply download file from url)
    cache_directory : string
        cache directory (full path)
    append : string
        append to output file (ex: '.nii.gz')
    verbose : bool
        print statements?

    Returns
    -------
    data_path : string
        data file name (full path)

    Examples
    --------
    >>> from mindboggle.mio.fetch_data import fetch_check_data
    >>> from mindboggle.mio.fetch_data import cache_hashes
    >>> # osf.io URL for OASIS-30_Atropos_template_to_MNI152_affine.txt
    >>> data_file = 'OASIS-30_Atropos_template_to_MNI152_affine.txt'
    >>> url = 'https://osf.io/ufydw/?action=download&version=1'
    >>> hashes = cache_hashes()
    >>> cache_directory = ''
    >>> append = ''
    >>> verbose = False
    >>> data_path = fetch_check_data(data_file, url, hashes, cache_directory,
    ...                              append, verbose) # doctest: +SKIP

    """
    import os
    import shutil

    from mindboggle.mio.fetch_data import fetch_data, fetch_hash

    # ------------------------------------------------------------------------
    # Set temporary cache directory if not specified:
    # ------------------------------------------------------------------------
    if not cache_directory:
        cache_directory = os.path.join(os.environ['HOME'], 'hash_temp')

    # ------------------------------------------------------------------------
    # Check hash table for file name, and store corresponding hash:
    # ------------------------------------------------------------------------
    if hashes and data_file in list(hashes):
        stored_hash = hashes[data_file]

        # --------------------------------------------------------------------
        # Create missing cache and hash directories:
        # --------------------------------------------------------------------
        if not os.path.exists(cache_directory):
            if verbose:
                print("Create missing cache directory: {0}".
                      format(cache_directory))
            os.mkdir(cache_directory)
        hash_dir = os.path.join(cache_directory, stored_hash)
        if not os.path.exists(hash_dir):
            if verbose:
                print("Create missing hash directory: {0}".format(hash_dir))
            os.mkdir(os.path.join(hash_dir))

        # --------------------------------------------------------------------
        # Check hash subdirectory for file:
        # --------------------------------------------------------------------
        data_path = os.path.join(hash_dir, data_file)
        if os.path.exists(data_path):
            if verbose:
                print("File already exists and matches hash: {0}".format(url))
            return data_path

        # --------------------------------------------------------------------
        # If file not in cache, download, compute hash, and verify:
        # --------------------------------------------------------------------
        else:
            if verbose:
                print("Retrieve file from URL: {0}".format(url))

            # Download file as a temporary file:
            temp_file = fetch_data(url)

            # Compute the file's hash:
            data_hash = fetch_hash(temp_file)

            # If hash matches name of the hash directory, save file:
            if os.path.join(cache_directory, data_hash) == hash_dir:
                # Add append:
                if append:
                    data_path += append
                if verbose:
                    print("Copy file to cache: {0}".format(data_path))
                shutil.copyfile(temp_file, data_path)
                return data_path
            else:
                raise IOError("Retrieved hash does not match stored hash.")
    else:
        raise IOError("Data file '{0}' not in hash table.".
                      format(data_file))


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
    use_ants_transforms : bool
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
    >>> segmented_file = 'ants/OASIS-TRT-20-1/tmpBrainSegmentation.nii.gz'
    >>> use_ants_transforms = True
    >>> m, s, a_s2t, w_s2t, a_t2s, w_t2s = fetch_ants_data(segmented_file,
    ...     use_ants_transforms) # doctest: +SKIP

    """
    import os

    prefix = segmented_file.split('BrainSegmentation.nii.gz', 1)[0]

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
    for ants_file in files:
        if not os.path.exists(ants_file):
            raise IOError('antsCorticalThickness.sh output ' + ants_file +
                          ' does not exist.')

    return mask, segments, affine_subject2template, warp_subject2template, \
                           affine_template2subject, warp_template2subject


# ============================================================================
# Doctests
# ============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)  # py.test --doctest-modules