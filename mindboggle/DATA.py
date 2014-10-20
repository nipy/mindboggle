#!/usr/bin/env python
"""
Functions related to Mindboggle data.


Authors:
    - Arno Klein, 2014  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2014,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def hashes_url():
    """
    Hashes and URL to verify retrieved data used by the Mindboggle software.

    Data include atlases and templates required for registration and labeling.
    See get_hash() in io_uri.py to create new hashes.

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
    url = 'http://mindboggle.info/data/cache/'
    cache_env = 'MINDBOGGLE_CACHE'
    cache = os.path.join(os.environ['HOME'], 'mindboggle_cache')

    #-------------------------------------------------------------------------
    # Atlases and templates:
    #-------------------------------------------------------------------------
    hashes['OASIS-30_Atropos_template.nii.gz'] = 'f95dbe37ab40e8ad59c1b1eabc7f230c'
    hashes['OASIS-30_Atropos_template_in_MNI152.nii.gz'] = '2efbe9c0755f4aa4506ad78e60ec319f'
    hashes['OASIS-30_Atropos_template_to_MNI152_affine.txt'] = 'f36e3d5d99f7c4a9bb70e2494ed7340b'
    hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_OASIS-30.nii.gz'] = '706dd2515b111b8eb9d7a7a0d76e833a'
    hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152.nii.gz'] = 'f0b48c3d717a27c022cd62f1f5c18e8f'
    hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_label_probabilities_in_OASIS-30.nii.gz'] = 'a27a66e8c3955ed2b3833f8dae9a637b'
    hashes['OASIS-TRT-20_jointfusion_DKT31_CMA_label_probabilities_in_MNI152.nii.gz'] = '529af0aa3254f2d496dc902cb0b93aa4'
    hashes['MNI152_T1_1mm_brain.nii.gz'] = '7c47f203858f57c50ee93bcf7b0662fc'
    hashes['lh.DKTatlas40.gcs'] = 'b28a600a713c269f4c561f66f64337b2'
    hashes['rh.DKTatlas40.gcs'] = '50150021c84b0e829aed3a68709404d7'
    hashes['depth_curv_border_nonborder_parameters.pkl'] = 'a943a0f46c2c2b3bfdfdf5a64176c6c3'

    return hashes, url, cache_env, cache
