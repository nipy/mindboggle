#!/usr/bin/env python
"""
Functions related to Mindboggle data.


Authors:
    - Arno Klein, 2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def hashes_url():
    """
    Hashes and URL to verify retrieved data used by the mindboggle software.

    Data include atlases and templates necessary for registration and labeling.

    Returns
    -------
    hashes : dictionary
        dictionary of data file names and hashes for those files
    url : string
        URL to Mindboggle data
    cache_env : string
        cache environment variable
    cache : string
        path to cache directory

    """
    import os

    url = 'http://mindboggle.info/data/cache/'
    cache_env = 'MINDBOGGLE_CACHE'
    cache = os.path.join(os.environ['HOME'], 'mindboggle_cache')

    hashes = {}

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
    #-------------------------------------------------------------------------
    # Example Mindboggle output volumes and transforms:
    #-------------------------------------------------------------------------
    # hashes['mask.nii.gz'] = '294ae2739ea594ba84a8ca198a8d0885'
    # hashes['brain.nii.gz'] = '04bca82acadc65f192ed5f2e196d95fb'
    # hashes['brain_to_OASIS-TRT-20_template.nii.gz'] = '7fc71d10106da964f6974ac54bb13427'
    # hashes['brain_to_OASIS-TRT-20_templateAffine.txt'] = '321b09c4f2d440c10ca58c0243fdecde'
    # hashes['brain_to_OASIS-TRT-20_templateInverseWarp.nii.gz'] = '27f8ec6f8ed310217124673e5f69410d'
    # hashes['brain_to_OASIS-TRT-20_templateWarp.nii.gz'] = '85e819ee425be725238fd2e812fb4c62'
    #-------------------------------------------------------------------------
    # Example Mindboggle output labels:
    #-------------------------------------------------------------------------
    # hashes['lh.DKTatlas100.gcs.vtk'] = '6b209cd662fbe2ac886af9d5112bc1e8'
    # hashes['relabeled_lh.DKTatlas100.gcs.vtk'] = 'b01d363d66b922b631540430de09ddf0'
    # hashes['propagated_labels.nii.gz'] = '9255d39e8d3dbeef39fa7d13e483e788'
    # hashes['propagated_labels_and_brain_to_OASIS-TRT-20_template.nii.gz'] = '328d993cda53868358597d24db7b3d41'
    #-------------------------------------------------------------------------
    # Example Mindboggle output shapes:
    #-------------------------------------------------------------------------
    # hashes['lh.pial.area.vtk'] = '33141ab5de64d10464dcf24397df6143'
    # hashes['lh.pial.mean_curvature.vtk'] = 'c33d6c56863246f42329c9a85fe90137'
    # hashes['lh.pial.travel_depth.vtk'] = 'a583bcfc5054ff5bda30bc5492ac6783'
    # hashes['lh.pial.geodesic_depth.vtk'] = '0391224282c143d771aed3e9455e4ddc'
    # hashes['lh.sulc.vtk'] = '00f2b3d8bad993493d061e3d9608b1cf'
    # hashes['lh.thickness.vtk'] = 'a368898314c15046bb28bd89f52cf604'
    # hashes['likelihoods.vtk'] = 'a381c49b3b5f77dda027c42af0d93ec4'
    #-------------------------------------------------------------------------
    # Example Mindboggle output features:
    #-------------------------------------------------------------------------
    # hashes['folds.vtk'] = '2f1f3c9b20f735a203efc60221dc294f'
    # hashes['sulci.vtk'] = 'b70a240b5d1a4b91a546ffddd08cb41e'
    # hashes['fundi.vtk'] = '22704fa945b73c954693ab617d23a953'
    # hashes['smooth_skeletons.vtk'] = 'e8eb4751f3ee1eb47ac01453590aeb73'
    #-------------------------------------------------------------------------
    # Example Mindboggle output tables:
    #-------------------------------------------------------------------------
    # hashes['label_volumes.csv'] = '0b7aaf2d081515b7a268b87e3aee7602'
    # hashes['label_shapes.csv'] = 'a7652e0accd0c886b47b657bf36d54b1'
    # hashes['sulcus_shapes.csv'] = '0917d987b7d9fb9f81fbd8192e1ecdc2'
    # hashes['fundus_shapes.csv'] = '4329842930b178ba90540a9a4d308d5c'
    # hashes['vertex_shapes.csv'] = '09c6a9495cfff8671c4d4616f49fbf0b'

    return hashes, url, cache_env, cache
