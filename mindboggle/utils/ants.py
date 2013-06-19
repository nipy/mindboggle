#!/usr/bin/python
"""
Functions that call the legacy ANTs package created by PICSL at UPenn.

Authors:
Arno Klein, 2011-2013  .  arno@mindboggle.info  .  www.binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def register_volume(source, target, iterations='30x99x11', output_stem=''):
    """
    Use ANTs to register a source image volume to a target image volume.

    This program uses the ANTs SyN registration method.

    Parameters
    ----------
    source : string
        file name of source image volume
    target : string
        file name of target image volume
    iterations : string
        number of iterations ("0" for affine, "30x99x11" default)
    output_stem : string
        file name stem for output transform matrix

    Returns
    -------
    affine_transform : string
        file name for affine transform matrix
    nonlinear_transform : string
        file name for nonlinear transform nifti file

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.ants import register_volume
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> source = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> target = os.path.join(path, 'atlases', 'MNI152_T1_1mm_brain.nii.gz')
    >>> iterations = "0"
    >>> output_stem = ""
    >>> #
    >>> register_volume(source, target, iterations, output_stem)

    """
    import os

    if not output_stem:
        src = os.path.basename(source).split('.')[0]
        tgt = os.path.basename(target).split('.')[0]
        output_stem = os.path.join(os.getcwd(), src+'_to_'+tgt)

    args = ['ANTS 3 -m', 'CC[' + target + ',' + source + ',1,2]',
            '-r Gauss[2,0] -t SyN[0.5] -i', iterations,
            '-o', output_stem, '--use-Histogram-Matching',
            '--number-of-affine-iterations 10000x10000x10000x10000x10000']
    cmd = ' '.join(args)
    print(cmd); os.system(cmd); # p = Popen(args);

    affine_transform = output_stem + 'Affine.txt'
    nonlinear_transform = output_stem + 'Warp.nii.gz'

    return affine_transform, nonlinear_transform

"""
def transform_volume(source, target, interp='', output_stem=''):
    ""
    Use ANTs to transform a source image volume to a target image volume.

    This program uses the ANTs WarpImageMultiTransform function.

    Parameters
    ----------
    source : string
        file name of source image volume
    target : string
        file name of target image volume
    interp : string
        interpolation type ("NN" for nearest neighbor)
    output_stem : string
        file name stem for output transform matrix

    Returns
    -------
    output_stem : string
        file name stem for output transform matrix
"""