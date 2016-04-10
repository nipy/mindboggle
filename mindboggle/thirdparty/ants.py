#!/usr/bin/python
"""
Functions that call ANTs (UPenn's PICSL group) commands.

Mindboggle functions call the following ANTs functions ::

    ImageMath:
        'PropagateLabelsThroughMask' option in PropagateLabelsThroughMask()
        if modify_surface_labels set to True:
            PropagateLabelsThroughMask() called and '+' option in mindboggle

    ThresholdImage:
        PropagateLabelsThroughMask()

    antsApplyTransformsToPoints:
        write_shape_stats(), write_vertex_measures() in mindboggle

Authors:
    - Arno Klein, 2012-2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def antsApplyTransformsToPoints(points, transform_files,
                                inverse_booleans=[0]):
    """
    Run ANTs antsApplyTransformsToPoints function to transform points.
    (Creates pre- and post-transformed .csv points files for ANTs.)

    Parameters
    ----------
    points : list of lists of three integers
        point coordinate data
    transform_files : list
        transform file names
    inverse_booleans : list
        for each transform, one to apply inverse of transform (otherwise zero)

    Returns
    -------
    transformed_points : list of lists of three integers
        transformed point coordinate data

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.thirdparty.ants import antsApplyTransformsToPoints
    >>> from mindboggle.mio.vtks import read_points
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> xfm1 = fetch_data(urls['ants_affine_template2subject'])
    >>> xfm2 = fetch_data(urls['ants_warp_template2subject'])
    >>> xfm3 = fetch_data(urls['OASIS-30_Atropos_template_to_MNI152_affine'])
    >>> transform_files = [xfm1, xfm2, xfm3]
    >>> vtk_file = fetch_data(urls['left_pial'])
    >>> points  = read_points(vtk_file)
    >>> inverse_booleans = [0,0,1]
    >>> transformed_points = antsApplyTransformsToPoints(points,
    ...     transform_files, inverse_booleans) # doctest: +SKIP
    >>> print(np.array_str(np.array(transformed_points[0:5]),
    ...       precision=5, suppress_small=True)) # doctest: +SKIP
    [[-11.23189 -46.78223 -39.88869]
     [-11.71384 -46.87075 -40.13328]
     [-12.56237 -46.99126 -40.04564]
     [ -9.66693 -46.0446  -41.36334]
     [-10.67998 -46.45458 -40.7572 ]]

    """
    import os
    from io import open

    from mindboggle.guts.utilities import execute

    #-------------------------------------------------------------------------
    # Write points (x,y,z,1) to a .csv file:
    #-------------------------------------------------------------------------
    points_file = os.path.join(os.getcwd(), 'points.csv')
    fid = open(points_file, 'w')
    fid.write('x,y,z,t\n')
    fid.close()
    fid = open(points_file, 'a')
    for point in points:
        string_of_zeros = (4 - len(point)) * ',0'
        fid.write(','.join([str(x) for x in point]) + string_of_zeros + '\n')
    fid.close()

    #-------------------------------------------------------------------------
    # Apply transforms to points in .csv file:
    #-------------------------------------------------------------------------
    transformed_points_file = os.path.join(os.getcwd(),
                                           'transformed_points.csv')
    transform_string = ''
    for ixfm, transform_file in enumerate(transform_files):
        transform_string += " --t [{0},{1}]".\
            format(transform_file, str(inverse_booleans[ixfm]))
    cmd = ['antsApplyTransformsToPoints', '-d', '3', '-i', points_file,
           '-o', transformed_points_file, transform_string]
    try:
        execute(cmd, 'os')
    except:
        raise Exception("Cannot find antsApplyTransformsToPoints command.")

    if not os.path.exists(transformed_points_file):
        raise IOError("antsApplyTransformsToPoints did not create {0}.".
                      format(transformed_points_file))

    #-------------------------------------------------------------------------
    # Return transformed points:
    #-------------------------------------------------------------------------
    fid = open(transformed_points_file, 'r')
    lines = fid.readlines()
    fid.close()
    transformed_points = []
    for iline, line in enumerate(lines):
        if iline > 0:
            point_xyz1 = [float(x) for x in line.split(',')]
            transformed_points.append(point_xyz1[0:3])

    return transformed_points


def ImageMath(volume1, volume2, operator='m', output_file=''):
    """
    Use the ImageMath function in ANTs to perform operation on two volumes::

        m         : Multiply ---  use vm for vector multiply
        +         : Add ---  use v+ for vector add
        -         : Subtract ---  use v- for vector subtract
        /         : Divide
        ^         : Power
        exp       : Take exponent exp(imagevalue*value)
        addtozero : add image-b to image-a only over points where image-a has zero values
        overadd   : replace image-a pixel with image-b pixel if image-b pixel is non-zero
        abs       : absolute value
        total     : Sums up values in an image or in image1*image2 (img2 is the probability mask)
        mean      :  Average of values in an image or in image1*image2 (img2 is the probability mask)
        vtotal    : Sums up volumetrically weighted values in an image or in image1*image2 (img2 is the probability mask)
        Decision  : Computes result=1./(1.+exp(-1.0*( pix1-0.25)/pix2))
        Neg       : Produce image negative


    Parameters
    ----------
    volume1 : string
        nibabel-readable image volume
    volume2 : string
        nibabel-readable image volume
    operator : string
        ImageMath string corresponding to mathematical operator
    output_file : string
        nibabel-readable image volume

    Returns
    -------
    output_file : string
        name of output nibabel-readable image volume

    Examples
    --------
    >>> # Mask head with brain mask:
    >>> import os
    >>> from mindboggle.thirdparty.ants import ImageMath
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> volume1 = fetch_data(urls['T1_001'])
    >>> volume2 = fetch_data(urls['ants_mask'])
    >>> operator = 'm'
    >>> output_file = ''
    >>> output_file = ImageMath(volume1, volume2, operator, output_file) # doctest: +SKIP

    View result (skip test):

    >>> from mindboggle.mio.plots import plot_volumes
    >>> plot_volumes(output_file) # doctest: +SKIP

    """
    import os
    from mindboggle.guts.utilities import execute

    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   os.path.basename(volume1) + '_' +
                                   os.path.basename(volume2))
    cmd = ['ImageMath', '3', output_file, operator, volume1, volume2]
    execute(cmd, 'os')
    if not os.path.exists(output_file):
        raise IOError("ImageMath did not create " + output_file + ".")

    return output_file


def ThresholdImage(volume, output_file='', threshlo=1, threshhi=10000):
    """
    Use the ThresholdImage function in ANTs to threshold image volume.

    Usage: ThresholdImage ImageDimension ImageIn.ext outImage.ext
           threshlo threshhi <insideValue> <outsideValue>

    Parameters
    ----------
    volume : string
        nibabel-readable image volume
    output_file : string
        nibabel-readable image volume
    threshlo : integer
        lower threshold
    threshhi : integer
        upper threshold

    Returns
    -------
    output_file : string
        name of output nibabel-readable image volume

    Examples
    --------
    >>> import os
    >>> from mindboggle.thirdparty.ants import ThresholdImage
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> volume = fetch_data(urls['T1_001'])
    >>> os.rename(volume, volume + '.nii.gz')
    >>> volume += '.nii.gz'
    >>> output_file = ''
    >>> threshlo = 500
    >>> threshhi = 10000
    >>> output_file = ThresholdImage(volume, output_file, threshlo, threshhi) # doctest: +SKIP

    View result (skip test):

    >>> from mindboggle.mio.plots import plot_volumes
    >>> plot_volumes(output_file) # doctest: +SKIP

    """
    import os
    from mindboggle.guts.utilities import execute

    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   'threshold_' + os.path.basename(volume))
    cmd = ['ThresholdImage', '3', volume, output_file,
           str(threshlo), str(threshhi)]
    execute(cmd, 'os')
    if not os.path.exists(output_file):
        raise IOError("ThresholdImage did not create " + output_file + ".")

    return output_file


def PropagateLabelsThroughMask(mask, labels, mask_index=None, output_file='',
                               binarize=True, stopvalue=''):
    """
    Use ANTs to fill a binary volume mask with initial labels.

    This program uses ThresholdImage and the ImageMath
    PropagateLabelsThroughMask functions in ANTs.

    ThresholdImage ImageDimension ImageIn.ext outImage.ext
        threshlo threshhi <insideValue> <outsideValue>

    PropagateLabelsThroughMask: Final output is the propagated label image.
        ImageMath ImageDimension Out.ext PropagateLabelsThroughMask
        speed/binaryimagemask.nii.gz initiallabelimage.nii.gz ...

    Parameters
    ----------
    mask : string
        nibabel-readable image volume
    labels : string
        nibabel-readable image volume with integer labels
    mask_index : integer (optional)
        mask with just voxels having this value
    output_file : string (optional)
        nibabel-readable labeled image volume
    binarize : bool (optional)
        binarize mask?
    stopvalue : integer (optional)
        stopping value

    Returns
    -------
    output_file : string
        name of labeled output nibabel-readable image volume

    Examples
    --------
    >>> # Propagate FreeSurfer labels through brain mask:
    >>> import os
    >>> from mindboggle.thirdparty.ants import PropagateLabelsThroughMask
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> labels = fetch_data(urls['freesurfer_labels'])
    >>> mask = fetch_data(urls['ants_mask'])
    >>> mask_index = None
    >>> output_file = ''
    >>> binarize = True
    >>> stopvalue = None
    >>> output_file = PropagateLabelsThroughMask(mask, labels, mask_index,
    ...     output_file, binarize, stopvalue) # doctest: +SKIP

    View result (skip test):

    >>> from mindboggle.mio.plots import plot_volumes
    >>> plot_volumes(output_file) # doctest: +SKIP

    """
    import os
    from mindboggle.guts.utilities import execute

    if not output_file:
        #output_file = os.path.join(os.getcwd(),
        #                           'PropagateLabelsThroughMask.nii.gz')
        output_file = os.path.join(os.getcwd(),
                                   os.path.basename(labels) + '_through_' +
                                   os.path.basename(mask))

    print('mask: {0}, labels: {1}'.format(mask, labels))

    # Binarize image volume:
    if binarize:
        temp_file = os.path.join(os.getcwd(),
                                 'PropagateLabelsThroughMask.nii.gz')
        cmd = ['ThresholdImage', '3', mask, temp_file, '0 1 0 1']
        execute(cmd, 'os')
        mask = temp_file

    # Mask with just voxels having mask_index value:
    if mask_index:
        mask2 = os.path.join(os.getcwd(), 'temp.nii.gz')

        cmd = ['ThresholdImage', '3', mask, mask2,
               str(mask_index), str(mask_index)]
        execute(cmd, 'os')
    else:
        mask2 = mask

    # Propagate labels:

    cmd = ['ImageMath', '3', output_file, 'PropagateLabelsThroughMask',
           mask2, labels]
    if stopvalue:
        cmd.extend(str(stopvalue))
    execute(cmd, 'os')
    if not os.path.exists(output_file):
        raise IOError("ImageMath did not create " + output_file + ".")

    return output_file


def ResampleImageBySpacing(volume, output_file='', outxspc=1, outyspc=1,
                           outzspc=1, dosmooth=0, addvox=0, nninterp=1):
    """
    Use the ResampleImageBySpacing function in ANTs to resample image volume.

    Usage: ResampleImageBySpacing ImageDimension ImageIn.ext outImage.ex
           outxspc outyspc <outzspc> <dosmooth?> <addvox> <nninterp?>

    Parameters
    ----------
    volume : string
        nibabel-readable image volume
    output_file : string
        nibabel-readable image volume
    outxspc : integer
        output x-spacing
    outyspc : integer
        output y-spacing
    outzspc : integer
        output z-spacing
    dosmooth : bool
        smooth?
    addvox : integer
        pad each dimension by addvox
    nninterp : bool
        nearest-neighbor interpolation?

    Returns
    -------
    output_file : string
        name of output nibabel-readable image volume

    Examples
    --------
    >>> # Resample image so that 1mm voxels are 0.5mm voxels:
    >>> import os
    >>> from mindboggle.thirdparty.ants import ResampleImageBySpacing
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> volume = fetch_data(urls['T1_001'])
    >>> os.rename(volume, volume + '.nii.gz')
    >>> volume += '.nii.gz'
    >>> output_file = ''
    >>> outxspc = 1/2.0
    >>> outyspc = 1/2.0
    >>> outzspc = 1/2.0
    >>> dosmooth = 0
    >>> addvox = 0
    >>> nninterp = 1
    >>> output_file = ResampleImageBySpacing(volume, output_file, outxspc,
    ...     outxspc, outzspc, dosmooth, addvox, nninterp) # doctest: +SKIP

    View result (skip test):

    >>> from mindboggle.mio.plots import plot_volumes
    >>> plot_volumes(output_file) # doctest: +SKIP

    """
    import os
    from mindboggle.guts.utilities import execute

    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   'resampled_' + os.path.basename(volume))
    cmd = ['ResampleImageBySpacing', '3', volume, output_file,
           str(outxspc), str(outyspc), str(outzspc)]
    execute(cmd, 'os')
    if not os.path.exists(output_file):
        raise IOError("ResampleImageBySpacing did not create {0).".
                      format(output_file))

    return output_file


# def fill_volume_with_surface_labels(hemi, left_mask, right_mask,
#                                     surface_files, mask_index=None,
#                                     output_file='', binarize=False):
#     """
#     Use ANTs to fill a volume mask with surface mesh labels.
#
#     Note ::
#
#         - This uses PropagateLabelsThroughMask in the ANTs ImageMath function.
#
#         - Partial volume information is lost when mapping surface to volume.
#
#         - Check to make sure transform_to_volume() is operational.
#
#     Parameters
#     ----------
#     hemi : string
#         either 'lh' or 'rh', indicating brain's left or right hemisphere
#     left_mask : string
#         nibabel-readable left brain image mask volume
#     right_mask : string
#         nibabel-readable left brain image mask volume
#     surface_files : string or list of strings
#         VTK file(s) containing surface mesh(es) with labels as scalars
#     mask_index : integer (optional)
#         mask with just voxels having this value
#     output_file : string
#         name of output file
#     binarize : bool
#         binarize mask?
#
#     Returns
#     -------
#     output_file : string
#         name of labeled output nibabel-readable image volume
#
#     Examples
#     --------
#     >>> from mindboggle.thirdparty.ants import fill_volume_with_surface_labels
#     >>> from mindboggle.mio.fetch_data import prep_tests
#     >>> urls, fetch_data = prep_tests()
#     >>> labels_left = fetch_data(urls['left_freesurfer_labels'])
#     >>> labels_right = fetch_data(urls['right_freesurfer_labels'])
#     >>> surface_files = [labels_left, labels_left]
#     >>> left_mask = fetch_data(urls['T1_001'])
#     >>> right_mask = fetch_data(urls['T1_001'])
#     >>> # For a quick test, simply mask with whole brain:
#     >>> hemi = 'rh'
#     >>> mask_index = None
#     >>> output_file = ''
#     >>> binarize = True
#     >>> output_file = fill_volume_with_surface_labels(hemi, left_mask,
#     ...     right_mask, surface_files, mask_index, output_file, binarize)
#
#     >>> # View
#     >>> from mindboggle.mio.plots import plot_volumes
#     >>> plot_volumes(output_file) # doctest: +SKIP
#
#     """
#     import os
#
#     from mindboggle.mio.vtks import transform_to_volume
#     from mindboggle.guts.relabel import overwrite_volume_labels
#     from mindboggle.thirdparty.ants import PropagateLabelsThroughMask
#
#     if isinstance(surface_files, str):
#         surface_files = [surface_files]
#
#     if hemi == 'lh':
#         mask = left_mask
#     elif hemi == 'rh':
#         mask = right_mask
#
#     # Transform vtk coordinates to voxel index coordinates in a target
#     # volume by using the header transformation:
#     surface_in_volume = transform_to_volume(surface_files[0], mask)
#
#     # Do the same for additional vtk surfaces:
#     if len(surface_files) == 2:
#         surfaces_in_volume = os.path.join(os.getcwd(), 'surfaces.nii.gz')
#         surface_in_volume2 = transform_to_volume(surface_files[1], mask)
#
#         overwrite_volume_labels(surface_in_volume, surface_in_volume2,
#                                 surfaces_in_volume, ignore_labels=[0])
#         surface_in_volume = surfaces_in_volume
#
#     # Use ANTs to fill a binary volume mask with initial labels:
#     output_file = PropagateLabelsThroughMask(mask, surface_in_volume,
#                                              mask_index, output_file,
#                                              binarize)
#     if not os.path.exists(output_file):
#         raise IOError("PropagateLabelsThroughMask() did not create {0}.".
#                       format(output_file))
#
#     return output_file  # surface_in_volume


# def ANTS(source, target, iterations='30x99x11', output_stem=''):
#     """
#     Use ANTs to register a source image volume to a target image volume.
#
#     This program uses the ANTs SyN registration method.
#
#     Parameters
#     ----------
#     source : string
#         file name of source image volume
#     target : string
#         file name of target image volume
#     iterations : string
#         number of iterations ("0" for affine, "30x99x11" default)
#     output_stem : string
#         file name stem for output transform matrix
#
#     Returns
#     -------
#     affine_transform : string
#         file name for affine transform matrix
#     nonlinear_transform : string
#         file name for nonlinear transform nifti file
#     nonlinear_inverse_transform : string
#         file name for nonlinear inverse transform nifti file
#     output_stem : string
#         file name stem for output transform matrix
#
#     Examples
#     --------
#     >>> import os
#     >>> from mindboggle.thirdparty.ants import ANTS
#     >>> path = os.environ['MINDBOGGLE_DATA']
#     >>> source = os.path.join(path, 'mri', 't1weighted_brain.nii.gz')
#     >>> target = os.path.join(path, 'atlases', 'MNI152_T1_1mm_brain.nii.gz')
#     >>> iterations = "0"
#     >>> output_stem = ""
#     >>> ANTS(source, target, iterations, output_stem)
#
#     """
#     import os
#     from mindboggle.guts.utilities import execute
#
#     if not output_stem:
#         src = os.path.basename(source).split('.')[0]
#         tgt = os.path.basename(target).split('.')[0]
#         output_stem = os.path.join(os.getcwd(), src+'_to_'+tgt)
#
#     cmd = ['ANTS', '3', '-m CC[' + target + ',' + source + ',1,2]',
#             '-r Gauss[2,0]', '-t SyN[0.5] -i', iterations,
#             '-o', output_stem, '--use-Histogram-Matching',
#             '--number-of-affine-iterations 10000x10000x10000x10000x10000']
#     execute(cmd, 'os')
#
#     affine_transform = output_stem + 'Affine.txt'
#     nonlinear_transform = output_stem + 'Warp.nii.gz'
#     nonlinear_inverse_transform = output_stem + 'InverseWarp.nii.gz'
#
#     if not os.path.exists(affine_transform):
#         raise IOError("ANTs did not create " + affine_transform + ".")
#     if not os.path.exists(nonlinear_transform):
#         raise IOError("ANTs did not create " + nonlinear_transform + ".")
#     if not os.path.exists(nonlinear_inverse_transform):
#         raise IOError("ANTs did not create " + nonlinear_inverse_transform + ".")
#
#     return affine_transform, nonlinear_transform,\
#            nonlinear_inverse_transform, output_stem


# def WarpImageMultiTransform(source, target, output='',
#                             interp='--use-NN', xfm_stem='',
#                             affine_transform='', nonlinear_transform='',
#                             inverse=False, affine_only=False):
#     """
#     Use ANTs to transform a source image volume to a target image volume.
#
#     This program uses the ANTs WarpImageMultiTransform function.
#
#     Parameters
#     ----------
#     source : string
#         file name of source image volume
#     target : string
#         file name of target (reference) image volume
#     output : string
#         file name of output image volume
#     interp : string
#         interpolation type ("--use-NN" for nearest neighbor)
#     xfm_stem : string
#         file name stem for output transform
#     affine_transform : string
#         file containing affine transform
#     nonlinear_transform : string
#         file containing nonlinear transform
#     inverse : bool
#         apply inverse transform?
#     affine_only : bool
#         apply only affine transform?
#
#     Returns
#     -------
#     output : string
#         output label file name
#
#     """
#     import os
#     import sys
#     from mindboggle.guts.utilities import execute
#
#     if xfm_stem:
#         affine_transform = xfm_stem + 'Affine.txt'
#         if inverse:
#             nonlinear_transform = xfm_stem + 'InverseWarp.nii.gz'
#         else:
#             nonlinear_transform = xfm_stem + 'Warp.nii.gz'
#     elif not affine_transform and not nonlinear_transform:
#         sys.exit('Provide either xfm_stem or affine_transform and '
#                  'nonlinear_transform.')
#
#     if not output:
#         output = os.path.join(os.getcwd(), 'WarpImageMultiTransform.nii.gz')
#
#     if not os.path.exists(nonlinear_transform):
#         affine_only = True
#
#     if affine_only:
#         if inverse:
#             cmd = ['WarpImageMultiTransform', '3', source, output, '-R',
#                    target, interp, '-i', affine_transform]
#         else:
#             cmd = ['WarpImageMultiTransform', '3', source, output, '-R',
#                    target, interp, affine_transform]
#     else:
#         if inverse:
#             cmd = ['WarpImageMultiTransform', '3', source, output, '-R',
#                    target, interp, '-i', affine_transform, nonlinear_transform]
#         else:
#             cmd = ['WarpImageMultiTransform', '3', source, output, '-R',
#                    target, interp, nonlinear_transform, affine_transform]
#     execute(cmd, 'os')
#
#     if not os.path.exists(output):
#         raise IOError("WarpImageMultiTransform did not create " + output + ".")
#
#     return output


# def ComposeMultiTransform(transform_files, inverse_Booleans,
#                           output_transform_file='', ext='.txt'):
#     """
#     Run ANTs ComposeMultiTransform function to create a single transform.
#
#     Parameters
#     ----------
#     transform_files : list of strings
#         transform file names
#     inverse_Booleans : list of Booleans
#         Booleans to indicate which transforms to take the inverse of
#     output_transform_file : string
#         transform file name
#     ext : string
#         '.txt' to save transform file as text, '.mat' for data file
#
#     Returns
#     -------
#     output_transform_file : string
#         single composed transform file name
#
#     Examples
#     --------
#     >>> import os
#     >>> from mindboggle.thirdparty.ants import ComposeMultiTransform
#     >>> from mindboggle.mio.fetch_data import prep_tests
#     >>> urls, fetch_data = prep_tests()
#     >>> xfm1 = fetch_data(urls['ants_affine_template2subject'])
#     >>> xfm2 = fetch_data(urls['ants_affine_subject2template'])
#     >>> transform_files = [xfm1, xfm2]
#     >>> inverse_Booleans = [False, False]
#     >>> output_transform_file = ''
#     >>> ext = '.txt'
#     >>> ComposeMultiTransform(transform_files, inverse_Booleans,
#     ...                       output_transform_file, ext)
#
#     """
#     import os
#
#     from mindboggle.guts.utilities import execute
#
#     if not output_transform_file:
#         output_transform_file = os.path.join(os.getcwd(), 'affine' + ext)
#
#     xfms = []
#     for ixfm, xfm in enumerate(transform_files):
#         if inverse_Booleans[ixfm]:
#             xfms.append('-i')
#         xfms.append(xfm)
#
#     cmd = ['ComposeMultiTransform', '3', output_transform_file, ' '.join(xfms)]
#     execute(cmd, 'os')
#     #if not os.path.exists(output_transform_file):
#     #    raise IOError(output_transform_file + " not found")
#
#     return output_transform_file

#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod()
