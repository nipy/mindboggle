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
        write_shape_stats(), write_vertex_measures()

Authors:
Arno Klein, 2011-2014  .  arno@mindboggle.info  .  www.binarybottle.com

Copyright 2014,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def antsApplyTransformsToPoints(points, transform_files, inverse_booleans=[0]):
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
    >>> from mindboggle.thirdparty.ants import antsApplyTransformsToPoints
    >>> from mindboggle.io.vtk import read_vtk
    >>> transform_files = ['/Users/arno/Data/antsCorticalThickness/Twins-2-1/antsTemplateToSubject1GenericAffine.mat','/Users/arno/Data/antsCorticalThickness/Twins-2-1/antsTemplateToSubject0Warp.nii.gz','/Users/arno/Data/antsCorticalThickness/Twins-2-1/antsSubjectToTemplate0GenericAffine.mat','/Users/arno/Data/antsCorticalThickness/Twins-2-1/antsSubjectToTemplate1Warp.nii.gz']
    >>> transform_files = [transform_files[0],transform_files[1],'/Users/arno/Data/mindboggle_cache/f36e3d5d99f7c4a9bb70e2494ed7340b/OASIS-30_Atropos_template_to_MNI152_affine.txt']
    >>> vtk_file = '/Users/arno/mindboggle_working/Twins-2-1/Mindboggle/_hemi_lh/Surface_to_vtk/lh.pial.vtk'
    >>> faces, lines, indices, points, npoints, scalars, name, foo1 = read_vtk(vtk_file)
    >>> inverse_booleans = [0,0,1]
    >>> transformed_points = antsApplyTransformsToPoints(points, transform_files, inverse_booleans)

    """
    import os

    from mindboggle.guts.utilities import execute

    #-------------------------------------------------------------------------
    # Write points (x,y,z,1) to a .csv file:
    #-------------------------------------------------------------------------
    points_file = os.path.join(os.getcwd(), 'points.csv')
    fid = open(points_file, 'wa')
    fid.write('x,y,z,t\n')
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
    execute(cmd, 'os')
    if not os.path.exists(transformed_points_file):
        str1 = "antsApplyTransformsToPoints did not create "
        raise(IOError(str1 + transformed_points_file + "."))

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
    >>> import os
    >>> from mindboggle.thirdparty.ants import ImageMath
    >>> from mindboggle.io.plot import plot_volumes
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> volume1 = os.path.join(path, 'arno', 'mri', 't1weighted.nii.gz')
    >>> volume2 = os.path.join(path, 'arno', 'mri', 'mask.nii.gz')
    >>> operator = 'm'
    >>> output_file = ''
    >>> output_file = ImageMath(volume1, volume2, operator, output_file)
    >>> # View
    >>> plot_volumes(output_file)

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
        raise(IOError("ImageMath did not create " + output_file + "."))

    return output_file


def ThresholdImage(volume, output_file='', threshlo=1, threshhi=10000):
    """
    Use the ThresholdImage function in ANTs to threshold image volume::

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
    >>> from mindboggle.io.plot import plot_volumes
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> volume = os.path.join(path, 'arno', 'mri', 't1weighted.nii.gz')
    >>> output_file = ''
    >>> threshlo = 500
    >>> threshhi = 10000
    >>> output_file = ThresholdImage(volume, output_file, threshlo, threshhi)
    >>> # View
    >>> plot_volumes(output_file)

    """
    import os
    from mindboggle.guts.utilities import execute

    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   'threshold_' + os.path.basename(volume))

    cmd = 'ThresholdImage 3 {0} {1} {2} {3}'.format(volume, output_file,
                                                    threshlo, threshhi)
    execute(cmd, 'os')
    if not os.path.exists(output_file):
        raise(IOError("ThresholdImage did not create " + output_file + "."))

    return output_file


def PropagateLabelsThroughMask(mask, labels, mask_index=None,
                               output_file='', binarize=True, stopvalue=''):
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
    output_file : string
        nibabel-readable labeled image volume
    binarize : Boolean
        binarize mask?
    stopvalue : integer
        stopping value

    Returns
    -------
    output_file : string
        name of labeled output nibabel-readable image volume

    Examples
    --------
    >>> import os
    >>> from mindboggle.thirdparty.ants import PropagateLabelsThroughMask
    >>> from mindboggle.io.plot import plot_volumes
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> labels = os.path.join(path, 'arno', 'labels', 'labels.DKT25.manual.nii.gz')
    >>> mask = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> mask_index = None
    >>> output_file = ''
    >>> binarize = True
    >>> stopvalue = None
    >>> output_file = PropagateLabelsThroughMask(mask, labels, mask_index, output_file, binarize, stopvalue)
    >>> # View
    >>> plot_volumes(output_file)

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
        cmd = 'ThresholdImage 3 {0} {1} {2} {3} 1 0'.format(mask, mask2,
               mask_index, mask_index)
        execute(cmd)
    else:
        mask2 = mask

    # Propagate labels:
    cmd = ['ImageMath', '3', output_file, 'PropagateLabelsThroughMask',
            mask2, labels]
    if stopvalue:
        cmd.extend(stopvalue)
    execute(cmd, 'os')
    if not os.path.exists(output_file):
        raise(IOError("ImageMath did not create " + output_file + "."))

    return output_file


def fill_volume_with_surface_labels(hemi, left_mask, right_mask,
                                    surface_files, mask_index=None,
                                    output_file='', binarize=False):
    """
    Use ANTs to fill a volume mask with surface mesh labels.

    Note ::

        - This uses PropagateLabelsThroughMask in the ANTs ImageMath function.

        - Partial volume information is lost when mapping surface to volume.

    Parameters
    ----------
    hemi : string
        either 'lh' or 'rh', indicating brain's left or right hemisphere
    left_mask : string
        nibabel-readable left brain image mask volume
    right_mask : string
        nibabel-readable left brain image mask volume
    surface_files : string or list of strings
        VTK file(s) containing surface mesh(es) with labels as scalars
    mask_index : integer (optional)
        mask with just voxels having this value
    output_file : string
        name of output file
    binarize : Boolean
        binarize mask?

    Returns
    -------
    output_file : string
        name of labeled output nibabel-readable image volume

    Examples
    --------
    >>> import os
    >>> from mindboggle.thirdparty.ants import fill_volume_with_surface_labels
    >>> from mindboggle.io.plot import plot_volumes
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> surface_files = [os.path.join(path, 'arno', 'labels',
    >>>     'lh.labels.DKT25.manual.vtk'), os.path.join(path, 'arno', 'labels',
    >>>     'rh.labels.DKT25.manual.vtk')]
    >>> # For a quick test, simply mask with whole brain:
    >>> hemi = 'rh'
    >>> left_mask = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> right_mask = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> mask_index = None
    >>> output_file = ''
    >>> binarize = True
    >>> output_file = fill_volume_with_surface_labels(hemi, left_mask, right_mask, surface_files, mask_index, output_file, binarize)
    >>> # View
    >>> plot_volumes(output_file)

    """
    import os

    from mindboggle.io.vtk import transform_to_volume
    from mindboggle.guts.relabel import overwrite_volume_labels
    from mindboggle.thirdparty.ants import PropagateLabelsThroughMask

    if isinstance(surface_files, str):
        surface_files = [surface_files]

    if hemi == 'lh':
        mask = left_mask
    elif hemi == 'rh':
        mask = right_mask

    # Transform vtk coordinates to voxel index coordinates in a target
    # volume by using the header transformation:
    surface_in_volume = transform_to_volume(surface_files[0], mask)

    # Do the same for additional vtk surfaces:
    if len(surface_files) == 2:
        surfaces_in_volume = os.path.join(os.getcwd(), 'surfaces.nii.gz')
        surface_in_volume2 = transform_to_volume(surface_files[1], mask)

        overwrite_volume_labels(surface_in_volume, surface_in_volume2,
                                surfaces_in_volume, ignore_labels=[0])
        surface_in_volume = surfaces_in_volume

    # Use ANTs to fill a binary volume mask with initial labels:
    output_file = PropagateLabelsThroughMask(mask, surface_in_volume,
                                             mask_index, output_file,
                                             binarize)
    if not os.path.exists(output_file):
        str1 = "PropagateLabelsThroughMask() did not create "
        raise(IOError(str1 + output_file + "."))

    return output_file  # surface_in_volume


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
#     >>> source = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
#     >>> target = os.path.join(path, 'atlases', 'MNI152_T1_1mm_brain.nii.gz')
#     >>> iterations = "0"
#     >>> output_stem = ""
#     >>> #
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
#         raise(IOError("ANTs did not create " + affine_transform + "."))
#     if not os.path.exists(nonlinear_transform):
#         raise(IOError("ANTs did not create " + nonlinear_transform + "."))
#     if not os.path.exists(nonlinear_inverse_transform):
#         raise(IOError("ANTs did not create " + nonlinear_inverse_transform + "."))
#
#     return affine_transform, nonlinear_transform,\
#            nonlinear_inverse_transform, output_stem
#
#
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
#     inverse : Boolean
#         apply inverse transform?
#     affine_only : Boolean
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
#         raise(IOError("WarpImageMultiTransform did not create " + output + "."))
#
#     return output
#
#
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
#     >>> from mindboggle.thirdparty.ants import ComposeMultiTransform
#     >>> transform_files = ['affine1.mat', 'affine2.mat']
#     >>> transform_files = ['/data/Brains/Mindboggle101/antsCorticalThickness/OASIS-TRT-20_volumes/OASIS-TRT-20-1/antsTemplateToSubject0GenericAffine.mat','/data/Brains/Mindboggle101/antsCorticalThickness/OASIS-TRT-20_volumes/OASIS-TRT-20-1/antsTemplateToSubject0GenericAffine.mat']
#     >>> inverse_Booleans = [False, False]
#     >>> output_transform_file = ''
#     >>> ext = '.txt'
#     >>> ComposeMultiTransform(transform_files, inverse_Booleans, output_transform_file, ext)
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
#     cmd = ['ComposeMultiTransform 3', output_transform_file, ' '.join(xfms)]
#     print(cmd)
#     execute(cmd, 'os')
#     #if not os.path.exists(output_transform_file):
#     #    raise(IOError(output_transform_file + " not found"))
#
#     return output_transform_file


